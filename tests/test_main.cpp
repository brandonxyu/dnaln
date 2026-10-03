// Unit and integration tests (no external framework). Run: make test
#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

#include <unistd.h>

#include "align.hpp"
#include "chain.hpp"
#include "common.hpp"
#include "evaluate.hpp"
#include "index.hpp"
#include "mapper.hpp"
#include "minimizer.hpp"
#include "reference.hpp"
#include "seqio.hpp"
#include "simulate.hpp"

using namespace dnaln;

// ---------------------------------------------------------------- mini framework
namespace {
int g_checks = 0, g_failures = 0;
struct TestCase {
  const char* name;
  void (*fn)();
};
std::vector<TestCase>& registry() {
  static std::vector<TestCase> r;
  return r;
}
struct Registrar {
  Registrar(const char* n, void (*f)()) { registry().push_back({n, f}); }
};
void report_failure(const char* expr, const char* file, int line, const std::string& detail) {
  if (++g_failures <= 25) std::fprintf(stderr, "    FAILED %s at %s:%d %s\n", expr, file, line, detail.c_str());
}
}  // namespace

#define TEST(name)                                \
  static void name();                             \
  static Registrar registrar_##name(#name, name); \
  static void name()
#define CHECK(cond)                                                   \
  do {                                                                \
    ++g_checks;                                                       \
    if (!(cond)) report_failure(#cond, __FILE__, __LINE__, "");       \
  } while (0)
#define CHECK_EQ(a, b)                                                                                    \
  do {                                                                                                    \
    ++g_checks;                                                                                           \
    const auto va_ = (a);                                                                                 \
    const auto vb_ = (b);                                                                                 \
    if (!(va_ == vb_))                                                                                    \
      report_failure(#a " == " #b, __FILE__, __LINE__,                                                    \
                     "(" + std::to_string((long long)va_) + " vs " + std::to_string((long long)vb_) + ")"); \
  } while (0)

// ---------------------------------------------------------------- helpers
namespace {

std::vector<uint8_t> random_codes(Rng& rng, size_t n, double n_rate = 0.0) {
  std::vector<uint8_t> v(n);
  for (auto& c : v) c = rng.uniform() < n_rate ? 4 : static_cast<uint8_t>(rng.below(4));
  return v;
}

// Mutates a sequence with substitutions and short indels.
std::vector<uint8_t> mutate(Rng& rng, const std::vector<uint8_t>& s, double sub, double indel) {
  std::vector<uint8_t> o;
  for (size_t i = 0; i < s.size(); ++i) {
    if (rng.uniform() < indel) {
      if (rng.next() & 1) continue;                                  // deletion
      o.push_back(static_cast<uint8_t>(rng.below(4)));               // insertion
    }
    o.push_back(rng.uniform() < sub ? static_cast<uint8_t>(rng.below(4)) : s[i]);
  }
  return o;
}

// Straightforward full-matrix DP with the same band semantics (reference for the kernels).
int naive_banded_sw(const std::vector<uint8_t>& q, const std::vector<uint8_t>& t_in, const ScoreParams& sp) {
  const int ql = static_cast<int>(q.size()), W = sp.band_cols(), tl = ql + W - 1;
  std::vector<uint8_t> t(static_cast<size_t>(std::max(tl, 0)), 4);
  std::copy(t_in.begin(), t_in.begin() + std::min<long>(static_cast<long>(t_in.size()), tl), t.begin());
  const int NEG = -1000000, oe = sp.gap_open + sp.gap_ext, ge = sp.gap_ext;
  std::vector<std::vector<int>> H(ql + 1, std::vector<int>(tl + 1, 0)), E(H), F(H);
  for (auto* m : {&E, &F})
    for (auto& row : *m) std::fill(row.begin(), row.end(), NEG);
  int best = 0;
  for (int i = 1; i <= ql; ++i) {
    for (int j = 1; j <= tl; ++j) {
      const int d = j - i;
      if (d < 0 || d >= W) continue;  // out of band: H = 0, E = F = -inf
      E[i][j] = std::max(H[i - 1][j] - oe, E[i - 1][j] - ge);
      F[i][j] = std::max(H[i][j - 1] - oe, F[i][j - 1] - ge);
      H[i][j] = std::max({0, H[i - 1][j - 1] + sp.sub(q[i - 1], t[j - 1]), E[i][j], F[i][j]});
      best = std::max(best, H[i][j]);
    }
  }
  return best;
}

struct TaskSet {
  std::vector<std::vector<uint8_t>> q, t;
  std::vector<AlignTask> tasks;
};

// Random banded alignment problems: related and unrelated pairs, Ns, variable lengths,
// short (padded) targets.
TaskSet random_tasks(Rng& rng, size_t n, const ScoreParams& sp, int max_len = 300) {
  TaskSet s;
  for (size_t i = 0; i < n; ++i) {
    const int ql = 1 + static_cast<int>(rng.below(static_cast<uint64_t>(max_len)));
    const double nrate = rng.below(4) == 0 ? 0.02 : 0.0;
    std::vector<uint8_t> q = random_codes(rng, static_cast<size_t>(ql), nrate);
    std::vector<uint8_t> t;
    const int tl = ql + 2 * sp.bw;
    if (rng.below(5) == 0) {
      t = random_codes(rng, static_cast<size_t>(tl), nrate);
    } else {
      // Target = left flank + mutated query + right flank, shifted around the band center.
      const int shift = static_cast<int>(rng.below(static_cast<uint64_t>(2 * sp.bw + 1)));
      t = random_codes(rng, static_cast<size_t>(shift));
      // Separate statements: function-argument evaluation order differs between compilers.
      const double sub_rate = 0.03 * rng.uniform();
      const double indel_rate = 0.02 * rng.uniform();
      const std::vector<uint8_t> m = mutate(rng, q, sub_rate, indel_rate);
      t.insert(t.end(), m.begin(), m.end());
      const std::vector<uint8_t> tail = random_codes(rng, static_cast<size_t>(std::max(0, tl - (int)t.size())));
      t.insert(t.end(), tail.begin(), tail.end());
    }
    if (rng.below(10) == 0) t.resize(rng.below(t.size() + 1));  // exercise N padding
    s.q.push_back(std::move(q));
    s.t.push_back(std::move(t));
  }
  for (size_t i = 0; i < n; ++i)
    s.tasks.push_back({s.q[i].data(), static_cast<int32_t>(s.q[i].size()), s.t[i].data(),
                       static_cast<int32_t>(s.t[i].size())});
  return s;
}

const ScoreParams kScoringSchemes[] = {
    ScoreParams{},                                       // default 2/4/4/2, bw 16
    ScoreParams{1, 4, 6, 1, 1, 8},                       // BWA-like
    ScoreParams{2, 8, 12, 2, 1, 5},                      // minimap2 -x sr like
    ScoreParams{3, 2, 0, 3, 1, 21},                      // zero gap-open (linear gaps)
};

}  // namespace

// ---------------------------------------------------------------- tests
TEST(encode_and_revcomp) {
  const std::string s = "ACGTNacgtx";
  std::vector<uint8_t> c(s.size()), rc(s.size());
  encode(s.data(), s.size(), c.data());
  const uint8_t expect[] = {0, 1, 2, 3, 4, 0, 1, 2, 3, 4};
  for (size_t i = 0; i < s.size(); ++i) CHECK_EQ(c[i], expect[i]);
  revcomp_codes(c.data(), c.size(), rc.data());
  for (size_t i = 0; i < s.size(); ++i) CHECK_EQ(rc[i], comp_code(c[c.size() - 1 - i]));
  std::string r;
  revcomp_append("AACGTN", r);
  CHECK(r == "NACGTT");
}

TEST(cigar_roundtrip) {
  std::vector<uint32_t> c;
  CHECK(parse_cigar("5S10M2I3D20M", 12, c));
  CHECK(cigar_to_string(c) == "5S10M2I3D20M");
  CHECK(!parse_cigar("10M5", 4, c));
  CHECK(!parse_cigar("M", 1, c));
}

TEST(minimizers_match_brute_force) {
  Rng rng(7);
  for (int trial = 0; trial < 200; ++trial) {
    const int k = 5 + static_cast<int>(rng.below(20)), w = 1 + static_cast<int>(rng.below(20));
    const std::vector<uint8_t> s = random_codes(rng, 50 + rng.below(400), trial % 3 == 0 ? 0.01 : 0.0);
    std::vector<Minimizer> got;
    sketch(s.data(), static_cast<uint32_t>(s.size()), k, w, 1000, got);

    // Brute force: canonical hash of every k-mer, leftmost minimum of each full window
    // of w k-mers that lies within an N-free run.
    const uint64_t mask = (1ULL << (2 * k)) - 1;
    struct K {
      int64_t start;
      uint64_t h;
      uint32_t strand;
      bool valid;
    };
    std::vector<std::vector<K>> runs(1);
    for (size_t i = 0; i < s.size(); ++i) {
      if (s[i] > 3) {
        runs.emplace_back();
        continue;
      }
      if (i + 1 < static_cast<size_t>(k)) continue;
      bool hasn = false;
      for (int j = 0; j < k; ++j) hasn |= s[i - j] > 3;
      if (hasn) continue;
      uint64_t f = 0, r = 0;
      for (int j = 0; j < k; ++j) {
        f = f << 2 | s[i + 1 - k + j];
        r = r << 2 | (3 - s[i - j]);
      }
      runs.back().push_back({static_cast<int64_t>(i + 1 - k), hash64(std::min(f, r), mask), r < f ? 1u : 0u, f != r});
    }
    std::vector<std::pair<int64_t, uint64_t>> expect;
    for (const auto& run : runs) {
      for (size_t st = 0; st + static_cast<size_t>(w) <= run.size(); ++st) {
        const K* best = nullptr;
        for (size_t j = st; j < st + static_cast<size_t>(w); ++j)
          if (run[j].valid && (!best || run[j].h < best->h)) best = &run[j];
        if (best && (expect.empty() || expect.back().first != best->start + 1000))
          expect.push_back({best->start + 1000, best->h});
      }
    }
    CHECK_EQ(got.size(), expect.size());
    for (size_t i = 0; i < std::min(got.size(), expect.size()); ++i) {
      CHECK_EQ(got[i].pos, expect[i].first);
      CHECK(got[i].hash == expect[i].second);
    }
  }
}

TEST(minimizers_are_strand_symmetric) {
  // A sequence and its reverse complement select the same canonical minimizer hashes.
  Rng rng(11);
  const std::vector<uint8_t> s = random_codes(rng, 5000);
  std::vector<uint8_t> rc(s.size());
  revcomp_codes(s.data(), s.size(), rc.data());
  std::vector<Minimizer> a, b;
  sketch(s.data(), static_cast<uint32_t>(s.size()), 15, 10, 0, a);
  sketch(rc.data(), static_cast<uint32_t>(rc.size()), 15, 10, 0, b);
  std::vector<uint64_t> ha, hb;
  for (auto& m : a) ha.push_back(m.hash);
  for (auto& m : b) hb.push_back(m.hash);
  std::sort(ha.begin(), ha.end());
  std::sort(hb.begin(), hb.end());
  size_t common = 0;
  for (size_t i = 0, j = 0; i < ha.size() && j < hb.size();) {
    if (ha[i] == hb[j]) ++common, ++i, ++j;
    else if (ha[i] < hb[j]) ++i;
    else ++j;
  }
  CHECK(common >= 0.95 * ha.size());
}

TEST(index_finds_every_reference_minimizer) {
  Rng rng(3);
  Reference ref;
  for (int c = 0; c < 3; ++c) {
    const std::vector<uint8_t> s = random_codes(rng, 20000 + rng.below(5000), 0.001);
    ref.contigs.push_back({"c" + std::to_string(c), ref.seq.size(), static_cast<uint32_t>(s.size())});
    ref.seq.insert(ref.seq.end(), s.begin(), s.end());
  }
  MinimizerIndex idx;
  idx.build(ref, 15, 10);
  std::vector<Minimizer> mins;
  for (const Contig& c : ref.contigs)
    sketch(ref.seq.data() + c.offset, c.len, 15, 10, static_cast<uint32_t>(c.offset), mins);
  CHECK_EQ(idx.n_positions(), mins.size());
  size_t found = 0;
  for (const Minimizer& m : mins) {
    uint32_t n;
    const uint32_t* p = idx.lookup(m.hash, n);
    for (uint32_t i = 0; i < n; ++i) found += p[i] == (m.pos << 1 | m.strand);
  }
  CHECK_EQ(found, mins.size());
  uint32_t n = 1;
  CHECK(idx.lookup(1ULL << 40, n) == nullptr);  // k=15 hashes are < 2^30, so this key is absent
  CHECK_EQ(n, 0u);
}

TEST(scalar_kernel_matches_naive_dp) {
  Rng rng(5);
  for (const ScoreParams& sp : kScoringSchemes) {
    const TaskSet ts = random_tasks(rng, 300, sp, 120);
    ScalarAligner sa(sp);
    std::vector<AlignResult> r(ts.tasks.size());
    sa.align(ts.tasks.data(), ts.tasks.size(), false, r.data());
    for (size_t i = 0; i < ts.tasks.size(); ++i) CHECK_EQ(r[i].score, naive_banded_sw(ts.q[i], ts.t[i], sp));
  }
}

TEST(traceback_is_consistent) {
  Rng rng(9);
  for (const ScoreParams& sp : kScoringSchemes) {
    const TaskSet ts = random_tasks(rng, 1000, sp);
    ScalarAligner sa(sp);
    std::vector<AlignResult> tb(ts.tasks.size()), so(ts.tasks.size());
    sa.align(ts.tasks.data(), ts.tasks.size(), true, tb.data());
    sa.align(ts.tasks.data(), ts.tasks.size(), false, so.data());
    for (size_t i = 0; i < ts.tasks.size(); ++i) {
      const AlignTask& t = ts.tasks[i];
      const AlignResult& r = tb[i];
      CHECK_EQ(r.score, so[i].score);
      if (r.score == 0) continue;
      CHECK_EQ(score_cigar(t.q, t.t, r.qb, r.tb, r.cigar, sp), r.score);
      int qspan = 0, tspan = 0;
      for (uint32_t c : r.cigar) {
        if ((c & 0xf) != kCigD) qspan += static_cast<int>(c >> 4);
        if ((c & 0xf) != kCigI) tspan += static_cast<int>(c >> 4);
      }
      CHECK_EQ(qspan, r.qe - r.qb);
      CHECK_EQ(tspan, r.te - r.tb);
      CHECK((r.cigar.front() & 0xf) == kCigM && (r.cigar.back() & 0xf) == kCigM);
    }
  }
}

TEST(simd_kernel_bit_identical_to_scalar) {
  Rng rng(13);
  for (const ScoreParams& sp : kScoringSchemes) {
    for (size_t n : {size_t{1}, size_t{5}, static_cast<size_t>(SimdAligner::lanes()), size_t{1003}}) {
      const TaskSet ts = random_tasks(rng, n, sp);
      ScalarAligner sa(sp);
      SimdAligner va(sp);
      for (bool tb : {false, true}) {
        std::vector<AlignResult> a(n), b(n);
        sa.align(ts.tasks.data(), n, tb, a.data());
        va.align(ts.tasks.data(), n, tb, b.data());
        for (size_t i = 0; i < n; ++i) {
          CHECK_EQ(a[i].score, b[i].score);
          if (!tb) continue;
          CHECK_EQ(a[i].qb, b[i].qb);
          CHECK_EQ(a[i].qe, b[i].qe);
          CHECK_EQ(a[i].tb, b[i].tb);
          CHECK_EQ(a[i].te, b[i].te);
          CHECK_EQ(a[i].nm, b[i].nm);
          CHECK(a[i].cigar == b[i].cigar);
        }
      }
    }
  }
}

TEST(perfect_match_alignment) {
  ScoreParams sp;
  Rng rng(17);
  const std::vector<uint8_t> q = random_codes(rng, 150);
  std::vector<uint8_t> t = random_codes(rng, static_cast<size_t>(sp.bw));
  t.insert(t.end(), q.begin(), q.end());
  const std::vector<uint8_t> tail = random_codes(rng, static_cast<size_t>(sp.bw));
  t.insert(t.end(), tail.begin(), tail.end());
  const AlignTask task{q.data(), 150, t.data(), static_cast<int32_t>(t.size())};
  for (int simd = 0; simd < 2; ++simd) {
    AlignResult r;
    if (simd) SimdAligner(sp).align(&task, 1, true, &r);
    else ScalarAligner(sp).align(&task, 1, true, &r);
    CHECK_EQ(r.score, 150 * sp.match);
    CHECK_EQ(r.qb, 0);
    CHECK_EQ(r.qe, 150);
    CHECK_EQ(r.tb, sp.bw);
    CHECK_EQ(r.nm, 0);
    CHECK(cigar_to_string(r.cigar) == "150M");
  }
}

TEST(chaining_finds_colinear_anchors) {
  // Anchors on one diagonal with a 3 bp deletion plus scattered noise.
  std::vector<Anchor> a;
  for (int i = 0; i < 10; ++i) a.push_back({static_cast<uint32_t>(1000 + i * 12 + (i >= 5 ? 3 : 0)), i * 12, 0});
  a.push_back({50000, 30, 0});
  a.push_back({90000, 60, 1});
  Chainer ch;
  std::vector<Chain> out;
  ch.run(a, 15, ChainParams{}, out);
  CHECK(!out.empty());
  CHECK_EQ(out[0].n_anchors, 10);
  CHECK_EQ(out[0].dmin, 1000);
  CHECK_EQ(out[0].dmax, 1003);
}

TEST(end_to_end_mapping_accuracy) {
  namespace fs = std::filesystem;
  const fs::path dir = fs::temp_directory_path() / ("dnaln_test_" + std::to_string(::getpid()));
  fs::create_directories(dir);
  const std::string ref_path = (dir / "ref.fa").string(), prefix = (dir / "sim").string();
  // Full E. coli-sized genome: the fixed repeat families are then ~1% of it, as in E. coli
  // (on a tiny genome they would dominate and make many reads unplaceable).
  random_reference(4641652, 0.508, 21, ref_path);
  const Reference ref = Reference::load(ref_path);
  SimParams sim;
  sim.n = 4000;
  simulate_reads(ref, sim, prefix);
  const std::vector<SeqRecord> reads = read_all_records(prefix + ".fq");
  const TruthMap truth = load_truth(prefix + ".truth.tsv");
  CHECK_EQ(reads.size(), sim.n);
  CHECK_EQ(truth.size(), sim.n);

  MinimizerIndex idx;
  MapParams mp;
  idx.build(ref, mp.k, mp.w);
  double acc[2];
  for (int simd = 0; simd < 2; ++simd) {
    mp.simd = simd;
    Mapper m(ref, idx, mp);
    std::vector<MapRecord> recs;
    m.map(reads.data(), reads.size(), recs);
    Evaluator ev(truth);
    for (size_t i = 0; i < reads.size(); ++i)
      ev.add(reads[i].name, recs[i].mapped(), recs[i].mapped() ? ref.contigs[recs[i].contig].name : "*",
             recs[i].pos, recs[i].rev, recs[i].mapq, recs[i].cigar);
    acc[simd] = ev.accuracy();
    CHECK(ev.accuracy() > 0.98);
    CHECK(ev.base_concordance() > 0.98);
    CHECK(ev.wrong_at(60) <= ev.mapped_at(60) / 200);  // MAPQ 60 must be reliable (<= 0.5% wrong)
  }
  CHECK(acc[0] == acc[1]);  // kernels are bit-identical, so the mapping must be too
  std::printf("    (end-to-end accuracy on 4,000 simulated reads: %.2f%%)\n", 100.0 * acc[1]);
  fs::remove_all(dir);
}

int main() {
  std::printf("dnaln unit tests (SIMD: %s, %d lanes)\n", SimdAligner::isa(), SimdAligner::lanes());
  int failed_tests = 0;
  for (const TestCase& t : registry()) {
    const int before = g_failures;
    Timer timer;
    t.fn();
    const bool ok = g_failures == before;
    failed_tests += !ok;
    std::printf("  [%s] %-40s %7.1f ms\n", ok ? " OK " : "FAIL", t.name, timer.sec() * 1e3);
  }
  std::printf("%s: %zu tests, %d checks, %d failures\n", failed_tests ? "FAILED" : "PASSED", registry().size(),
              g_checks, g_failures);
  return failed_tests ? 1 : 0;
}
