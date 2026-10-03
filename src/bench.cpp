// bench.cpp — automated benchmark harness (`dnaln bench`).
//
//   1. Kernel: scalar baseline vs SIMD on one banded alignment per simulated read,
//      placed at the read's true locus (score-only and traceback modes).
//   2. Correctness on the same ground-truth alignments: SIMD == scalar bit-for-bit,
//      CIGARs re-score to the reported score, SW score >= score of the true path
//      (optimality), and base-level agreement with the true alignment.
//   3. End-to-end mapping throughput (single thread) and accuracy vs ground truth,
//      SIMD vs scalar kernel, software prefetch on vs off.
//   4. Cache sweep: kernel throughput as the band (and so the working set) grows
//      past L1 and L2.
#include <algorithm>
#include <cstdio>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#if defined(__APPLE__)
#include <sys/sysctl.h>
#endif
#include <unistd.h>

#include "cli.hpp"
#include "evaluate.hpp"
#include "groundtruth.hpp"

namespace dnaln {

namespace {

struct HwInfo {
  std::string cpu = "unknown CPU";
  long long l1d = 0, l2 = 0;
};

HwInfo hw_info() {
  HwInfo h;
#if defined(__APPLE__)
  char buf[256];
  size_t len = sizeof buf;
  if (sysctlbyname("machdep.cpu.brand_string", buf, &len, nullptr, 0) == 0) h.cpu = buf;
  auto get = [](const char* a, const char* b) {
    long long v = 0;
    size_t n = sizeof v;
    if (sysctlbyname(a, &v, &n, nullptr, 0) == 0 && v > 0) return v;
    v = 0;
    n = sizeof v;
    return sysctlbyname(b, &v, &n, nullptr, 0) == 0 ? v : 0LL;
  };
  h.l1d = get("hw.perflevel0.l1dcachesize", "hw.l1dcachesize");  // performance cores
  h.l2 = get("hw.perflevel0.l2cachesize", "hw.l2cachesize");
#elif defined(__linux__)
  if (FILE* fp = std::fopen("/proc/cpuinfo", "r")) {
    char line[512];
    while (std::fgets(line, sizeof line, fp)) {
      if (std::strncmp(line, "model name", 10) == 0 || std::strncmp(line, "Model", 5) == 0) {
        const char* c = std::strchr(line, ':');
        if (c) {
          h.cpu = c + 2;
          if (!h.cpu.empty() && h.cpu.back() == '\n') h.cpu.pop_back();
        }
        break;
      }
    }
    std::fclose(fp);
  }
#if defined(_SC_LEVEL1_DCACHE_SIZE)
  h.l1d = sysconf(_SC_LEVEL1_DCACHE_SIZE);
  h.l2 = sysconf(_SC_LEVEL2_CACHE_SIZE);
#endif
#endif
  return h;
}

std::string human_bytes(double b) {
  char buf[32];
  if (b >= 1048576) std::snprintf(buf, sizeof buf, "%.1f MB", b / 1048576);
  else if (b >= 1024) std::snprintf(buf, sizeof buf, "%.1f KB", b / 1024);
  else std::snprintf(buf, sizeof buf, "%.0f B", b);
  return buf;
}

std::string with_commas(unsigned long long v) {
  std::string s = std::to_string(v);
  for (int i = static_cast<int>(s.size()) - 3; i > 0; i -= 3) s.insert(static_cast<size_t>(i), ",");
  return s;
}

template <class F>
double best_of(int reps, F&& f) {
  double best = 1e300;
  for (int r = 0; r < reps; ++r) {
    Timer t;
    f();
    best = std::min(best, t.sec());
  }
  return best;
}

// Walks the true CIGAR inside the window: its score, whether it stays in the band, and
// the true window column of every query base (-1 for inserted bases).
bool truth_path(const AlignTask& tk, const TruthRec& tr, const ScoreParams& sp, int& score,
                std::vector<int>& cols) {
  int i = 0, j = sp.bw, s = 0;
  bool in_band = true;
  const int W = sp.band_cols();
  cols.clear();
  for (uint32_t c : tr.cigar) {
    const uint32_t op = c & 0xf;
    const int len = static_cast<int>(c >> 4);
    if (op == kCigM) {
      for (int k = 0; k < len; ++k, ++i, ++j) {
        if (j - i < 0 || j - i >= W || i >= tk.qlen || j >= tk.tlen) {
          in_band = false;
          continue;
        }
        s += sp.sub(tk.q[i], tk.t[j]);
        cols.push_back(j);
      }
    } else if (op == kCigI) {
      s -= sp.gap_open + len * sp.gap_ext;
      cols.insert(cols.end(), static_cast<size_t>(len), -1);
      i += len;
    } else if (op == kCigD) {
      s -= sp.gap_open + len * sp.gap_ext;
      j += len;
    }
  }
  score = s;
  return in_band;
}

void predicted_cols(const AlignResult& r, int qlen, std::vector<int>& cols) {
  cols.assign(static_cast<size_t>(qlen), -1);
  int i = r.qb, j = r.tb;
  for (uint32_t c : r.cigar) {
    const uint32_t op = c & 0xf;
    const int len = static_cast<int>(c >> 4);
    if (op == kCigM) {
      for (int k = 0; k < len; ++k) cols[static_cast<size_t>(i++)] = j++;
    } else if (op == kCigI) {
      i += len;
    } else if (op == kCigD) {
      j += len;
    }
  }
}

bool same_result(const AlignResult& a, const AlignResult& b, bool tb) {
  if (a.score != b.score) return false;
  return !tb || (a.qb == b.qb && a.qe == b.qe && a.tb == b.tb && a.te == b.te && a.nm == b.nm && a.cigar == b.cigar);
}

struct KernelTiming {
  double scalar_sec = 0, simd_sec = 0;
  size_t identical = 0;
  double speedup() const { return simd_sec > 0 ? scalar_sec / simd_sec : 0; }
};

KernelTiming time_kernels(const ScoreParams& sp, const std::vector<AlignTask>& tasks, bool tb, int reps,
                          std::vector<AlignResult>& rs, std::vector<AlignResult>& rv) {
  ScalarAligner scalar(sp);
  SimdAligner simd(sp);
  rs.resize(tasks.size());
  rv.resize(tasks.size());
  KernelTiming k;
  k.scalar_sec = best_of(reps, [&] { scalar.align(tasks.data(), tasks.size(), tb, rs.data()); });
  k.simd_sec = best_of(reps, [&] { simd.align(tasks.data(), tasks.size(), tb, rv.data()); });
  for (size_t i = 0; i < tasks.size(); ++i) k.identical += same_result(rs[i], rv[i], tb);
  return k;
}

struct MapRun {
  double sec = 0;
  MapStats st;
  std::vector<MapRecord> recs;
};

MapRun run_mapping(const Reference& ref, const MinimizerIndex& idx, const MapParams& p,
                   const std::vector<SeqRecord>& reads, std::string* sam) {
  Mapper m(ref, idx, p);
  MapRun r;
  r.recs.resize(reads.size());
  std::vector<MapRecord> batch;
  std::string buf;
  const size_t B = 8192;
  Timer t;
  for (size_t i = 0; i < reads.size(); i += B) {
    const size_t n = std::min(B, reads.size() - i);
    m.map(&reads[i], n, batch);
    buf.clear();
    for (size_t j = 0; j < n; ++j) m.append_sam(reads[i + j], batch[j], buf);  // SAM formatting is timed
    if (sam) sam->append(buf);
    for (size_t j = 0; j < n; ++j) std::swap(r.recs[i + j], batch[j]);
  }
  r.sec = t.sec();
  r.st = m.stats();
  return r;
}

const char* verdict(bool ok) { return ok ? "PASS" : "FAIL"; }

}  // namespace

int cmd_bench(int argc, char** argv) {
  std::set<std::string> values = kMapValueOpts;
  values.insert({"--reps", "--limit", "--sweep-limit", "--sam", "--tol"});
  const Args a(argc, argv, values, {"--skip-map", "--skip-sweep", "--kernel-only", "--scalar", "--no-prefetch",
                                    "-h", "--help"});
  if (a.pos.size() != 3 || a.has("-h") || a.has("--help")) {
    std::fprintf(stderr,
                 "Usage: dnaln bench [options] <ref.fa> <reads.fq> <reads.truth.tsv>\n\n"
                 "  --reps INT         timing repetitions, best is reported [3]\n"
                 "  --limit INT        max ground-truth alignments for the kernel benchmark [all reads]\n"
                 "  --sweep-limit INT  alignments per point of the cache sweep [20000]\n"
                 "  --sam FILE         also write the SIMD mapper's SAM output to FILE\n"
                 "  --kernel-only      only run sections 1-2 (useful under a profiler)\n"
                 "  --skip-map / --skip-sweep\n%s",
                 kMapOptsHelp);
    return a.has("-h") || a.has("--help") ? 0 : 1;
  }
  const MapParams mp = map_params_from_args(a);
  const ScoreParams& sp = mp.sc;
  const int reps = static_cast<int>(std::max(1LL, a.num("--reps", 3)));
  const bool kernel_only = a.has("--kernel-only");

  const HwInfo hw = hw_info();
  std::printf("dnaln benchmark harness\n");
  std::printf("  CPU: %s | L1d %s, L2 %s | SIMD: %s, %d x int16 lanes per alignment group\n", hw.cpu.c_str(),
              hw.l1d ? human_bytes(static_cast<double>(hw.l1d)).c_str() : "?",
              hw.l2 ? human_bytes(static_cast<double>(hw.l2)).c_str() : "?", SimdAligner::isa(),
              SimdAligner::lanes());
  std::printf("  scoring: match %d, mismatch -%d, gap open %d + extend %d/base, band +/-%d (W=%d)\n", sp.match,
              sp.mismatch, sp.gap_open, sp.gap_ext, sp.bw, sp.band_cols());

  Timer t_load;
  const Reference ref = Reference::load(a.pos[0]);
  const std::vector<SeqRecord> reads = read_all_records(a.pos[1]);
  const TruthMap truth = load_truth(a.pos[2]);
  std::printf("  data: %zu bp reference (%zu contigs), %s reads, %s truth records (loaded in %.2f s)\n\n",
              ref.seq.size(), ref.contigs.size(), with_commas(reads.size()).c_str(),
              with_commas(truth.size()).c_str(), t_load.sec());

  // ------------------------------------------------------------------ 1. kernel
  const size_t limit = static_cast<size_t>(a.num("--limit", static_cast<long long>(reads.size())));
  const GtSet gt = build_gt(ref, reads, truth, sp.bw, limit);
  const size_t n = gt.tasks.size();
  if (n == 0) die("no reads matched the truth file");
  uint64_t cells = 0;
  for (int q : gt.qlen) cells += static_cast<uint64_t>(q) * sp.band_cols();

  std::printf("== 1. Banded Smith-Waterman kernel: scalar baseline vs SIMD ==\n");
  std::printf("   %s alignment tasks (one per read, at its true locus), %s DP cells each pass, best of %d\n\n",
              with_commas(n).c_str(), with_commas(cells).c_str(), reps);
  std::printf("   %-11s %-18s %10s %9s %9s %22s\n", "mode", "kernel", "time (ms)", "GCUPS", "speedup",
              "identical to scalar");
  std::vector<AlignResult> rs, rv;
  KernelTiming kt[2];
  for (int tb = 0; tb < 2; ++tb) {
    kt[tb] = time_kernels(sp, gt.tasks, tb, reps, rs, rv);
    const char* mode = tb ? "traceback" : "score-only";
    std::printf("   %-11s %-18s %10.1f %9.2f %9s %22s\n", mode, "scalar", kt[tb].scalar_sec * 1e3,
                cells / kt[tb].scalar_sec * 1e-9, "1.00x", "-");
    char sx[16], id[48];
    std::snprintf(sx, sizeof sx, "%.2fx", kt[tb].speedup());
    std::snprintf(id, sizeof id, "%s / %s", with_commas(kt[tb].identical).c_str(), with_commas(n).c_str());
    std::printf("   %-11s %-18s %10.1f %9.2f %9s %22s\n", mode, SimdAligner(sp).name().c_str(),
                kt[tb].simd_sec * 1e3, cells / kt[tb].simd_sec * 1e-9, sx, id);
  }
  std::printf("\n");

  // ------------------------------------------------------- 2. correctness vs truth
  // rs / rv now hold the traceback-mode results of the last timing pass.
  size_t rescored_ok = 0, in_band = 0, optimal_ok = 0, majority_ok = 0;
  uint64_t bases = 0, bases_ok = 0;
  std::vector<int> tcols, pcols;
  for (size_t i = 0; i < n; ++i) {
    const AlignTask& tk = gt.tasks[i];
    const AlignResult& r = rv[i];
    rescored_ok += score_cigar(tk.q, tk.t, r.qb, r.tb, r.cigar, sp) == r.score;
    int tscore = 0;
    if (truth_path(tk, *gt.truth[i], sp, tscore, tcols)) {
      ++in_band;
      optimal_ok += r.score >= tscore;
    }
    predicted_cols(r, tk.qlen, pcols);
    size_t ok = 0, tot = 0;
    for (size_t b = 0; b < tcols.size() && b < pcols.size(); ++b) {
      if (tcols[b] < 0) continue;
      ++tot;
      ok += pcols[b] == tcols[b];
    }
    bases += tot;
    bases_ok += ok;
    majority_ok += tot && ok >= 0.95 * tot;
  }
  const bool simd_eq = kt[0].identical == n && kt[1].identical == n;
  std::printf("== 2. Correctness on %s ground-truth simulated alignments ==\n\n", with_commas(n).c_str());
  std::printf("   %-58s %21s  %s\n", "SIMD == scalar (score, coordinates, CIGAR, NM; both modes)",
              (with_commas(std::min(kt[0].identical, kt[1].identical)) + " / " + with_commas(n)).c_str(),
              verdict(simd_eq));
  std::printf("   %-58s %21s  %s\n", "CIGAR re-scores to the reported SW score",
              (with_commas(rescored_ok) + " / " + with_commas(n)).c_str(), verdict(rescored_ok == n));
  std::printf("   %-58s %21s  %s\n", "Optimality: SW score >= score of the true alignment path",
              (with_commas(optimal_ok) + " / " + with_commas(in_band)).c_str(), verdict(optimal_ok == in_band));
  std::printf("   %-58s %21s\n", "  (true path fully inside the +/-bw band)",
              (with_commas(in_band) + " / " + with_commas(n)).c_str());
  std::printf("   %-58s %20.3f%%\n", "Read bases aligned to their true reference position",
              bases ? 100.0 * bases_ok / bases : 0.0);
  std::printf("   %-58s %20.3f%%\n", "Alignments agreeing with truth on >= 95% of bases",
              100.0 * majority_ok / n);
  std::printf("\n");

  if (kernel_only) return simd_eq && rescored_ok == n && optimal_ok == in_band ? 0 : 2;

  // ------------------------------------------------------------- 3. mapping
  double best_rps = 0, accuracy = 0;
  if (!a.has("--skip-map")) {
    std::printf("== 3. End-to-end mapping, single thread (%s reads; excludes index build and file I/O) ==\n\n",
                with_commas(reads.size()).c_str());
    Timer t_idx;
    MinimizerIndex idx;
    idx.build(ref, mp.k, mp.w);
    std::printf("   index: k=%d w=%d | %s distinct minimizers, %s positions, max occurrence %u | %s | built in "
                "%.2f s\n\n",
                mp.k, mp.w, with_commas(idx.n_keys()).c_str(), with_commas(idx.n_positions()).c_str(),
                idx.max_count(), human_bytes(static_cast<double>(idx.memory_bytes())).c_str(), t_idx.sec());

    struct Config {
      const char* label;
      bool simd, prefetch;
    };
    const Config configs[] = {{"SIMD kernel, prefetch", true, true},
                              {"SIMD kernel, no prefetch", true, false},
                              {"scalar kernel, prefetch", false, true}};
    std::printf("   %-26s %12s %10s %10s %10s %9s %10s\n", "configuration", "reads/s", "seed (s)", "align (s)",
                "other (s)", "GCUPS", "accuracy");
    std::string sam;
    std::unique_ptr<Evaluator> main_eval;
    for (const Config& c : configs) {
      MapParams p = mp;
      p.simd = c.simd;
      p.prefetch = c.prefetch;
      MapRun best;
      for (int r = 0; r < reps; ++r) {
        const bool keep_sam = a.has("--sam") && c.simd && c.prefetch && r == 0;
        MapRun run = run_mapping(ref, idx, p, reads, keep_sam ? &sam : nullptr);
        if (r == 0 || run.sec < best.sec) best = std::move(run);
      }
      Evaluator ev(truth, static_cast<int>(a.num("--tol", 10)));
      for (size_t i = 0; i < reads.size(); ++i) {
        const MapRecord& m = best.recs[i];
        ev.add(reads[i].name, m.mapped(), m.mapped() ? ref.contigs[m.contig].name : "*", m.pos, m.rev, m.mapq,
               m.cigar);
      }
      const double rps = reads.size() / best.sec;
      std::printf("   %-26s %12s %10.3f %10.3f %10.3f %9.2f %9.3f%%\n", c.label,
                  with_commas(static_cast<unsigned long long>(rps)).c_str(), best.st.t_seed, best.st.t_align,
                  best.sec - best.st.t_seed - best.st.t_align, best.st.cells / best.st.t_align * 1e-9,
                  100.0 * ev.accuracy());
      if (c.simd && c.prefetch) {
        best_rps = rps;
        accuracy = ev.accuracy();
        main_eval = std::make_unique<Evaluator>(ev);
      }
    }
    std::printf("\n   seed = minimizers + index lookups + chaining; align = banded SW (ranking + traceback);\n"
                "   other = SAM formatting and bookkeeping. Best of %d runs.\n\n", reps);
    std::printf("   Accuracy of the SIMD run against ground truth:\n\n");
    print_eval_report(stdout, {"dnaln"}, {main_eval.get()});
    std::printf("\n");
    if (a.has("--sam")) {
      const std::string path = a.str("--sam", "");
      FILE* fp = std::fopen(path.c_str(), "w");
      if (!fp) die("cannot write %s", path.c_str());
      Mapper(ref, idx, mp).write_sam_header(fp, "dnaln bench");
      std::fwrite(sam.data(), 1, sam.size(), fp);
      std::fclose(fp);
    }
  }

  // ------------------------------------------------------------ 4. cache sweep
  if (!a.has("--skip-sweep")) {
    const size_t sl = std::min(n, static_cast<size_t>(a.num("--sweep-limit", 20000)));
    std::printf("== 4. Cache behavior: kernel throughput vs band width / working set (%s alignments) ==\n\n",
                with_commas(sl).c_str());
    std::printf("   %5s %5s %16s %18s %13s %11s %11s %11s\n", "bw", "W", "SIMD WS (score)", "SIMD WS (trace)",
                "scalar GCUPS", "SIMD GCUPS", "x (score)", "x (trace)");
    const int median_len = gt.qlen[gt.qlen.size() / 2];
    for (int bw : {4, 8, 16, 32, 64, 128}) {
      ScoreParams s2 = sp;
      s2.bw = bw;
      const GtSet g2 = build_gt(ref, reads, truth, bw, sl);
      uint64_t c2 = 0;
      for (int q : g2.qlen) c2 += static_cast<uint64_t>(q) * s2.band_cols();
      std::vector<AlignResult> x, y;
      const KernelTiming k0 = time_kernels(s2, g2.tasks, false, reps, x, y);
      const KernelTiming k1 = time_kernels(s2, g2.tasks, true, reps, x, y);
      SimdAligner probe(s2);
      const double ws0 = static_cast<double>(probe.working_set_bytes(median_len, false));
      const double ws1 = static_cast<double>(probe.working_set_bytes(median_len, true));
      auto where = [&](double ws) {
        return std::string(hw.l1d && ws <= hw.l1d ? " L1" : hw.l2 && ws <= hw.l2 ? " L2" : " >L2");
      };
      std::printf("   %5d %5d %16s %18s %13.2f %11.2f %10.2fx %10.2fx\n", bw, s2.band_cols(),
                  (human_bytes(ws0) + where(ws0)).c_str(), (human_bytes(ws1) + where(ws1)).c_str(),
                  c2 / k0.scalar_sec * 1e-9, c2 / k0.simd_sec * 1e-9, k0.speedup(), k1.speedup());
    }
    std::printf("\n   WS = kernel scratch per %d-lane group for a %d bp read (H/E rows, transposed sequences,\n"
                "   and 1 byte/cell of traceback directions); fits-in level is vs this machine's caches.\n\n",
                SimdAligner::lanes(), median_len);
  }

  // ----------------------------------------------------------------- summary
  const double sx = std::min(kt[0].speedup(), kt[1].speedup());
  std::printf("== Summary vs project targets ==\n\n");
  if (best_rps > 0)
    std::printf("   %-52s %14s  %s\n", "Mapping throughput, 1 thread (target >= 20,000 reads/s)",
                with_commas(static_cast<unsigned long long>(best_rps)).c_str(), verdict(best_rps >= 20000));
  char buf[64];
  std::snprintf(buf, sizeof buf, "%.2fx / %.2fx", kt[0].speedup(), kt[1].speedup());
  std::printf("   %-52s %14s  %s\n", "SIMD kernel speedup, score-only / traceback (>= 6x)", buf, verdict(sx >= 6.0));
  std::snprintf(buf, sizeof buf, "%s", with_commas(n).c_str());
  std::printf("   %-52s %14s  %s\n", "Ground-truth alignments validated (SIMD==scalar, optimal)", buf,
              verdict(simd_eq && rescored_ok == n && optimal_ok == in_band));
  if (best_rps > 0) {
    std::snprintf(buf, sizeof buf, "%.3f%%", 100.0 * accuracy);
    std::printf("   %-52s %14s  %s\n", "Mapping accuracy (compare to minimap2 with `dnaln eval`)", buf, "-");
  }
  std::printf("\n");
  return 0;
}

}  // namespace dnaln
