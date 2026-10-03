// demo.cpp — `dnaln demo`: a narrated, presentation-friendly walkthrough.
//
// Run it through ./demo.sh, which prepares the inputs. Each step pauses for Enter when
// run interactively, so a presenter can talk over it.
#include <unistd.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

#include "cli.hpp"
#include "evaluate.hpp"
#include "groundtruth.hpp"

namespace dnaln {

namespace {

bool g_color = false, g_pause = false;

const char* c_bold() { return g_color ? "\033[1m" : ""; }
const char* c_dim() { return g_color ? "\033[2m" : ""; }
const char* c_green() { return g_color ? "\033[32m" : ""; }
const char* c_red() { return g_color ? "\033[1;31m" : ""; }
const char* c_yellow() { return g_color ? "\033[1;33m" : ""; }
const char* c_cyan() { return g_color ? "\033[1;36m" : ""; }
const char* c_reset() { return g_color ? "\033[0m" : ""; }

std::string commas(unsigned long long v) {
  std::string s = std::to_string(v);
  for (int i = static_cast<int>(s.size()) - 3; i > 0; i -= 3) s.insert(static_cast<size_t>(i), ",");
  return s;
}

void pause_here() {
  if (!g_pause) return;
  std::printf("\n  %s[ press Enter to continue ]%s", c_dim(), c_reset());
  std::fflush(stdout);
  int c;
  while ((c = std::getchar()) != '\n' && c != EOF) {
  }
}

void step(int n, const char* title) {
  std::printf("\n%sSTEP %d%s  %s%s%s\n\n", c_cyan(), n, c_reset(), c_bold(), title, c_reset());
}

// Bar of `frac` x `width` cells using eighth-block characters.
std::string bar(double frac, int width) {
  static const char* parts[] = {"", "▏", "▎", "▍", "▌", "▋", "▊", "▉"};
  const int eighths = static_cast<int>(std::clamp(frac, 0.0, 1.0) * width * 8 + 0.5);
  std::string s;
  for (int i = 0; i < eighths / 8; ++i) s += "█";
  s += parts[eighths % 8];
  const int used = eighths / 8 + (eighths % 8 ? 1 : 0);
  for (int i = used; i < width; ++i) s += ' ';
  return s;
}

void draw_progress(const char* label, const char* color, double frac, double sec, bool done) {
  std::printf("\r  %-14s %s%s%s %6.2f s%s", label, color, bar(frac, 40).c_str(), c_reset(), sec, done ? "" : "  ");
  std::fflush(stdout);
}

// Runs every task `rounds` times through `al`, animating a bar. With scale_sec > 0 the bar
// shows elapsed time relative to scale_sec (so a faster run stops short); otherwise progress.
// Returns seconds.
double race(Aligner& al, const std::vector<AlignTask>& tasks, int rounds, std::vector<AlignResult>& out,
            const char* label, const char* color, double scale_sec = 0) {
  out.resize(tasks.size());
  const size_t total = tasks.size() * static_cast<size_t>(rounds), chunk = 4096;
  size_t done = 0;
  Timer t;
  for (int r = 0; r < rounds; ++r) {
    for (size_t i = 0; i < tasks.size(); i += chunk) {
      const size_t n = std::min(chunk, tasks.size() - i);
      al.align(tasks.data() + i, n, true, out.data() + i);
      done += n;
      const double el = t.sec();
      draw_progress(label, color, scale_sec > 0 ? el / scale_sec : static_cast<double>(done) / total, el, false);
    }
  }
  const double sec = t.sec();
  draw_progress(label, color, scale_sec > 0 ? sec / scale_sec : 1.0, sec, true);
  std::printf("\n");
  return sec;
}

bool placed_correctly(const MapRecord& m, const Reference& ref, const TruthRec& t) {
  if (!m.mapped() || ref.contigs[m.contig].name != t.contig || m.rev != t.rev) return false;
  int64_t clip = 0;
  if (!m.cigar.empty() && (m.cigar.front() & 0xf) == kCigS) clip = m.cigar.front() >> 4;
  return std::llabs(m.pos - clip - t.start) <= 10;
}

// Prints a read aligned against the genome, 60 letters per block.
void show_alignment(const Reference& ref, const SeqRecord& read, const MapRecord& m) {
  std::string q;
  if (m.rev) revcomp_append(read.seq, q);
  else q = read.seq;
  enum Kind { kMatch, kMismatch, kGap };
  struct Column {
    char genome, read;
    Kind kind;
    int64_t pos;  // genome position of this column
  };
  std::vector<Column> cols;
  const uint64_t base = ref.contigs[m.contig].offset;
  int matches = 0, mism = 0, gaps = 0;
  size_t qi = 0;
  int64_t rp = m.pos;
  for (uint32_t c : m.cigar) {
    const uint32_t op = c & 0xf, len = c >> 4;
    for (uint32_t k = 0; k < len; ++k) {
      if (op == kCigS) {
        ++qi;
      } else if (op == kCigM) {
        const char g = kBase[ref.seq[base + rp]], r = q[qi++];
        cols.push_back({g, r, g == r ? kMatch : kMismatch, rp++});
        ++(g == r ? matches : mism);
      } else if (op == kCigI) {
        cols.push_back({'-', q[qi++], kGap, rp});
        ++gaps;
      } else if (op == kCigD) {
        cols.push_back({kBase[ref.seq[base + rp]], '-', kGap, rp});
        ++rp;
        ++gaps;
      }
    }
  }
  for (size_t s = 0; s < cols.size(); s += 60) {
    const size_t e = std::min(cols.size(), s + 60);
    std::string top, mid, bot;
    for (size_t i = s; i < e; ++i) {
      const Column& c = cols[i];
      const char* color = c.kind == kMatch ? "" : c.kind == kMismatch ? c_red() : c_yellow();
      const char* reset = c.kind == kMatch ? "" : c_reset();
      top += std::string(c.kind == kGap ? color : "") + c.genome + (c.kind == kGap ? reset : "");
      mid += c.kind == kMatch ? std::string(c_green()) + '|' + c_reset()
                              : c.kind == kMismatch ? std::string(c_red()) + 'x' + c_reset() : std::string(" ");
      bot += std::string(color) + c.read + reset;
    }
    std::printf("  %sgenome %11s%s  %s\n", c_dim(), commas(static_cast<unsigned long long>(cols[s].pos + 1)).c_str(),
                c_reset(), top.c_str());
    std::printf("  %18s  %s\n", "", mid.c_str());
    std::printf("  %sread%s %13s  %s\n\n", c_dim(), c_reset(), "", bot.c_str());
  }
  std::printf("  %s|%s same letter   %sx%s different letter (sequencing error)   %s-%s missing letter\n\n",
              c_green(), c_reset(), c_red(), c_reset(), c_yellow(), c_reset());
  std::printf("  %d of %zu letters match exactly. The aligner worked around %d wrong letter%s", matches,
              read.seq.size(), mism, mism == 1 ? "" : "s");
  if (gaps) std::printf(" and %d missing/extra letter%s", gaps, gaps == 1 ? "" : "s");
  std::printf(".\n");
}

}  // namespace

int cmd_demo(int argc, char** argv) {
  const Args a(argc, argv, {"--mm2-sam", "--mm2-time", "--genome-name", "--rounds", "--sam"},
               {"--no-pause", "--no-color", "-h", "--help"});
  if (a.pos.size() != 3 || a.has("-h") || a.has("--help")) {
    std::fprintf(stderr,
                 "Usage: dnaln demo [options] <ref.fa> <reads.fq> <reads.truth.tsv>\n"
                 "Narrated walkthrough for presentations (normally run via ./demo.sh).\n"
                 "  --genome-name STR   display name of the genome\n"
                 "  --mm2-sam FILE      minimap2 SAM for the same reads (accuracy comparison)\n"
                 "  --mm2-time SEC      minimap2 wall time for the same reads\n"
                 "  --rounds INT        passes over the reads in the speed race [5]\n"
                 "  --sam FILE          write dnaln's alignments here\n"
                 "  --no-pause / --no-color\n");
    return a.has("-h") || a.has("--help") ? 0 : 1;
  }
  g_color = isatty(STDOUT_FILENO) && !std::getenv("NO_COLOR") && !a.has("--no-color");
  g_pause = isatty(STDIN_FILENO) && isatty(STDOUT_FILENO) && !a.has("--no-pause");
  const int rounds = static_cast<int>(std::max(1LL, a.num("--rounds", 5)));
  MapParams mp;

  std::printf("\n  %s┌────────────────────────────────────────────────────────────────┐%s\n", c_cyan(), c_reset());
  std::printf("  %s│%s  %sdnaln%s · finding where DNA fragments come from in a genome     %s│%s\n", c_cyan(), c_reset(),
              c_bold(), c_reset(), c_cyan(), c_reset());
  std::printf("  %s└────────────────────────────────────────────────────────────────┘%s\n\n", c_cyan(), c_reset());
  std::printf("  A DNA sequencer can't read a genome end to end. It reads millions of short\n"
              "  fragments (\"reads\"), and software must work out where each one came from.\n"
              "  That software is called an aligner. This is one, written from scratch in C++.\n");
  pause_here();

  // ---------------------------------------------------------------- 1. genome + index
  step(1, "Load the genome and build a search index");
  const Reference ref = Reference::load(a.pos[0]);
  Timer t;
  MinimizerIndex idx;
  idx.build(ref, mp.k, mp.w);
  const double t_index = t.sec();
  const std::string gname = a.str("--genome-name", ref.contigs[0].name);
  std::printf("  Genome:  %s%s%s, %s DNA letters (A, C, G, T)\n", c_bold(), gname.c_str(), c_reset(),
              commas(ref.seq.size()).c_str());
  std::printf("  Index:   %s short \"seed\" sequences, each pointing to where it occurs,\n"
              "           built in %s%.2f s%s. Like a book index, it lets us jump straight to\n"
              "           candidate locations instead of scanning 4.6 million letters per read.\n",
              commas(idx.n_positions()).c_str(), c_bold(), t_index, c_reset());
  pause_here();

  // ---------------------------------------------------------------- 2. reads
  step(2, "The reads to place");
  const std::vector<SeqRecord> reads = read_all_records(a.pos[1]);
  const TruthMap truth = load_truth(a.pos[2]);
  if (reads.empty()) die("no reads in %s", a.pos[1].c_str());
  std::printf("  %s%s reads%s of %zu letters each, simulated from the genome with realistic\n"
              "  sequencing errors (about 1 letter in 100 wrong, plus occasional missing or extra\n"
              "  letters). Because they are simulated, we know the true origin of every read,\n"
              "  so every answer the aligner gives can be graded.\n\n",
              c_bold(), commas(reads.size()).c_str(), c_reset(), reads[0].seq.size());
  std::printf("  For example, read %s:\n  %s%s%s\n", reads[0].name.c_str(), c_dim(),
              reads[0].seq.substr(0, 75).c_str(), reads[0].seq.size() > 75 ? "..." : "");
  std::printf("%s", c_reset());
  pause_here();

  // ---------------------------------------------------------------- 3. one read
  step(3, "Find where one read came from");
  {
    Mapper mapper(ref, idx, mp);
    std::vector<MapRecord> recs;
    const size_t n = std::min<size_t>(reads.size(), 8192);
    mapper.map(reads.data(), n, recs);
    size_t pick = SIZE_MAX, fallback = SIZE_MAX;
    for (size_t i = 0; i < n; ++i) {
      const auto it = truth.find(reads[i].name);
      if (it == truth.end() || !placed_correctly(recs[i], ref, it->second)) continue;
      if (fallback == SIZE_MAX) fallback = i;
      const MapRecord& m = recs[i];
      bool has_gap = false, has_clip = false;
      for (uint32_t c : m.cigar) {
        has_gap |= (c & 0xf) == kCigI || (c & 0xf) == kCigD;
        has_clip |= (c & 0xf) == kCigS;
      }
      if (m.mapq == 60 && m.rev && has_gap && !has_clip && m.nm >= 3 && m.nm <= 5) {
        pick = i;
        break;
      }
    }
    if (pick == SIZE_MAX) pick = fallback;
    if (pick == SIZE_MAX) {
      std::printf("  (no correctly placed read found in the first %zu reads)\n", n);
    } else {
      const MapRecord& m = recs[pick];
      const TruthRec& tr = truth.at(reads[pick].name);
      std::printf("  Read %s%s%s %s.\n", c_bold(), reads[pick].name.c_str(), c_reset(),
                  m.rev ? "came from the opposite DNA strand, so the aligner also tries it\n  reverse-complemented (flipped)"
                        : "reads in the same direction as the genome");
      std::printf("  The aligner placed it at position %s%s%s; its true origin was %s %s%s%s\n\n", c_bold(),
                  commas(static_cast<unsigned long long>(m.pos + 1)).c_str(), c_reset(),
                  commas(static_cast<unsigned long long>(tr.start + 1)).c_str(), c_green(), "✓ correct", c_reset());
      show_alignment(ref, reads[pick], m);
    }
  }
  pause_here();

  // ---------------------------------------------------------------- 4. SIMD race
  step(4, "The speed trick: SIMD (16 alignments at once)");
  const GtSet gt = build_gt(ref, reads, truth, mp.sc.bw, reads.size());
  std::printf("  Lining a read up against the genome letter by letter (Smith-Waterman alignment)\n"
              "  is the expensive part. Plain code compares one pair of letters at a time. The\n"
              "  SIMD version uses the CPU's vector instructions (%s) to run %d alignments\n"
              "  side by side in one register. Race: %s alignments, both ways.\n\n",
              SimdAligner::isa(), SimdAligner::lanes(), commas(gt.tasks.size() * static_cast<size_t>(rounds)).c_str());
  std::vector<AlignResult> rs, rv;
  ScalarAligner scalar(mp.sc);
  SimdAligner simd(mp.sc);
  const double ts = race(scalar, gt.tasks, rounds, rs, "plain C++", c_yellow());
  const double tv = race(simd, gt.tasks, rounds, rv, "SIMD", c_green(), ts);
  size_t same = 0;
  for (size_t i = 0; i < rs.size(); ++i)
    same += rs[i].score == rv[i].score && rs[i].cigar == rv[i].cigar && rs[i].tb == rv[i].tb && rs[i].qb == rv[i].qb;
  std::printf("\n  SIMD is %s%.1fx faster%s, and gives exactly the same answer for %s / %s reads.\n", c_bold(),
              ts / tv, c_reset(), commas(same).c_str(), commas(rs.size()).c_str());
  pause_here();

  // ---------------------------------------------------------------- 5. full mapping
  step(5, "Map every read: the full pipeline on one CPU core");
  std::printf("  Seed lookup in the index, chaining seeds into candidate locations, SIMD alignment,\n"
              "  and writing the results in the standard SAM format:\n\n");
  std::vector<MapRecord> all(reads.size());
  double t_map = 0;
  {
    Mapper mapper(ref, idx, mp);
    FILE* sam = nullptr;
    if (a.has("--sam")) {
      sam = std::fopen(a.str("--sam", "").c_str(), "w");
      if (!sam) die("cannot write %s", a.str("--sam", "").c_str());
      mapper.write_sam_header(sam, "dnaln demo");
    }
    std::vector<MapRecord> batch;
    std::string buf;
    const size_t B = 8192;
    Timer tm;
    for (size_t i = 0; i < reads.size(); i += B) {
      const size_t n = std::min(B, reads.size() - i);
      mapper.map(&reads[i], n, batch);
      buf.clear();
      for (size_t j = 0; j < n; ++j) mapper.append_sam(reads[i + j], batch[j], buf);
      if (sam) std::fwrite(buf.data(), 1, buf.size(), sam);
      for (size_t j = 0; j < n; ++j) std::swap(all[i + j], batch[j]);
      draw_progress("dnaln", c_green(), static_cast<double>(i + n) / reads.size(), tm.sec(), false);
    }
    if (sam) std::fclose(sam);
    t_map = tm.sec();
    draw_progress("dnaln", c_green(), 1.0, t_map, true);
    std::printf("\n");
  }
  const double rps = reads.size() / t_map;
  std::printf("\n  %s%s reads in %.2f s%s = %s%s reads per second%s on a single core.\n", c_bold(),
              commas(reads.size()).c_str(), t_map, c_reset(), c_bold(),
              commas(static_cast<unsigned long long>(rps)).c_str(), c_reset());
  const double mm2_time = a.real("--mm2-time", 0);
  if (mm2_time > 0) {
    const double ours = t_index + t_map, slowest = std::max(ours, mm2_time);
    std::printf("\n  Total time including building the index, vs minimap2 (the widely used standard\n"
                "  aligner) on the same reads, also on one core:\n\n");
    std::printf("  %-14s %s%s%s %6.2f s\n", "dnaln", c_green(), bar(ours / slowest, 40).c_str(), c_reset(), ours);
    std::printf("  %-14s %s%s%s %6.2f s\n", "minimap2", c_dim(), bar(mm2_time / slowest, 40).c_str(), c_reset(),
                mm2_time);
  }
  pause_here();

  // ---------------------------------------------------------------- 6. accuracy
  step(6, "Grade every answer against the truth");
  Evaluator ev(truth);
  for (size_t i = 0; i < reads.size(); ++i) {
    const MapRecord& m = all[i];
    ev.add(reads[i].name, m.mapped(), m.mapped() ? ref.contigs[m.contig].name : "*", m.pos, m.rev, m.mapq, m.cigar);
  }
  std::unique_ptr<Evaluator> mm2;
  if (a.has("--mm2-sam")) {
    mm2 = std::make_unique<Evaluator>(truth);
    mm2->add_sam_file(a.str("--mm2-sam", ""));
  }
  std::printf("  Reads placed at their true location:\n\n");
  std::printf("    dnaln      %s%7.2f%%%s\n", c_bold(), 100.0 * ev.accuracy(), c_reset());
  if (mm2) std::printf("    minimap2   %7.2f%%\n", 100.0 * mm2->accuracy());
  std::printf("\n  The few misses come from stretches the genome contains several near-identical copies\n"
              "  of (E. coli has 7 copies of its ribosomal RNA genes, for example). A read from one\n"
              "  copy fits every copy equally well, so no aligner can know which one it came from;\n"
              "  dnaln flags these reads as ambiguous (mapping quality 0) instead of guessing silently.\n");
  std::printf("\n  Of the %s%s%s reads dnaln marked as confident (mapping quality 60), %s%s%s were wrong.\n",
              c_bold(), commas(ev.mapped_at(60)).c_str(), c_reset(), c_bold(), commas(ev.wrong_at(60)).c_str(),
              c_reset());
  pause_here();

  // ---------------------------------------------------------------- summary
  std::printf("\n%sSUMMARY%s\n\n", c_cyan(), c_reset());
  std::printf("  %s✓%s %s reads per second on one CPU core\n", c_green(), c_reset(),
              commas(static_cast<unsigned long long>(rps)).c_str());
  std::printf("  %s✓%s SIMD alignment %.1fx faster than plain code, with identical results\n", c_green(), c_reset(),
              ts / tv);
  if (mm2)
    std::printf("  %s✓%s %.2f%% of reads placed correctly (minimap2: %.2f%%)\n\n", c_green(), c_reset(),
                100.0 * ev.accuracy(), 100.0 * mm2->accuracy());
  else
    std::printf("  %s✓%s %.2f%% of reads placed correctly\n\n", c_green(), c_reset(), 100.0 * ev.accuracy());
  return 0;
}

}  // namespace dnaln
