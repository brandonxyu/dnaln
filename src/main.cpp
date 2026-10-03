// main.cpp — `dnaln` command-line entry point.
#include <cstdio>
#include <cstring>
#include <string>

#include "cli.hpp"
#include "evaluate.hpp"
#include "simulate.hpp"

namespace dnaln {

int cmd_map(int argc, char** argv) {
  std::set<std::string> values = kMapValueOpts;
  values.insert({"-o", "--batch"});
  const Args a(argc, argv, values, {"--scalar", "--no-prefetch", "-h", "--help"});
  if (a.pos.size() != 2 || a.has("-h") || a.has("--help")) {
    std::fprintf(stderr,
                 "Usage: dnaln map [options] <ref.fa[.gz]> <reads.fq[.gz]>  > out.sam\n\n"
                 "  -o FILE        write SAM to FILE instead of stdout\n"
                 "  --batch INT    reads per batch [8192]\n%s",
                 kMapOptsHelp);
    return a.has("-h") || a.has("--help") ? 0 : 1;
  }
  const MapParams p = map_params_from_args(a);
  const size_t batch = static_cast<size_t>(a.num("--batch", 8192));
  std::string cmdline = "dnaln map";
  for (int i = 0; i < argc; ++i) cmdline += std::string(" ") + argv[i];

  Timer t_total;
  Timer t;
  const Reference ref = Reference::load(a.pos[0]);
  const double t_load = t.sec();
  t.reset();
  MinimizerIndex idx;
  idx.build(ref, p.k, p.w);
  const double t_index = t.sec();
  log_msg("reference: %zu contig(s), %zu bp (loaded in %.2f s)", ref.contigs.size(), ref.seq.size(), t_load);
  log_msg("index: k=%d w=%d, %zu distinct minimizers, %zu positions, %.1f MB, built in %.2f s", p.k, p.w,
          idx.n_keys(), idx.n_positions(), idx.memory_bytes() / 1048576.0, t_index);

  FILE* out = stdout;
  if (a.has("-o")) {
    out = std::fopen(a.str("-o", "").c_str(), "w");
    if (!out) die("cannot open '%s' for writing", a.str("-o", "").c_str());
  }
  Mapper mapper(ref, idx, p);
  mapper.write_sam_header(out, cmdline);
  log_msg("alignment kernel: %s", mapper.kernel_name().c_str());

  SeqReader reader(a.pos[1]);
  std::vector<SeqRecord> reads;
  std::vector<MapRecord> recs;
  std::string buf;
  Timer t_map;
  size_t n;
  while ((n = reader.next_batch(reads, batch)) > 0) {
    mapper.map(reads.data(), n, recs);
    buf.clear();
    for (size_t i = 0; i < n; ++i) mapper.append_sam(reads[i], recs[i], buf);
    std::fwrite(buf.data(), 1, buf.size(), out);
  }
  if (out != stdout) std::fclose(out);
  else std::fflush(out);
  const double sec = t_map.sec();
  const MapStats& st = mapper.stats();
  log_msg("mapped %llu / %llu reads (%.2f%%) in %.3f s -> %.0f reads/s single-threaded (excl. index build)",
          static_cast<unsigned long long>(st.mapped), static_cast<unsigned long long>(st.reads),
          st.reads ? 100.0 * st.mapped / st.reads : 0.0, sec, st.reads / sec);
  log_msg("time split: seed+chain %.3f s, alignment %.3f s (%.2f GCUPS), parse+output %.3f s; total wall %.2f s",
          st.t_seed, st.t_align, st.t_align > 0 ? st.cells / st.t_align * 1e-9 : 0.0,
          sec - st.t_seed - st.t_align, t_total.sec());
  return 0;
}

int cmd_simulate(int argc, char** argv) {
  const Args a(argc, argv, {"-n", "-l", "-e", "-i", "--indel-ext", "-s", "-o"}, {"-h", "--help"});
  if (a.pos.size() != 1 || a.has("-h") || a.has("--help")) {
    std::fprintf(stderr,
                 "Usage: dnaln simulate [options] <ref.fa>\n\n"
                 "Simulates reads with known true alignments. Writes <prefix>.fq and <prefix>.truth.tsv.\n"
                 "  -n INT        number of reads [100000]\n"
                 "  -l INT        read length [150]\n"
                 "  -e FLOAT      substitution rate per base [0.01]\n"
                 "  -i FLOAT      indel rate per base [0.001]\n"
                 "  --indel-ext F indel extension probability [0.3]\n"
                 "  -s INT        random seed [42]\n"
                 "  -o PREFIX     output prefix [sim]\n");
    return a.has("-h") || a.has("--help") ? 0 : 1;
  }
  SimParams p;
  p.n = static_cast<uint64_t>(a.num("-n", static_cast<long long>(p.n)));
  p.len = static_cast<int>(a.num("-l", p.len));
  p.sub = a.real("-e", p.sub);
  p.indel = a.real("-i", p.indel);
  p.indel_ext = a.real("--indel-ext", p.indel_ext);
  p.seed = static_cast<uint64_t>(a.num("-s", static_cast<long long>(p.seed)));
  const std::string prefix = a.str("-o", "sim");
  const Reference ref = Reference::load(a.pos[0]);
  Timer t;
  simulate_reads(ref, p, prefix);
  log_msg("simulated %llu reads of %d bp (sub %.4f, indel %.4f) in %.2f s -> %s.fq, %s.truth.tsv",
          static_cast<unsigned long long>(p.n), p.len, p.sub, p.indel, t.sec(), prefix.c_str(), prefix.c_str());
  return 0;
}

int cmd_randref(int argc, char** argv) {
  const Args a(argc, argv, {"-L", "--gc", "-s", "-o"}, {"-h", "--help"});
  if (!a.pos.empty() || a.has("-h") || a.has("--help")) {
    std::fprintf(stderr,
                 "Usage: dnaln randref [-L 4641652] [--gc 0.508] [-s 1] -o synthetic.fa\n\n"
                 "Writes a random genome with E. coli-like size, GC content and repeat families\n"
                 "(for offline testing when the real reference is not available).\n");
    return a.has("-h") || a.has("--help") ? 0 : 1;
  }
  const uint64_t len = static_cast<uint64_t>(a.num("-L", 4641652));
  const std::string out = a.str("-o", "synthetic.fa");
  random_reference(len, a.real("--gc", 0.508), static_cast<uint64_t>(a.num("-s", 1)), out);
  log_msg("wrote %llu bp synthetic genome to %s", static_cast<unsigned long long>(len), out.c_str());
  return 0;
}

int cmd_eval(int argc, char** argv) {
  const Args a(argc, argv, {"--tol", "--max-gap"}, {"-h", "--help"});
  if (a.pos.size() < 2 || a.has("-h") || a.has("--help")) {
    std::fprintf(stderr,
                 "Usage: dnaln eval [--tol 10] [--max-gap 2.0] <truth.tsv> <a.sam> [b.sam ...]\n\n"
                 "Scores SAM files against the simulator's ground truth, side by side.\n"
                 "With several SAMs, the last one is the reference for the accuracy-gap check\n"
                 "(e.g. `dnaln eval truth.tsv dnaln.sam minimap2.sam`).\n");
    return a.has("-h") || a.has("--help") ? 0 : 1;
  }
  const TruthMap truth = load_truth(a.pos[0]);
  std::vector<Evaluator> evs;
  evs.reserve(a.pos.size() - 1);
  std::vector<std::string> labels;
  std::vector<const Evaluator*> ptrs;
  for (size_t i = 1; i < a.pos.size(); ++i) {
    evs.emplace_back(truth, static_cast<int>(a.num("--tol", 10)));
    evs.back().add_sam_file(a.pos[i]);
    const std::string& path = a.pos[i];
    const size_t slash = path.find_last_of('/');
    labels.push_back(slash == std::string::npos ? path : path.substr(slash + 1));
  }
  for (const auto& e : evs) ptrs.push_back(&e);
  print_eval_report(stdout, labels, ptrs, a.real("--max-gap", 2.0));
  return 0;
}

}  // namespace dnaln

namespace {
void usage() {
  std::fprintf(stderr,
               "dnaln %s — fast DNA short-read aligner (minimizer seeds + SIMD banded Smith-Waterman)\n\n"
               "Usage: dnaln map [options] <genome.fa> <reads.fq> > alignments.sam\n"
               "       dnaln <command> [options]\n\n"
               "Commands:\n"
               "  map       align reads to a reference, write SAM\n"
               "  simulate  simulate reads with ground-truth alignments\n"
               "  randref   generate a synthetic E. coli-like genome (offline testing)\n"
               "  eval      score SAM file(s) against ground truth (e.g. dnaln vs minimap2)\n"
               "  bench     benchmark harness: kernel speedup, correctness, throughput, cache sweep\n"
               "  demo      narrated walkthrough for presentations (run ./demo.sh)\n\n"
               "Run `dnaln <command> --help` for options, `dnaln --version` for the version.\n",
               dnaln::kVersion);
}
}  // namespace

int main(int argc, char** argv) {
  if (argc < 2) {
    usage();
    return 1;
  }
  const std::string cmd = argv[1];
  if (cmd == "map") return dnaln::cmd_map(argc - 2, argv + 2);
  if (cmd == "simulate") return dnaln::cmd_simulate(argc - 2, argv + 2);
  if (cmd == "randref") return dnaln::cmd_randref(argc - 2, argv + 2);
  if (cmd == "eval") return dnaln::cmd_eval(argc - 2, argv + 2);
  if (cmd == "bench") return dnaln::cmd_bench(argc - 2, argv + 2);
  if (cmd == "demo") return dnaln::cmd_demo(argc - 2, argv + 2);
  if (cmd == "-v" || cmd == "--version" || cmd == "version") {
    std::printf("dnaln %s (%s SIMD)\n", dnaln::kVersion, dnaln::SimdAligner::isa());
    return 0;
  }
  if (cmd == "-h" || cmd == "--help" || cmd == "help") {
    usage();
    return 0;
  }
  std::fprintf(stderr, "unknown command '%s'\n\n", cmd.c_str());
  usage();
  return 1;
}
