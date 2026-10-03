// cli.hpp — tiny command-line parser and the subcommand entry points.
#pragma once

#include <cstdlib>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "common.hpp"
#include "mapper.hpp"

namespace dnaln {

// Options take a value unless listed in `flags`; everything else is positional.
class Args {
 public:
  Args(int argc, char** argv, const std::set<std::string>& values, const std::set<std::string>& flags) {
    for (int i = 0; i < argc; ++i) {
      const std::string a = argv[i];
      if (a.size() > 1 && a[0] == '-') {
        if (flags.count(a)) {
          opts_[a] = "1";
        } else if (values.count(a)) {
          if (i + 1 >= argc) die("option %s needs a value", a.c_str());
          opts_[a] = argv[++i];
        } else {
          die("unknown option '%s' (see --help)", a.c_str());
        }
      } else {
        pos.push_back(a);
      }
    }
  }
  std::vector<std::string> pos;
  bool has(const std::string& k) const { return opts_.count(k) != 0; }
  std::string str(const std::string& k, const std::string& def) const {
    const auto it = opts_.find(k);
    return it == opts_.end() ? def : it->second;
  }
  long long num(const std::string& k, long long def) const {
    const auto it = opts_.find(k);
    if (it == opts_.end()) return def;
    char* end = nullptr;
    const long long v = std::strtoll(it->second.c_str(), &end, 10);
    if (*end) die("option %s expects an integer, got '%s'", k.c_str(), it->second.c_str());
    return v;
  }
  double real(const std::string& k, double def) const {
    const auto it = opts_.find(k);
    if (it == opts_.end()) return def;
    char* end = nullptr;
    const double v = std::strtod(it->second.c_str(), &end);
    if (*end) die("option %s expects a number, got '%s'", k.c_str(), it->second.c_str());
    return v;
  }

 private:
  std::map<std::string, std::string> opts_;
};

// Mapping/scoring options shared by `map` and `bench`.
inline const std::set<std::string> kMapValueOpts = {"-k", "-w", "--bw", "-A", "-B", "-O", "-E",
                                                    "--max-occ", "--min-score", "--max-cand", "--min-chain"};
inline const char* kMapOptsHelp =
    "  -k INT         minimizer k-mer size [15]\n"
    "  -w INT         minimizer window size [10]\n"
    "  --bw INT       band half-width of banded Smith-Waterman [16]\n"
    "  -A/-B INT      match score / mismatch penalty [2/4]\n"
    "  -O/-E INT      gap open / extension penalty; gap of length L costs O+L*E [4/2]\n"
    "  --max-occ INT  ignore minimizers with more occurrences [500]\n"
    "  --max-cand INT candidate loci aligned per read [10]\n"
    "  --min-chain INT minimum chain score [20]\n"
    "  --min-score INT minimum alignment score to report [40]\n"
    "  --scalar       use the scalar alignment kernel instead of SIMD\n"
    "  --no-prefetch  disable software prefetching of index lookups\n";

inline MapParams map_params_from_args(const Args& a) {
  MapParams p;
  p.k = static_cast<int>(a.num("-k", p.k));
  p.w = static_cast<int>(a.num("-w", p.w));
  p.sc.bw = static_cast<int>(a.num("--bw", p.sc.bw));
  p.sc.match = static_cast<int>(a.num("-A", p.sc.match));
  p.sc.mismatch = static_cast<int>(a.num("-B", p.sc.mismatch));
  p.sc.gap_open = static_cast<int>(a.num("-O", p.sc.gap_open));
  p.sc.gap_ext = static_cast<int>(a.num("-E", p.sc.gap_ext));
  p.max_occ = static_cast<uint32_t>(a.num("--max-occ", p.max_occ));
  p.max_cand = static_cast<int>(a.num("--max-cand", p.max_cand));
  p.ch.min_score = static_cast<int>(a.num("--min-chain", p.ch.min_score));
  p.min_aln_score = static_cast<int>(a.num("--min-score", p.min_aln_score));
  p.simd = !a.has("--scalar");
  p.prefetch = !a.has("--no-prefetch");
  if (p.sc.bw < 1 || p.sc.bw > 1000) die("--bw must be in [1, 1000]");
  if (p.sc.match < 1 || p.sc.mismatch < 0 || p.sc.gap_open < 0 || p.sc.gap_ext < 1)
    die("scores must satisfy A >= 1, B >= 0, O >= 0, E >= 1");
  return p;
}

int cmd_map(int argc, char** argv);
int cmd_simulate(int argc, char** argv);
int cmd_randref(int argc, char** argv);
int cmd_eval(int argc, char** argv);
int cmd_bench(int argc, char** argv);
int cmd_demo(int argc, char** argv);

}  // namespace dnaln
