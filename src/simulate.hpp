// simulate.hpp — read simulator with ground-truth alignments, and a synthetic genome generator.
#pragma once

#include <cstdint>
#include <string>

#include "reference.hpp"

namespace dnaln {

struct SimParams {
  uint64_t n = 100000;
  int len = 150;
  double sub = 0.01;       // per-base substitution rate
  double indel = 0.001;    // per-base indel rate (half insertions, half deletions)
  double indel_ext = 0.3;  // probability an indel extends by one more base
  uint64_t seed = 42;
};

// Writes <prefix>.fq and <prefix>.truth.tsv (name, contig, start, end, strand, CIGAR in
// reference orientation, 0-based half-open coordinates).
void simulate_reads(const Reference& ref, const SimParams& p, const std::string& prefix);

// Random genome with E. coli-like GC content and repeat families (rRNA-operon-like,
// IS-element-like and short REP-like copies) for offline testing.
void random_reference(uint64_t len, double gc, uint64_t seed, const std::string& path);

}  // namespace dnaln
