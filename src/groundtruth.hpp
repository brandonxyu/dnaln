// groundtruth.hpp — alignment tasks placed at each simulated read's true locus
// (shared by the benchmark harness and the demo).
#pragma once

#include <cstddef>
#include <vector>

#include "align.hpp"
#include "evaluate.hpp"
#include "reference.hpp"
#include "seqio.hpp"

namespace dnaln {

// One alignment per simulated read: the read in reference orientation against the
// window [true_start - bw, true_start + len + bw), i.e. the true diagonal at d = bw.
struct GtSet {
  std::vector<uint8_t> q, t;
  std::vector<size_t> qoff, toff;
  std::vector<int> qlen;
  std::vector<const TruthRec*> truth;
  std::vector<AlignTask> tasks;
};

GtSet build_gt(const Reference& ref, const std::vector<SeqRecord>& reads, const TruthMap& truth, int bw,
               size_t limit);

}  // namespace dnaln
