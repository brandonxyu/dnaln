// mapper.hpp — the read-mapping pipeline.
//
// Reads are processed in batches so the alignment kernel always sees full SIMD
// groups:
//   1. seed:   minimizers of the read -> index lookups (prefetched) -> anchors
//   2. chain:  colinear DP chaining -> candidate loci (strand + diagonal)
//   3. rank:   reads with >1 candidate: score-only banded SW on every candidate
//   4. align:  banded SW with traceback on the best candidate -> CIGAR, NM
//   5. MAPQ from the best/second-best SW scores; SAM output
#pragma once

#include <memory>
#include <string>
#include <vector>

#include "align.hpp"
#include "chain.hpp"
#include "index.hpp"
#include "minimizer.hpp"
#include "reference.hpp"
#include "seqio.hpp"

namespace dnaln {

struct MapParams {
  int k = 15, w = 10;
  uint32_t max_occ = 500;  // ignore minimizers occurring more often than this
  ScoreParams sc;
  ChainParams ch;
  int max_cand = 10;         // candidate loci aligned per read
  double pri_ratio = 0.5;    // keep chains scoring >= pri_ratio * best chain
  int min_aln_score = 40;    // report as unmapped below this SW score
  bool simd = true;
  bool prefetch = true;
};

struct MapRecord {
  int32_t contig = -1;  // -1: unmapped
  int64_t pos = 0;      // 0-based leftmost aligned reference position within the contig
  uint8_t rev = 0, mapq = 0;
  int32_t score = 0, nm = 0, n_cand = 0;
  std::vector<uint32_t> cigar;  // including soft clips
  bool mapped() const { return contig >= 0; }
};

struct MapStats {
  uint64_t reads = 0, mapped = 0, anchors = 0, candidates = 0, cells = 0;
  double t_seed = 0, t_align = 0;
};

class Mapper {
 public:
  Mapper(const Reference& ref, const MinimizerIndex& idx, const MapParams& p);
  void map(const SeqRecord* reads, size_t n, std::vector<MapRecord>& out);
  void append_sam(const SeqRecord& r, const MapRecord& m, std::string& out) const;
  void write_sam_header(FILE* fp, const std::string& cmdline) const;
  const MapStats& stats() const { return stats_; }
  std::string kernel_name() const { return aligner_->name(); }

 private:
  struct Cand {
    uint32_t read, rev;
    int32_t contig;
    int64_t diag, wstart;  // global ref position of the read start / band window start
    int32_t sw = 0;
    const uint8_t* t = nullptr;
  };
  void seed_and_chain(const uint8_t* fwd, int len, uint32_t rid);
  AlignTask task_for(const Cand& c) const;

  const Reference& ref_;
  const MinimizerIndex& idx_;
  MapParams p_;
  std::unique_ptr<Aligner> aligner_;
  Chainer chainer_;
  MapStats stats_;

  // Per-batch scratch, reused to avoid allocation in the hot path.
  struct Hit {
    const uint32_t* pos;
    uint32_t n, mi;
  };
  std::vector<Minimizer> mins_;
  std::vector<Hit> hits_;
  std::vector<Anchor> anchors_;
  std::vector<Chain> chains_;
  std::vector<uint8_t> codes_, pad_;
  std::vector<size_t> qoff_;
  std::vector<int32_t> qlen_, chosen_, second_;
  std::vector<Cand> cands_;
  std::vector<uint32_t> cbeg_, task_ref_;
  std::vector<AlignTask> tasks_;
  std::vector<AlignResult> res_;
};

}  // namespace dnaln
