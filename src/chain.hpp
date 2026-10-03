// chain.hpp — colinear chaining of minimizer anchors (minimap2-style DP).
//
// An anchor is a shared minimizer between read and reference. Anchors are put in
// "strand-oriented" read coordinates (positions on the reverse-complemented read
// for reverse-strand hits) so that a true hit is a set of anchors increasing in
// both coordinates. The DP links anchors i <- j with
//     f(i) = max(k, max_j f(j) + min(dq, dr, k) - gap_cost(|dr - dq|))
// over the previous `lookback` anchors, then extracts the best disjoint chains.
#pragma once

#include <cstdint>
#include <vector>

namespace dnaln {

struct Anchor {
  uint32_t rpos;    // global reference start of the k-mer
  int32_t qpos;     // start of the k-mer on the strand-oriented read
  uint32_t strand;  // 1 = read maps to the reverse strand
};

struct Chain {
  int32_t score;
  uint32_t strand;
  int32_t n_anchors;
  int64_t dmin, dmax;  // range of diagonals (rpos - qpos) covered by the chain
  uint32_t rpos;       // reference position of the chain's last anchor
};

struct ChainParams {
  int max_gap = 5000;  // max reference distance between linked anchors
  int band = 16;       // max diagonal drift (indel) between linked anchors
  int lookback = 50;   // predecessors examined per anchor
  int min_score = 20;
  int min_anchors = 2;
};

class Chainer {
 public:
  // Sorts `a` in place; writes chains sorted by descending score.
  void run(std::vector<Anchor>& a, int k, const ChainParams& p, std::vector<Chain>& out);

 private:
  std::vector<int32_t> f_, prev_, order_;
  std::vector<uint8_t> used_;
};

}  // namespace dnaln
