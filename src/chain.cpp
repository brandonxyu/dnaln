#include "chain.hpp"

#include <algorithm>
#include <numeric>

namespace dnaln {

namespace {
inline int ilog2(uint64_t x) { return 63 - __builtin_clzll(x); }
}  // namespace

void Chainer::run(std::vector<Anchor>& a, int k, const ChainParams& p, std::vector<Chain>& out) {
  out.clear();
  const int n = static_cast<int>(a.size());
  if (n == 0) return;
  std::sort(a.begin(), a.end(), [](const Anchor& x, const Anchor& y) {
    if (x.strand != y.strand) return x.strand < y.strand;
    if (x.rpos != y.rpos) return x.rpos < y.rpos;
    return x.qpos < y.qpos;
  });

  f_.resize(n);
  prev_.resize(n);
  for (int i = 0; i < n; ++i) {
    int best = k, best_j = -1;
    const int lo = std::max(0, i - p.lookback);
    for (int j = i - 1; j >= lo; --j) {
      if (a[j].strand != a[i].strand) break;
      const int64_t dr = static_cast<int64_t>(a[i].rpos) - a[j].rpos;
      if (dr > p.max_gap) break;
      const int64_t dq = static_cast<int64_t>(a[i].qpos) - a[j].qpos;
      if (dr <= 0 || dq <= 0) continue;
      const int64_t dd = dr > dq ? dr - dq : dq - dr;
      if (dd > p.band) continue;
      const int gain = static_cast<int>(std::min<int64_t>(std::min(dq, dr), k));
      const int cost = dd ? static_cast<int>(0.01 * k * static_cast<double>(dd) + 0.5 * ilog2(dd)) : 0;
      const int sc = f_[j] + gain - cost;
      if (sc > best) {
        best = sc;
        best_j = j;
      }
    }
    f_[i] = best;
    prev_[i] = best_j;
  }

  // Extract disjoint chains, best end-points first (as in minimap2's mg_chain_backtrack).
  order_.resize(n);
  std::iota(order_.begin(), order_.end(), 0);
  std::sort(order_.begin(), order_.end(), [&](int x, int y) { return f_[x] != f_[y] ? f_[x] > f_[y] : x < y; });
  used_.assign(n, 0);
  for (int end : order_) {
    if (used_[end]) continue;
    int i = end, cnt = 0;
    int64_t dmin = INT64_MAX, dmax = INT64_MIN;
    while (i >= 0 && !used_[i]) {
      used_[i] = 1;
      ++cnt;
      const int64_t d = static_cast<int64_t>(a[i].rpos) - a[i].qpos;
      dmin = std::min(dmin, d);
      dmax = std::max(dmax, d);
      i = prev_[i];
    }
    const int score = f_[end] - (i >= 0 ? f_[i] : 0);
    if (score < p.min_score || cnt < p.min_anchors) continue;
    out.push_back({score, a[end].strand, cnt, dmin, dmax, a[end].rpos});
  }
  std::stable_sort(out.begin(), out.end(), [](const Chain& x, const Chain& y) { return x.score > y.score; });
}

}  // namespace dnaln
