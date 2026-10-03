#include "index.hpp"

#include <algorithm>

#include "common.hpp"
#include "minimizer.hpp"

namespace dnaln {

void MinimizerIndex::build(const Reference& ref, int k, int w) {
  if (k < 4 || k > 28) die("k must be in [4, 28] (got %d)", k);
  if (w < 1 || w > 255) die("w must be in [1, 255] (got %d)", w);
  k_ = k;
  w_ = w;

  std::vector<Minimizer> mins;
  mins.reserve(ref.seq.size() * 2 / (w + 1) + 1024);
  for (const Contig& c : ref.contigs)
    sketch(ref.seq.data() + c.offset, c.len, k, w, static_cast<uint32_t>(c.offset), mins);

  std::sort(mins.begin(), mins.end(), [](const Minimizer& a, const Minimizer& b) {
    return a.hash != b.hash ? a.hash < b.hash : a.pos < b.pos;
  });

  n_keys_ = 0;
  for (size_t i = 0; i < mins.size(); ++i)
    if (i == 0 || mins[i].hash != mins[i - 1].hash) ++n_keys_;

  size_t cap = 16;
  while (cap < n_keys_ * 2) cap <<= 1;
  table_.assign(cap, Bucket{kEmpty, 0, 0});
  mask_ = cap - 1;
  pos_.resize(mins.size());
  max_count_ = 0;

  for (size_t i = 0; i < mins.size();) {
    size_t j = i;
    for (; j < mins.size() && mins[j].hash == mins[i].hash; ++j) pos_[j] = mins[j].pos << 1 | mins[j].strand;
    size_t s = mins[i].hash & mask_;
    while (table_[s].key != kEmpty) s = (s + 1) & mask_;
    table_[s] = {mins[i].hash, static_cast<uint32_t>(i), static_cast<uint32_t>(j - i)};
    max_count_ = std::max(max_count_, static_cast<uint32_t>(j - i));
    i = j;
  }
}

}  // namespace dnaln
