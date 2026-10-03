// index.hpp — minimizer hash index of a reference genome.
//
// Layout (chosen for cache behavior):
//   * table_: open-addressing hash table, 16-byte buckets {key, offset, count}, load
//     factor <= 0.5, linear probing. Minimizer hashes are already uniformly mixed,
//     so the low bits pick the slot directly and a lookup is usually 1 cache line.
//   * pos_:   all occurrences, grouped by key, as (global_pos << 1 | strand) — the
//     hits of one minimizer are contiguous, so reading them is a sequential scan.
// Callers can prefetch() a key's bucket ahead of lookup() to overlap cache misses.
#pragma once

#include <cstdint>
#include <vector>

#include "reference.hpp"

namespace dnaln {

class MinimizerIndex {
 public:
  void build(const Reference& ref, int k, int w);

  // Returns the occurrence list of hash `h` (n = 0 if absent).
  const uint32_t* lookup(uint64_t h, uint32_t& n) const {
    size_t i = h & mask_;
    for (;;) {
      const Bucket& b = table_[i];
      if (b.key == h) {
        n = b.count;
        return pos_.data() + b.offset;
      }
      if (b.key == kEmpty) {
        n = 0;
        return nullptr;
      }
      i = (i + 1) & mask_;
    }
  }
  void prefetch(uint64_t h) const { __builtin_prefetch(&table_[h & mask_]); }

  int k() const { return k_; }
  int w() const { return w_; }
  size_t n_keys() const { return n_keys_; }
  size_t n_positions() const { return pos_.size(); }
  uint32_t max_count() const { return max_count_; }
  size_t memory_bytes() const { return table_.size() * sizeof(Bucket) + pos_.size() * sizeof(uint32_t); }

 private:
  struct Bucket {
    uint64_t key;
    uint32_t offset, count;
  };
  static constexpr uint64_t kEmpty = ~0ULL;  // hashes are < 2^(2k), so never collide with this
  std::vector<Bucket> table_;
  std::vector<uint32_t> pos_;
  uint64_t mask_ = 0;
  size_t n_keys_ = 0;
  uint32_t max_count_ = 0;
  int k_ = 0, w_ = 0;
};

}  // namespace dnaln
