// minimizer.hpp — (w,k)-minimizer sketching of DNA sequences.
//
// A k-mer's canonical form is min(forward, reverse-complement), hashed with an
// invertible integer mix so lexicographically similar k-mers (poly-A etc.) do not
// dominate. In every window of w consecutive k-mers the one with the smallest hash
// is selected; consecutive windows usually share it, giving a density of ~2/(w+1).
#pragma once

#include <cstdint>
#include <vector>

namespace dnaln {

struct Minimizer {
  uint64_t hash;
  uint32_t pos;     // start of the k-mer (plus caller-supplied offset)
  uint32_t strand;  // 0: forward k-mer was canonical, 1: reverse complement was
};

// Thomas Wang's 64-bit integer hash, restricted to `mask` bits (a bijection on [0, mask]).
inline uint64_t hash64(uint64_t key, uint64_t mask) {
  key = (~key + (key << 21)) & mask;
  key = key ^ key >> 24;
  key = ((key + (key << 3)) + (key << 8)) & mask;
  key = key ^ key >> 14;
  key = ((key + (key << 2)) + (key << 4)) & mask;
  key = key ^ key >> 28;
  key = (key + (key << 31)) & mask;
  return key;
}

// Appends the minimizers of seq[0, len) to `out`. Requires 1 <= k <= 31, 1 <= w <= 255.
// N bases break k-mers and restart the window. Positions are reported + `offset`.
void sketch(const uint8_t* seq, uint32_t len, int k, int w, uint32_t offset, std::vector<Minimizer>& out);

}  // namespace dnaln
