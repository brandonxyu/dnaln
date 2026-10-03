#include "minimizer.hpp"

namespace dnaln {

void sketch(const uint8_t* seq, uint32_t len, int k, int w, uint32_t offset, std::vector<Minimizer>& out) {
  const uint64_t mask = (1ULL << (2 * k)) - 1;
  const int shift = 2 * (k - 1);
  constexpr uint64_t kInvalid = ~0ULL;  // palindromic k-mer (strand ambiguous): never selected

  // Sliding-window minimum with the van Herk / Gil-Werman scheme, streamed: k-mers are
  // grouped into blocks of w. A window ending at offset j of block B is the suffix of
  // block B-1 starting at j+1 plus the prefix of block B up to j, so its minimum is
  // min(suffix_min[B-1][j+1], prefix_min[B][j]). Prefix minima are a running min;
  // suffix minima are computed once per completed block. Everything is select-based,
  // so the hot loop has no data-dependent branches (unlike a monotone deque, whose
  // pops mispredict constantly). Ties resolve to the leftmost k-mer, as both the
  // prefix (strict <) and suffix (<=) scans and the final merge prefer earlier k-mers.
  struct Cand {
    uint64_t hash;
    uint32_t pos, strand;
  };
  Cand block[256], suffix[256];
  Cand prefix{kInvalid, 0, 0};
  Minimizer buf[256];
  int nbuf = 0;
  int64_t last = -1;  // position of the last emitted minimizer

  uint64_t fwd = 0, rev = 0;
  int run = 0;        // consecutive non-N bases
  uint32_t kidx = 0;  // k-mer index within the current N-free run
  int j = 0;          // offset of the k-mer within its block
  for (uint32_t i = 0; i < len; ++i) {
    const uint8_t c = seq[i];
    if (c > 3) {  // N: restart k-mer and window
      run = 0;
      kidx = 0;
      j = 0;
      continue;
    }
    fwd = ((fwd << 2) | c) & mask;
    rev = (rev >> 2) | (static_cast<uint64_t>(3 - c) << shift);
    if (++run < k) continue;

    const uint32_t strand = rev < fwd;
    const Cand cur{fwd != rev ? hash64(strand ? rev : fwd, mask) : kInvalid, i - k + 1, strand};
    block[j] = cur;
    prefix = (j == 0 || cur.hash < prefix.hash) ? cur : prefix;
    if (kidx + 1 >= static_cast<uint32_t>(w)) {  // a full window ends at this k-mer
      Cand m = prefix;
      if (j + 1 < w) m = prefix.hash < suffix[j + 1].hash ? prefix : suffix[j + 1];
      const bool valid = m.hash != kInvalid;
      buf[nbuf] = {m.hash, m.pos + offset, m.strand};
      nbuf += valid && static_cast<int64_t>(m.pos) != last;  // branch-free emit
      last = valid ? m.pos : last;
      if (nbuf == 256) {
        out.insert(out.end(), buf, buf + nbuf);
        nbuf = 0;
      }
    }
    if (++j == w) {  // block complete: suffix minima, right to left
      suffix[w - 1] = block[w - 1];
      for (int r = w - 2; r >= 0; --r) suffix[r] = block[r].hash <= suffix[r + 1].hash ? block[r] : suffix[r + 1];
      j = 0;
    }
    ++kidx;
  }
  out.insert(out.end(), buf, buf + nbuf);
}

}  // namespace dnaln
