// align.cpp — scalar baseline kernel and the traceback shared by both kernels.
#include <algorithm>

#include "align.hpp"

// The scalar baseline must stay scalar: forbid the compiler from vectorizing it so the
// SIMD speedup is measured against genuinely scalar code. (The F dependency along the
// row already prevents loop vectorization; this just makes the intent explicit.)
#if defined(__clang__)
#define DNALN_NO_VECTORIZE _Pragma("clang loop vectorize(disable) interleave(disable)")
#elif defined(__GNUC__) && __GNUC__ >= 14
#define DNALN_NO_VECTORIZE _Pragma("GCC novector")
#else
#define DNALN_NO_VECTORIZE
#endif

namespace dnaln {

ScalarAligner::ScalarAligner(const ScoreParams& sp) : sp_(sp) {
  for (int a = 0; a < 5; ++a)
    for (int b = 0; b < 8; ++b) mat_[a * 8 + b] = sp.sub(static_cast<uint8_t>(a), static_cast<uint8_t>(b));
}

void ScalarAligner::align(const AlignTask* tasks, size_t n, bool traceback, AlignResult* out) {
  for (size_t i = 0; i < n; ++i) {
    if (traceback) align_one<true>(tasks[i], out[i]);
    else align_one<false>(tasks[i], out[i]);
  }
}

template <bool TB>
void ScalarAligner::align_one(const AlignTask& tk, AlignResult& r) {
  const int W = sp_.band_cols(), ql = tk.qlen, tmax = ql + W - 1;
  const int oe = sp_.gap_open + sp_.gap_ext, ge = sp_.gap_ext;

  uint8_t* t = t_.get(static_cast<size_t>(tmax));
  const int tl = std::min(tk.tlen, tmax);
  std::copy(tk.t, tk.t + tl, t);
  std::fill(t + tl, t + tmax, uint8_t{4});

  int32_t* H = h_.get(static_cast<size_t>(W) + 1);
  int32_t* E = e_.get(static_cast<size_t>(W) + 1);
  std::fill(H, H + W + 1, 0);  // row -1 and the out-of-band cell H[W] are 0
  std::fill(E, E + W + 1, 0);
  uint8_t* D = TB ? dir_.get(static_cast<size_t>(ql) * W) : nullptr;

  int best = 0, best_i = -1, best_d = -1;
  for (int i = 0; i < ql; ++i) {
    const int* sc = mat_ + tk.q[i] * 8;
    const uint8_t* ti = t + i;
    uint8_t* Di = TB ? D + static_cast<size_t>(i) * W : nullptr;
    int f = 0, hleft = 0, hdiag = H[0];
    DNALN_NO_VECTORIZE
    for (int d = 0; d < W; ++d) {
      const int hup = H[d + 1], eup = E[d + 1];
      const int hd = hdiag + sc[ti[d]];
      const int eo = hup - oe, ee = eup - ge;
      const int e = eo > ee ? eo : ee;
      const int fo = hleft - oe, fe = f - ge;
      f = fo > fe ? fo : fe;
      int h = hd > e ? hd : e;
      h = h > 0 ? h : 0;
      h = h > f ? h : f;
      H[d] = h;
      E[d] = e;
      if constexpr (TB) {
        Di[d] = static_cast<uint8_t>((h == hd ? kDirDiag : 0) | (h == e ? kDirE : 0) | (ee > eo ? kDirExtE : 0) |
                                     (fe > fo ? 0x1F : 0));
        if (h > best) {
          best = h;
          best_i = i;
          best_d = d;
        }
      } else {
        best = h > best ? h : best;
      }
      hleft = h;
      hdiag = hup;
    }
  }
  r.clear();
  r.score = best;
  if (TB && best > 0) traceback(D, 1, 0, W, best_i, best_d, best, tk.q, tk.t, sp_, r);
}


void traceback(const uint8_t* dir, int stride, int lane, int W, int i_end, int d_end, int score,
               const uint8_t* q, const uint8_t* t, const ScoreParams& sp, AlignResult& r) {
  std::vector<uint32_t>& c = r.cigar;
  c.clear();
  auto push = [&c](uint32_t op, int len) {
    if (!c.empty() && (c.back() & 0xf) == op) c.back() += static_cast<uint32_t>(len) << 4;
    else c.push_back(static_cast<uint32_t>(len) << 4 | op);
  };
  auto at = [=](int i, int d) { return dir[(static_cast<size_t>(i) * W + d) * stride + lane]; };
  const int oe = sp.gap_open + sp.gap_ext, ge = sp.gap_ext;
  // v is the DP value of the current cell on the optimal path: H in state H, E/F in the
  // gap states. Every step inverts the recurrence, so v hits exactly 0 where the local
  // alignment begins.
  int i = i_end, d = d_end, v = score, nm = 0;
  for (;;) {
    int run = 0;  // state H: follow the diagonal
    uint8_t x = 0;
    while (i >= 0 && v > 0 && ((x = at(i, d)) & kDirDiag)) {
      const uint8_t a = q[i], b = t[i + d];
      v -= sp.sub(a, b);
      nm += a != b || a > 3;
      ++run;
      --i;
    }
    if (run) push(kCigM, run);
    if (i < 0 || v <= 0) break;
    run = 0;
    if (x & kDirE) {  // insertion: E(i,d) <- (i-1, d+1)
      bool ext;
      do {
        ext = at(i, d) & kDirExtE;
        v += ext ? ge : oe;
        ++run;
        --i;
        ++d;
      } while (ext);
      push(kCigI, run);
    } else {  // deletion: F(i,d) <- (i, d-1)
      bool ext;
      do {
        ext = at(i, d) & kDirExtF;
        v += ext ? ge : oe;
        ++run;
        --d;
      } while (ext);
      push(kCigD, run);
    }
    nm += run;
  }
  std::reverse(c.begin(), c.end());
  // The walk always ends with a diagonal step, so (i+1, d) is the first aligned cell.
  r.qb = i + 1;
  r.tb = i + 1 + d;
  r.qe = i_end + 1;
  r.te = i_end + d_end + 1;
  r.nm = nm;
}

int score_cigar(const uint8_t* q, const uint8_t* t, int qb, int tb, const std::vector<uint32_t>& cigar,
                const ScoreParams& sp) {
  int s = 0, i = qb, j = tb;
  for (uint32_t c : cigar) {
    const uint32_t op = c & 0xf;
    const int len = static_cast<int>(c >> 4);
    if (op == kCigM || op == kCigEq || op == kCigX) {
      for (int k = 0; k < len; ++k) s += sp.sub(q[i++], t[j++]);
    } else if (op == kCigI) {
      s -= sp.gap_open + len * sp.gap_ext;
      i += len;
    } else if (op == kCigD) {
      s -= sp.gap_open + len * sp.gap_ext;
      j += len;
    }
  }
  return s;
}

}  // namespace dnaln
