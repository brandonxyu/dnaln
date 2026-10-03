// align_simd.cpp — inter-sequence vectorized banded Smith–Waterman.
//
// kLanes alignments run in lock-step, one per 16-bit lane. Sequences are first
// transposed into lane-interleaved byte rows (Q[i*L + lane], T[j*L + lane]) so every
// cell update is a handful of unit-stride vector loads/stores, and substitution
// scores come from one byte-table lookup per 16 lanes (see simd.hpp).
//
// Memory/cache design:
//   * Only one band row of H and E is kept ((W+1) x L int16 each — ~2 KB for bw=16),
//     updated in place; it lives in L1 for the whole group.
//   * Sequences are 1 byte per base per lane (~5 KB per group for 150 bp reads).
//   * The traceback matrix is 1 byte per cell (direction bits packed with NEON
//     shift-insert), i.e. qlen*W*L bytes ≈ 79 KB for 150 bp reads, bw=16, L=16 —
//     sized to stay within L1, and written strictly sequentially.
//   * Scratch buffers are 64-byte aligned and reused across groups (no allocation
//     in the hot loop).
#include <algorithm>
#include <cstring>

#include "align.hpp"
#include "simd.hpp"

namespace dnaln {

using namespace simd;

SimdAligner::SimdAligner(const ScoreParams& sp) : sp_(sp), fallback_(sp) {
  for (int a = 0; a < 4; ++a)
    for (int b = 0; b < 4; ++b) sub16_[a << 2 | b] = static_cast<int8_t>(sp.sub(static_cast<uint8_t>(a), static_cast<uint8_t>(b)));
  // Byte-table scoring needs every score (and score + ambig) to fit in int8.
  byte_scores_ok_ = sp.match + sp.ambig <= 127 && sp.mismatch <= 127 && sp.ambig <= 127;
}

int SimdAligner::lanes() { return kLanes; }
const char* SimdAligner::isa() { return kIsa; }

std::string SimdAligner::name() const {
  return std::string(kIsa) + " x" + std::to_string(kLanes) + " lanes";
}

size_t SimdAligner::working_set_bytes(int qlen, bool traceback) const {
  const size_t W = sp_.band_cols(), L = kLanes;
  size_t b = 2 * (W + 1) * L * sizeof(int16_t)  // H, E band rows
             + W * L * sizeof(int16_t)              // column-index constants
             + (2 * qlen + W - 1) * L;              // transposed query + target bytes
  if (traceback) b += qlen * W * L;                          // direction bytes
  return b;
}

void SimdAligner::align(const AlignTask* tasks, size_t n, bool traceback, AlignResult* out) {
  // Keep 16-bit scores far from saturation; longer queries use the scalar kernel.
  const int max_len = 30000 / std::max(1, sp_.match);
  uint32_t idx[kLanes];
  int m = 0;
  auto flush = [&] {
    if (m == 0) return;
    if (traceback) align_group<true>(tasks, idx, m, out);
    else align_group<false>(tasks, idx, m, out);
    m = 0;
  };
  for (size_t i = 0; i < n; ++i) {
    if (tasks[i].qlen > max_len || !byte_scores_ok_) {
      fallback_.align(&tasks[i], 1, traceback, &out[i]);
      continue;
    }
    idx[m++] = static_cast<uint32_t>(i);
    if (m == kLanes) flush();
  }
  flush();
}

template <bool TB>
void SimdAligner::align_group(const AlignTask* tasks, const uint32_t* idx, int n, AlignResult* out) {
  constexpr int L = kLanes;
  const int W = sp_.band_cols();
  int qmax = 0;
  for (int l = 0; l < n; ++l) qmax = std::max(qmax, tasks[idx[l]].qlen);
  const int tmax = qmax + W - 1;

  uint8_t* Q = q_.get(static_cast<size_t>(qmax) * L);
  uint8_t* T = t_.get(static_cast<size_t>(tmax) * L);
  int16_t* H = h_.get(static_cast<size_t>(W + 1) * L);
  int16_t* E = e_.get(static_cast<size_t>(W + 1) * L);
  int16_t* Dc = TB ? dconst_.get(static_cast<size_t>(W) * L) : nullptr;
  uint8_t* D = TB ? dir_.get(static_cast<size_t>(qmax) * W * L) : nullptr;

  // Transpose into lane-interleaved byte rows (see simd.hpp for the encoding). Unused
  // lanes and tails are N, which scores negatively and can never raise a lane's best.
  for (int l = 0; l < L; ++l) {
    const AlignTask* tk = l < n ? &tasks[idx[l]] : nullptr;
    const int ql = tk ? tk->qlen : 0, tl = tk ? std::min(tk->tlen, tmax) : 0;
    for (int i = 0; i < ql; ++i) {
      const uint8_t c = tk->q[i];
      Q[static_cast<size_t>(i) * L + l] = c < 4 ? static_cast<uint8_t>(c << 2) : kSeqN;
    }
    for (int i = ql; i < qmax; ++i) Q[static_cast<size_t>(i) * L + l] = kSeqN;
    for (int j = 0; j < tl; ++j) {
      const uint8_t c = tk->t[j];
      T[static_cast<size_t>(j) * L + l] = c < 4 ? c : kSeqN;
    }
    for (int j = tl; j < tmax; ++j) T[static_cast<size_t>(j) * L + l] = kSeqN;
  }
  std::memset(H, 0, sizeof(int16_t) * (W + 1) * L);
  std::memset(E, 0, sizeof(int16_t) * (W + 1) * L);
  if (TB)
    for (int d = 0; d < W; ++d) std::fill(Dc + static_cast<size_t>(d) * L, Dc + static_cast<size_t>(d + 1) * L, static_cast<int16_t>(d));

  const SubTable st = make_subtable(sub16_, static_cast<int8_t>(-sp_.ambig));
  const Vec zero = vset1(0);
  const Vec vOE = vset1(static_cast<int16_t>(sp_.gap_open + sp_.gap_ext));
  const Vec vGE = vset1(static_cast<int16_t>(sp_.gap_ext));
  Vec best = zero, best_i = vset1(-1), best_d = vset1(-1);

  for (int i = 0; i < qmax; ++i) {
    const QueryRow qrow = qload(Q + static_cast<size_t>(i) * L);
    const uint8_t* Ti = T + static_cast<size_t>(i) * L;
    uint8_t* Di = TB ? D + static_cast<size_t>(i) * W * L : nullptr;
    const Vec row_start_best = best;
    Vec f = zero, hleft = zero, hdiag = vload(H);
    for (int d = 0; d < W; ++d) {
      const Vec hup = vload(H + static_cast<size_t>(d + 1) * L);
      const Vec eup = vload(E + static_cast<size_t>(d + 1) * L);
      const Vec hd = vadds(hdiag, vscore(st, qrow, Ti + static_cast<size_t>(d) * L));
      const Vec eo = vsubs(hup, vOE), ee = vsubs(eup, vGE);
      const Vec e = vmax(eo, ee);
      const Vec fo = vsubs(hleft, vOE), fe = vsubs(f, vGE);
      f = vmax(fo, fe);
      // F last: it is the only term on the loop-carried critical path.
      const Vec h = vmax(vmax(vmax(hd, e), zero), f);
      vstore(H + static_cast<size_t>(d) * L, h);
      vstore(E + static_cast<size_t>(d) * L, e);
      if constexpr (TB) {
        vstore_masks8(Di + static_cast<size_t>(d) * L, veq(h, hd), veq(h, e), vgt(ee, eo), vgt(fe, fo));
        const Vec up = vgt(h, best);  // strict '>' in row-major order == scalar tie-breaking
        best = vmax(best, h);
        best_d = vblend(up, vload(Dc + static_cast<size_t>(d) * L), best_d);
      } else {
        best = vmax(best, h);
      }
      hleft = h;
      hdiag = hup;
    }
    if constexpr (TB) best_i = vblend(vgt(best, row_start_best), vset1(static_cast<int16_t>(i)), best_i);
  }

  alignas(64) int16_t sb[L], si[L], sd[L];
  vstore(sb, best);
  vstore(si, best_i);
  vstore(sd, best_d);
  for (int l = 0; l < n; ++l) {
    AlignResult& r = out[idx[l]];
    r.clear();
    r.score = sb[l];
    if (TB && sb[l] > 0) traceback(D, L, l, W, si[l], sd[l], sb[l], tasks[idx[l]].q, tasks[idx[l]].t, sp_, r);
  }
}

}  // namespace dnaln
