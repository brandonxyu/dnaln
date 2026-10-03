// simd.hpp — thin portability layer over 16-bit-lane integer SIMD.
//
//   AVX2   : __m256i   16 lanes / register
//   SSE4.1 : __m128i    8 lanes / register
//   NEON   : int16x8_t  8 lanes / register   (Apple Silicon, ARM64 Linux)
//   other  : plain arrays (keeps the code portable; the compiler may auto-vectorize)
//
// The alignment kernel operates on `Vec`, which bundles kUnroll native registers so
// each logical operation issues kUnroll independent instructions. The banded
// recurrence has a loop-carried dependency (H -> F -> H along the band row); having
// two independent register chains in flight hides that latency on wide out-of-order
// cores (e.g. Apple M-series with 4 NEON pipes).
#pragma once

#include <cstdint>

#if defined(DNALN_NO_SIMD)  // force the portable path (testing)
#define DNALN_ISA_PORTABLE 1
#elif defined(__AVX2__)
#include <immintrin.h>
#define DNALN_ISA_AVX2 1
#elif defined(__SSE4_1__)
#include <smmintrin.h>
#define DNALN_ISA_SSE41 1
#elif defined(__ARM_NEON) || defined(__ARM_NEON__)
#include <arm_neon.h>
#define DNALN_ISA_NEON 1
#else
#define DNALN_ISA_PORTABLE 1
#endif

#ifndef DNALN_SIMD_UNROLL
#if defined(DNALN_ISA_AVX2)
#define DNALN_SIMD_UNROLL 1
#else
#define DNALN_SIMD_UNROLL 2
#endif
#endif

#define DNALN_INLINE inline __attribute__((always_inline))

namespace dnaln::simd {

#if defined(DNALN_ISA_AVX2)
inline constexpr const char* kIsa = "AVX2";
inline constexpr int kNative = 16;
using Reg = __m256i;
DNALN_INLINE Reg r_set1(int16_t x) { return _mm256_set1_epi16(x); }
DNALN_INLINE Reg r_load(const int16_t* p) { return _mm256_loadu_si256(reinterpret_cast<const __m256i*>(p)); }
DNALN_INLINE void r_store(int16_t* p, Reg a) { _mm256_storeu_si256(reinterpret_cast<__m256i*>(p), a); }
DNALN_INLINE Reg r_adds(Reg a, Reg b) { return _mm256_adds_epi16(a, b); }
DNALN_INLINE Reg r_subs(Reg a, Reg b) { return _mm256_subs_epi16(a, b); }
DNALN_INLINE Reg r_max(Reg a, Reg b) { return _mm256_max_epi16(a, b); }
DNALN_INLINE Reg r_eq(Reg a, Reg b) { return _mm256_cmpeq_epi16(a, b); }
DNALN_INLINE Reg r_gt(Reg a, Reg b) { return _mm256_cmpgt_epi16(a, b); }
DNALN_INLINE Reg r_and(Reg a, Reg b) { return _mm256_and_si256(a, b); }
DNALN_INLINE Reg r_or(Reg a, Reg b) { return _mm256_or_si256(a, b); }
DNALN_INLINE Reg r_blend(Reg m, Reg a, Reg b) { return _mm256_blendv_epi8(b, a, m); }  // m ? a : b
DNALN_INLINE void r_store_lo8(uint8_t* p, Reg a) {  // low byte of each lane (values 0..255)
  const __m128i lo = _mm256_castsi256_si128(a), hi = _mm256_extracti128_si256(a, 1);
  _mm_storeu_si128(reinterpret_cast<__m128i*>(p), _mm_packus_epi16(lo, hi));
}
#elif defined(DNALN_ISA_SSE41)
inline constexpr const char* kIsa = "SSE4.1";
inline constexpr int kNative = 8;
using Reg = __m128i;
DNALN_INLINE Reg r_set1(int16_t x) { return _mm_set1_epi16(x); }
DNALN_INLINE Reg r_load(const int16_t* p) { return _mm_loadu_si128(reinterpret_cast<const __m128i*>(p)); }
DNALN_INLINE void r_store(int16_t* p, Reg a) { _mm_storeu_si128(reinterpret_cast<__m128i*>(p), a); }
DNALN_INLINE Reg r_adds(Reg a, Reg b) { return _mm_adds_epi16(a, b); }
DNALN_INLINE Reg r_subs(Reg a, Reg b) { return _mm_subs_epi16(a, b); }
DNALN_INLINE Reg r_max(Reg a, Reg b) { return _mm_max_epi16(a, b); }
DNALN_INLINE Reg r_eq(Reg a, Reg b) { return _mm_cmpeq_epi16(a, b); }
DNALN_INLINE Reg r_gt(Reg a, Reg b) { return _mm_cmpgt_epi16(a, b); }
DNALN_INLINE Reg r_and(Reg a, Reg b) { return _mm_and_si128(a, b); }
DNALN_INLINE Reg r_or(Reg a, Reg b) { return _mm_or_si128(a, b); }
DNALN_INLINE Reg r_blend(Reg m, Reg a, Reg b) { return _mm_blendv_epi8(b, a, m); }
DNALN_INLINE void r_store_lo8(uint8_t* p, Reg a) {
  _mm_storel_epi64(reinterpret_cast<__m128i*>(p), _mm_packus_epi16(a, a));
}
#elif defined(DNALN_ISA_NEON)
inline constexpr const char* kIsa = "NEON";
inline constexpr int kNative = 8;
using Reg = int16x8_t;
DNALN_INLINE Reg r_set1(int16_t x) { return vdupq_n_s16(x); }
DNALN_INLINE Reg r_load(const int16_t* p) { return vld1q_s16(p); }
DNALN_INLINE void r_store(int16_t* p, Reg a) { vst1q_s16(p, a); }
DNALN_INLINE Reg r_adds(Reg a, Reg b) { return vqaddq_s16(a, b); }
DNALN_INLINE Reg r_subs(Reg a, Reg b) { return vqsubq_s16(a, b); }
DNALN_INLINE Reg r_max(Reg a, Reg b) { return vmaxq_s16(a, b); }
DNALN_INLINE Reg r_eq(Reg a, Reg b) { return vreinterpretq_s16_u16(vceqq_s16(a, b)); }
DNALN_INLINE Reg r_gt(Reg a, Reg b) { return vreinterpretq_s16_u16(vcgtq_s16(a, b)); }
DNALN_INLINE Reg r_and(Reg a, Reg b) { return vandq_s16(a, b); }
DNALN_INLINE Reg r_or(Reg a, Reg b) { return vorrq_s16(a, b); }
DNALN_INLINE Reg r_blend(Reg m, Reg a, Reg b) { return vbslq_s16(vreinterpretq_u16_s16(m), a, b); }
DNALN_INLINE void r_store_lo8(uint8_t* p, Reg a) { vst1_u8(p, vmovn_u16(vreinterpretq_u16_s16(a))); }
#else
inline constexpr const char* kIsa = "portable";
inline constexpr int kNative = 8;
struct Reg {
  int16_t v[8];
};
DNALN_INLINE int16_t sat16(int x) { return static_cast<int16_t>(x > 32767 ? 32767 : (x < -32768 ? -32768 : x)); }
#define DNALN_PORTABLE_OP(NAME, EXPR)                      \
  DNALN_INLINE Reg NAME(Reg a, Reg b) {                     \
    Reg o;                                                  \
    for (int i = 0; i < 8; ++i) o.v[i] = (EXPR);            \
    return o;                                               \
  }
DNALN_PORTABLE_OP(r_adds, sat16(a.v[i] + b.v[i]))
DNALN_PORTABLE_OP(r_subs, sat16(a.v[i] - b.v[i]))
DNALN_PORTABLE_OP(r_max, a.v[i] > b.v[i] ? a.v[i] : b.v[i])
DNALN_PORTABLE_OP(r_eq, static_cast<int16_t>(a.v[i] == b.v[i] ? -1 : 0))
DNALN_PORTABLE_OP(r_gt, static_cast<int16_t>(a.v[i] > b.v[i] ? -1 : 0))
DNALN_PORTABLE_OP(r_and, static_cast<int16_t>(a.v[i] & b.v[i]))
DNALN_PORTABLE_OP(r_or, static_cast<int16_t>(a.v[i] | b.v[i]))
#undef DNALN_PORTABLE_OP
DNALN_INLINE Reg r_set1(int16_t x) { Reg o; for (auto& v : o.v) v = x; return o; }
DNALN_INLINE Reg r_load(const int16_t* p) { Reg o; for (int i = 0; i < 8; ++i) o.v[i] = p[i]; return o; }
DNALN_INLINE void r_store(int16_t* p, Reg a) { for (int i = 0; i < 8; ++i) p[i] = a.v[i]; }
DNALN_INLINE Reg r_blend(Reg m, Reg a, Reg b) { Reg o; for (int i = 0; i < 8; ++i) o.v[i] = m.v[i] ? a.v[i] : b.v[i]; return o; }
DNALN_INLINE void r_store_lo8(uint8_t* p, Reg a) { for (int i = 0; i < 8; ++i) p[i] = static_cast<uint8_t>(a.v[i]); }
#endif

// Packs four lane masks (all-ones / zero) into one byte per lane and stores it:
//   bit7 = m0, bit6 = m1, bit5 = m2, bits4..0 = m3 (0x1F when set).
DNALN_INLINE void r_store_masks8(uint8_t* p, Reg m0, Reg m1, Reg m2, Reg m3) {
#if defined(DNALN_ISA_NEON)
  // Shift-right-and-insert chains the masks into the top bits; one narrowing shift stores them.
  uint16x8_t x = vsriq_n_u16(vreinterpretq_u16_s16(m0), vreinterpretq_u16_s16(m1), 1);
  x = vsriq_n_u16(x, vreinterpretq_u16_s16(m2), 2);
  x = vsriq_n_u16(x, vreinterpretq_u16_s16(m3), 3);
  vst1_u8(p, vshrn_n_u16(x, 8));
#else
  const Reg a = r_or(r_and(m0, r_set1(0x80)), r_and(m1, r_set1(0x40)));
  r_store_lo8(p, r_or(a, r_or(r_and(m2, r_set1(0x20)), r_and(m3, r_set1(0x1F)))));
#endif
}

// ---------------------------------------------------------------------------
// Vec: kUnroll native registers treated as one wide vector of kLanes int16 lanes.
// ---------------------------------------------------------------------------
inline constexpr int kUnroll = DNALN_SIMD_UNROLL;
inline constexpr int kLanes = kNative * kUnroll;

struct Vec {
  Reg r[kUnroll];
};

#define DNALN_VEC_BINOP(NAME, OP)                                  \
  DNALN_INLINE Vec NAME(Vec a, Vec b) {                            \
    Vec o;                                                         \
    for (int k = 0; k < kUnroll; ++k) o.r[k] = OP(a.r[k], b.r[k]); \
    return o;                                                      \
  }
DNALN_VEC_BINOP(vadds, r_adds)
DNALN_VEC_BINOP(vsubs, r_subs)
DNALN_VEC_BINOP(vmax, r_max)
DNALN_VEC_BINOP(veq, r_eq)
DNALN_VEC_BINOP(vgt, r_gt)
DNALN_VEC_BINOP(vand, r_and)
DNALN_VEC_BINOP(vor, r_or)
#undef DNALN_VEC_BINOP

DNALN_INLINE Vec vset1(int16_t x) {
  Vec o;
  for (int k = 0; k < kUnroll; ++k) o.r[k] = r_set1(x);
  return o;
}
DNALN_INLINE Vec vload(const int16_t* p) {
  Vec o;
  for (int k = 0; k < kUnroll; ++k) o.r[k] = r_load(p + k * kNative);
  return o;
}
DNALN_INLINE void vstore(int16_t* p, Vec a) {
  for (int k = 0; k < kUnroll; ++k) r_store(p + k * kNative, a.r[k]);
}
DNALN_INLINE Vec vblend(Vec m, Vec a, Vec b) {
  Vec o;
  for (int k = 0; k < kUnroll; ++k) o.r[k] = r_blend(m.r[k], a.r[k], b.r[k]);
  return o;
}
DNALN_INLINE void vstore_masks8(uint8_t* p, Vec m0, Vec m1, Vec m2, Vec m3) {
  for (int k = 0; k < kUnroll; ++k) r_store_masks8(p + k * kNative, m0.r[k], m1.r[k], m2.r[k], m3.r[k]);
}

// ---------------------------------------------------------------------------
// Substitution scores via a 16-entry byte table lookup (NEON TBX / x86 PSHUFB).
//
// Sequences are stored as bytes, 16 lanes per 16-byte row: query codes pre-shifted
// (code << 2) and target codes as-is, with N encoded as 0x80 in both. Then
// idx = q | t is (a << 2 | b) for two real bases and has bit 7 set if either is N,
// which both TBX (index >= 16 keeps the fallback) and PSHUFB (bit 7 -> 0) handle.
// One OR + one lookup + widening scores 16 lanes, replacing ~5 ops per 8 lanes of
// compare/select logic.
// ---------------------------------------------------------------------------
static_assert(kLanes % 16 == 0, "lane count must be a multiple of 16");
inline constexpr int kByteRows = kLanes / 16;
inline constexpr uint8_t kSeqN = 0x80;

#if defined(DNALN_ISA_NEON)
using ByteRow = uint8x16_t;
struct SubTable {
  int8x16_t tab, amb;
};
DNALN_INLINE SubTable make_subtable(const int8_t* sub16, int8_t ambig_score) {
  return {vld1q_s8(sub16), vdupq_n_s8(ambig_score)};
}
DNALN_INLINE ByteRow bload(const uint8_t* p) { return vld1q_u8(p); }
DNALN_INLINE void score16(const SubTable& st, ByteRow q, const uint8_t* t, Reg& lo, Reg& hi) {
  const int8x16_t s = vqtbx1q_s8(st.amb, st.tab, vorrq_u8(q, vld1q_u8(t)));
  lo = vmovl_s8(vget_low_s8(s));
  hi = vmovl_high_s8(s);
}
#elif defined(DNALN_ISA_AVX2) || defined(DNALN_ISA_SSE41)
using ByteRow = __m128i;
struct SubTable {
  __m128i tab, amb;  // tab holds score + ambig so that PSHUFB's 0 (for N) maps to -ambig
};
DNALN_INLINE SubTable make_subtable(const int8_t* sub16, int8_t ambig_score) {
  alignas(16) int8_t shifted[16];
  for (int i = 0; i < 16; ++i) shifted[i] = static_cast<int8_t>(sub16[i] - ambig_score);
  return {_mm_load_si128(reinterpret_cast<const __m128i*>(shifted)), _mm_set1_epi8(static_cast<char>(-ambig_score))};
}
DNALN_INLINE ByteRow bload(const uint8_t* p) { return _mm_loadu_si128(reinterpret_cast<const __m128i*>(p)); }
DNALN_INLINE __m128i score16_bytes(const SubTable& st, ByteRow q, const uint8_t* t) {
  const __m128i idx = _mm_or_si128(q, _mm_loadu_si128(reinterpret_cast<const __m128i*>(t)));
  return _mm_sub_epi8(_mm_shuffle_epi8(st.tab, idx), st.amb);
}
#else
struct ByteRow {
  uint8_t v[16];
};
struct SubTable {
  int8_t tab[16], amb;
};
DNALN_INLINE SubTable make_subtable(const int8_t* sub16, int8_t ambig_score) {
  SubTable st;
  for (int i = 0; i < 16; ++i) st.tab[i] = sub16[i];
  st.amb = ambig_score;
  return st;
}
DNALN_INLINE ByteRow bload(const uint8_t* p) { ByteRow b; for (int i = 0; i < 16; ++i) b.v[i] = p[i]; return b; }
DNALN_INLINE void score16(const SubTable& st, ByteRow q, const uint8_t* t, Reg& lo, Reg& hi) {
  for (int i = 0; i < 16; ++i) {
    const uint8_t idx = q.v[i] | t[i];
    const int16_t s = idx & 0x80 ? st.amb : st.tab[idx & 15];
    (i < 8 ? lo : hi).v[i & 7] = s;
  }
}
#endif

struct QueryRow {
  ByteRow b[kByteRows];
};
DNALN_INLINE QueryRow qload(const uint8_t* p) {
  QueryRow q;
  for (int k = 0; k < kByteRows; ++k) q.b[k] = bload(p + 16 * k);
  return q;
}

// Substitution scores of one band cell for all kLanes lanes.
DNALN_INLINE Vec vscore(const SubTable& st, const QueryRow& q, const uint8_t* t) {
  Vec o;
#if defined(DNALN_ISA_AVX2)
  for (int k = 0; k < kByteRows; ++k) o.r[k] = _mm256_cvtepi8_epi16(score16_bytes(st, q.b[k], t + 16 * k));
#elif defined(DNALN_ISA_SSE41)
  for (int k = 0; k < kByteRows; ++k) {
    const __m128i s = score16_bytes(st, q.b[k], t + 16 * k);
    o.r[2 * k] = _mm_cvtepi8_epi16(s);
    o.r[2 * k + 1] = _mm_cvtepi8_epi16(_mm_unpackhi_epi64(s, s));
  }
#else
  for (int k = 0; k < kByteRows; ++k) score16(st, q.b[k], t + 16 * k, o.r[2 * k], o.r[2 * k + 1]);
#endif
  return o;
}

}  // namespace dnaln::simd
