// common.hpp — shared utilities: nucleotide encoding, CIGAR helpers, timing, RNG,
// cache-line-aligned scratch buffers.
#pragma once

#include <algorithm>
#include <array>
#include <charconv>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

namespace dnaln {

inline constexpr const char* kVersion = "0.1.0";

// ---------------------------------------------------------------------------
// Nucleotides: A=0 C=1 G=2 T=3; anything else (N, IUPAC codes) = 4.
// ---------------------------------------------------------------------------
inline constexpr std::array<uint8_t, 256> kNt4 = [] {
  std::array<uint8_t, 256> t{};
  for (auto& x : t) x = 4;
  t['A'] = t['a'] = 0;
  t['C'] = t['c'] = 1;
  t['G'] = t['g'] = 2;
  t['T'] = t['t'] = t['U'] = t['u'] = 3;
  return t;
}();
inline constexpr char kBase[] = "ACGTN";

inline uint8_t comp_code(uint8_t c) { return c < 4 ? 3 - c : 4; }

inline void encode(const char* s, size_t n, uint8_t* out) {
  for (size_t i = 0; i < n; ++i) out[i] = kNt4[static_cast<uint8_t>(s[i])];
}

inline void revcomp_codes(const uint8_t* in, size_t n, uint8_t* out) {
  for (size_t i = 0; i < n; ++i) out[i] = comp_code(in[n - 1 - i]);
}

// Appends the reverse complement of an ASCII sequence to `out`.
void revcomp_append(const std::string& s, std::string& out);

// ---------------------------------------------------------------------------
// CIGAR, BAM-style encoding: (length << 4) | op.
// ---------------------------------------------------------------------------
enum : uint32_t { kCigM = 0, kCigI = 1, kCigD = 2, kCigN = 3, kCigS = 4, kCigH = 5, kCigEq = 7, kCigX = 8 };
inline constexpr char kCigOps[] = "MIDNSHP=X";

std::string cigar_to_string(const std::vector<uint32_t>& cigar);
bool parse_cigar(const char* s, size_t n, std::vector<uint32_t>& out);

inline void append_int(std::string& s, long long v) {
  char buf[24];
  const auto r = std::to_chars(buf, buf + sizeof buf, v);
  s.append(buf, r.ptr);
}

// ---------------------------------------------------------------------------
// Logging / errors
// ---------------------------------------------------------------------------
[[noreturn]] void die(const char* fmt, ...) __attribute__((format(printf, 1, 2)));
void log_msg(const char* fmt, ...) __attribute__((format(printf, 1, 2)));

class Timer {
 public:
  Timer() : t0_(clock::now()) {}
  void reset() { t0_ = clock::now(); }
  double sec() const { return std::chrono::duration<double>(clock::now() - t0_).count(); }

 private:
  using clock = std::chrono::steady_clock;
  clock::time_point t0_;
};

// xoshiro256** — small, fast, and reproducible across platforms (unlike std::*_distribution).
class Rng {
 public:
  explicit Rng(uint64_t seed) {
    for (auto& x : s_) {  // splitmix64 seeding
      seed += 0x9e3779b97f4a7c15ULL;
      uint64_t z = seed;
      z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
      z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
      x = z ^ (z >> 31);
    }
  }
  uint64_t next() {
    const uint64_t r = rotl(s_[1] * 5, 7) * 9, t = s_[1] << 17;
    s_[2] ^= s_[0];
    s_[3] ^= s_[1];
    s_[1] ^= s_[2];
    s_[0] ^= s_[3];
    s_[2] ^= t;
    s_[3] = rotl(s_[3], 45);
    return r;
  }
  double uniform() { return static_cast<double>(next() >> 11) * 0x1.0p-53; }
  uint64_t below(uint64_t n) { return next() % n; }

 private:
  static uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
  uint64_t s_[4];
};

// Grow-only, 64-byte (cache-line) aligned scratch buffer for the alignment kernels.
template <typename T>
class AlignedBuffer {
 public:
  AlignedBuffer() = default;
  AlignedBuffer(const AlignedBuffer&) = delete;
  AlignedBuffer& operator=(const AlignedBuffer&) = delete;
  ~AlignedBuffer() { std::free(p_); }

  T* get(size_t n) {
    if (n > cap_) {
      std::free(p_);
      const size_t bytes = ((std::max<size_t>(n, 1) * sizeof(T) + 63) / 64) * 64;
      p_ = static_cast<T*>(std::aligned_alloc(64, bytes));
      if (!p_) die("out of memory allocating %zu bytes", bytes);
      cap_ = n;
    }
    return p_;
  }
  size_t capacity_bytes() const { return cap_ * sizeof(T); }

 private:
  T* p_ = nullptr;
  size_t cap_ = 0;
};

// FNV-1a; used for deterministic tie-breaking between equally good hits.
inline uint32_t fnv1a(const std::string& s) {
  uint32_t h = 2166136261u;
  for (unsigned char c : s) h = (h ^ c) * 16777619u;
  return h;
}

}  // namespace dnaln
