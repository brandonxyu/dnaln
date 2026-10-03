#include "common.hpp"

#include <cstdarg>

namespace dnaln {

namespace {
constexpr std::array<char, 256> kCompChar = [] {
  std::array<char, 256> t{};
  for (auto& x : t) x = 'N';
  t['A'] = 'T'; t['C'] = 'G'; t['G'] = 'C'; t['T'] = 'A';
  t['a'] = 't'; t['c'] = 'g'; t['g'] = 'c'; t['t'] = 'a';
  t['n'] = 'n';
  return t;
}();
}  // namespace

void revcomp_append(const std::string& s, std::string& out) {
  const size_t n = s.size(), o = out.size();
  out.resize(o + n);
  for (size_t i = 0; i < n; ++i) out[o + i] = kCompChar[static_cast<uint8_t>(s[n - 1 - i])];
}

std::string cigar_to_string(const std::vector<uint32_t>& cigar) {
  if (cigar.empty()) return "*";
  std::string s;
  for (uint32_t c : cigar) {
    append_int(s, c >> 4);
    s += kCigOps[c & 0xf];
  }
  return s;
}

bool parse_cigar(const char* s, size_t n, std::vector<uint32_t>& out) {
  out.clear();
  if (n == 1 && s[0] == '*') return true;
  uint32_t len = 0;
  bool have_len = false;
  for (size_t i = 0; i < n; ++i) {
    const char c = s[i];
    if (c >= '0' && c <= '9') {
      len = len * 10 + static_cast<uint32_t>(c - '0');
      have_len = true;
      continue;
    }
    const char* p = std::strchr(kCigOps, c);
    if (!p || c == '\0' || !have_len) return false;
    out.push_back(len << 4 | static_cast<uint32_t>(p - kCigOps));
    len = 0;
    have_len = false;
  }
  return !have_len;
}

void die(const char* fmt, ...) {
  std::fflush(stdout);
  std::fputs("[dnaln] error: ", stderr);
  va_list ap;
  va_start(ap, fmt);
  std::vfprintf(stderr, fmt, ap);
  va_end(ap);
  std::fputc('\n', stderr);
  std::exit(1);
}

void log_msg(const char* fmt, ...) {
  std::fputs("[dnaln] ", stderr);
  va_list ap;
  va_start(ap, fmt);
  std::vfprintf(stderr, fmt, ap);
  va_end(ap);
  std::fputc('\n', stderr);
}

}  // namespace dnaln
