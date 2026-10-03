#include "simulate.hpp"

#include <algorithm>
#include <cstdio>
#include <vector>

#include "common.hpp"

namespace dnaln {

namespace {

FILE* open_or_die(const std::string& path, const char* mode) {
  FILE* fp = std::fopen(path.c_str(), mode);
  if (!fp) die("cannot open '%s' for writing", path.c_str());
  return fp;
}

void push_op(std::vector<uint32_t>& c, uint32_t op, uint32_t len) {
  if (!c.empty() && (c.back() & 0xf) == op) c.back() += len << 4;
  else c.push_back(len << 4 | op);
}

// Builds one read (reference orientation) starting at global position `start`.
// Indels are never placed at the read ends, so the true CIGAR starts and ends with M.
bool make_read(const Reference& ref, const Contig& c, int64_t start, const SimParams& p, Rng& rng,
               std::string& read, std::vector<uint32_t>& cigar, int64_t& end) {
  read.clear();
  cigar.clear();
  const int64_t limit = static_cast<int64_t>(c.offset + c.len);
  int64_t rp = start;
  uint32_t last = kCigM;
  while (static_cast<int>(read.size()) < p.len) {
    if (rp >= limit) return false;
    if (!read.empty() && last == kCigM && rng.uniform() < p.indel) {
      int L = 1;
      while (L < 8 && rng.uniform() < p.indel_ext) ++L;
      if (rng.uniform() < 0.5) {
        if (static_cast<int>(read.size()) + L >= p.len) continue;
        for (int i = 0; i < L; ++i) read += kBase[rng.below(4)];
        push_op(cigar, kCigI, L);
        last = kCigI;
      } else {
        rp += L;
        push_op(cigar, kCigD, L);
        last = kCigD;
      }
      continue;
    }
    uint8_t b = ref.seq[static_cast<size_t>(rp)];
    if (b > 3) return false;  // avoid N regions
    if (rng.uniform() < p.sub) b = static_cast<uint8_t>((b + 1 + rng.below(3)) & 3);
    read += kBase[b];
    push_op(cigar, kCigM, 1);
    last = kCigM;
    ++rp;
  }
  end = rp;
  return true;
}

}  // namespace

void simulate_reads(const Reference& ref, const SimParams& p, const std::string& prefix) {
  if (p.len < 20) die("read length must be >= 20");
  // Sample read start positions uniformly over all contigs long enough to hold a read.
  std::vector<int> contig;
  std::vector<uint64_t> cum;
  uint64_t total = 0;
  for (size_t i = 0; i < ref.contigs.size(); ++i) {
    const uint64_t L = ref.contigs[i].len;
    if (L < static_cast<uint64_t>(p.len) + 100) continue;
    total += L - p.len - 100;
    contig.push_back(static_cast<int>(i));
    cum.push_back(total);
  }
  if (total == 0) die("reference has no contig long enough for %d bp reads", p.len);

  FILE* fq = open_or_die(prefix + ".fq", "w");
  FILE* tr = open_or_die(prefix + ".truth.tsv", "w");
  std::fprintf(tr, "#name\tcontig\tstart\tend\tstrand\tcigar\n");
  Rng rng(p.seed);
  std::string read, seq, qual(static_cast<size_t>(p.len), 'I');
  std::vector<uint32_t> cigar;
  for (uint64_t r = 0; r < p.n; ++r) {
    int64_t start = 0, end = 0;
    const Contig* c = nullptr;
    for (int attempt = 0;; ++attempt) {
      if (attempt > 10000) die("could not sample a read without Ns; is the reference mostly N?");
      const uint64_t x = rng.below(total);
      const size_t ci = static_cast<size_t>(std::upper_bound(cum.begin(), cum.end(), x) - cum.begin());
      c = &ref.contigs[contig[ci]];
      start = static_cast<int64_t>(c->offset + (x - (ci ? cum[ci - 1] : 0)));
      if (make_read(ref, *c, start, p, rng, read, cigar, end)) break;
    }
    const bool rev = rng.next() & 1;
    seq.clear();
    if (rev) revcomp_append(read, seq);
    else seq = read;
    std::fprintf(fq, "@sim_%llu\n%s\n+\n%s\n", static_cast<unsigned long long>(r + 1), seq.c_str(), qual.c_str());
    std::fprintf(tr, "sim_%llu\t%s\t%lld\t%lld\t%c\t%s\n", static_cast<unsigned long long>(r + 1), c->name.c_str(),
                 static_cast<long long>(start - static_cast<int64_t>(c->offset)),
                 static_cast<long long>(end - static_cast<int64_t>(c->offset)), rev ? '-' : '+',
                 cigar_to_string(cigar).c_str());
  }
  std::fclose(fq);
  std::fclose(tr);
}

void random_reference(uint64_t len, double gc, uint64_t seed, const std::string& path) {
  Rng rng(seed);
  auto base = [&] {
    const double u = rng.uniform();
    return u < gc / 2 ? 'G' : u < gc ? 'C' : u < gc + (1 - gc) / 2 ? 'A' : 'T';
  };
  std::string s(len, 'A');
  for (auto& ch : s) ch = base();

  // Repeat families loosely modelled on E. coli K-12: 7 rRNA operons (~5 kb, near
  // identical), ~10 IS elements (~1.3 kb), and many short REP-like elements.
  struct Family {
    int copies, unit_len;
    double divergence;
  };
  const Family families[] = {{7, 5000, 0.002}, {10, 1300, 0.001}, {40, 150, 0.03}};
  for (const Family& f : families) {
    if (static_cast<uint64_t>(f.unit_len) * 4 > len) continue;
    std::string unit(static_cast<size_t>(f.unit_len), 'A');
    for (auto& ch : unit) ch = base();
    for (int c = 0; c < f.copies; ++c) {
      std::string copy = unit;
      for (auto& ch : copy)
        if (rng.uniform() < f.divergence) ch = "ACGT"[(std::string("ACGT").find(ch) + 1 + rng.below(3)) & 3];
      if (rng.next() & 1) {
        std::string rc;
        revcomp_append(copy, rc);
        copy.swap(rc);
      }
      const uint64_t pos = rng.below(len - copy.size());
      s.replace(pos, copy.size(), copy);
    }
  }

  FILE* fp = open_or_die(path, "w");
  std::fprintf(fp, ">synthetic_ecoli_like len=%llu gc=%.3f seed=%llu\n", static_cast<unsigned long long>(len), gc,
               static_cast<unsigned long long>(seed));
  for (uint64_t i = 0; i < len; i += 80) {
    std::fwrite(s.data() + i, 1, std::min<uint64_t>(80, len - i), fp);
    std::fputc('\n', fp);
  }
  std::fclose(fp);
}

}  // namespace dnaln
