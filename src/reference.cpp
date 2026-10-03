#include "reference.hpp"

#include "common.hpp"
#include "seqio.hpp"

namespace dnaln {

Reference Reference::load(const std::string& path) {
  Reference ref;
  SeqReader rd(path);
  SeqRecord rec;
  while (rd.next(rec)) {
    if (rec.seq.size() > UINT32_MAX) die("contig '%s' is too long", rec.name.c_str());
    const size_t off = ref.seq.size();
    ref.contigs.push_back({rec.name, off, static_cast<uint32_t>(rec.seq.size())});
    ref.seq.resize(off + rec.seq.size());
    encode(rec.seq.data(), rec.seq.size(), ref.seq.data() + off);
  }
  if (ref.contigs.empty()) die("no sequences found in '%s'", path.c_str());
  // Index positions are stored as (pos << 1 | strand) in 32 bits.
  if (ref.seq.size() >= (1ULL << 31)) die("reference is too large (limit: 2^31 bp)");
  return ref;
}

int Reference::contig_of(uint64_t gpos) const {
  size_t lo = 0, hi = contigs.size();
  while (hi - lo > 1) {
    const size_t mid = (lo + hi) / 2;
    if (contigs[mid].offset <= gpos) lo = mid; else hi = mid;
  }
  return static_cast<int>(lo);
}

bool Reference::inside(int ci, int64_t gstart, int len) const {
  const Contig& c = contigs[ci];
  return gstart >= static_cast<int64_t>(c.offset) &&
         gstart + len <= static_cast<int64_t>(c.offset + c.len);
}

void Reference::fetch(int ci, int64_t gstart, int len, uint8_t* out) const {
  const Contig& c = contigs[ci];
  const int64_t b = static_cast<int64_t>(c.offset), e = b + c.len;
  for (int i = 0; i < len; ++i) {
    const int64_t p = gstart + i;
    out[i] = (p >= b && p < e) ? seq[static_cast<size_t>(p)] : 4;
  }
}

int Reference::find(const std::string& name) const {
  for (size_t i = 0; i < contigs.size(); ++i)
    if (contigs[i].name == name) return static_cast<int>(i);
  return -1;
}

}  // namespace dnaln
