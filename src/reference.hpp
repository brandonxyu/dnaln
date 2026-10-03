// reference.hpp — reference genome held as one concatenated array of 0..4 base codes.
#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace dnaln {

struct Contig {
  std::string name;
  uint64_t offset = 0;  // start in the concatenated sequence
  uint32_t len = 0;
};

class Reference {
 public:
  std::vector<Contig> contigs;
  std::vector<uint8_t> seq;  // concatenated codes (A=0 C=1 G=2 T=3 N=4)

  static Reference load(const std::string& path);

  // Index of the contig containing global position `gpos`.
  int contig_of(uint64_t gpos) const;
  // True if [gstart, gstart+len) lies entirely inside contig `ci`.
  bool inside(int ci, int64_t gstart, int len) const;
  // Copies [gstart, gstart+len) into `out`; positions outside contig `ci` become N.
  void fetch(int ci, int64_t gstart, int len, uint8_t* out) const;
  int find(const std::string& name) const;  // -1 if absent
};

}  // namespace dnaln
