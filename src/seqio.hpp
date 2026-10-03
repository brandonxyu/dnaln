// seqio.hpp — buffered FASTA/FASTQ reading (plain or gzip, via zlib).
#pragma once

#include <zlib.h>

#include <cstdint>
#include <string>
#include <vector>

namespace dnaln {

struct SeqRecord {
  std::string name, seq, qual;  // qual is empty for FASTA input
};

// Line reader over a (possibly gzipped) file; "-" reads stdin.
class GzLineReader {
 public:
  explicit GzLineReader(const std::string& path);
  ~GzLineReader();
  GzLineReader(const GzLineReader&) = delete;
  GzLineReader& operator=(const GzLineReader&) = delete;

  bool read_line(std::string& line);  // false at EOF; strips trailing \r\n
  int peek();                         // next byte, or -1 at EOF

 private:
  bool fill();
  gzFile fp_ = nullptr;
  std::vector<char> buf_;
  size_t pos_ = 0, len_ = 0;
  bool eof_ = false;
};

// Streams FASTA and FASTQ records (formats may be mixed; multi-line FASTA supported).
class SeqReader : public GzLineReader {
 public:
  using GzLineReader::GzLineReader;
  bool next(SeqRecord& rec);
  // Fills up to `max` records into `out` (reusing its strings); returns how many were read.
  size_t next_batch(std::vector<SeqRecord>& out, size_t max);

 private:
  std::string line_;
};

std::vector<SeqRecord> read_all_records(const std::string& path, size_t limit = SIZE_MAX);

}  // namespace dnaln
