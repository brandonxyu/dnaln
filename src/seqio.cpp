#include "seqio.hpp"

#include <cstdio>
#include <cstring>

#include "common.hpp"

namespace dnaln {

GzLineReader::GzLineReader(const std::string& path) : buf_(1 << 20) {
  fp_ = path == "-" ? gzdopen(fileno(stdin), "rb") : gzopen(path.c_str(), "rb");
  if (!fp_) die("cannot open '%s'", path.c_str());
  gzbuffer(fp_, 1 << 18);
}

GzLineReader::~GzLineReader() {
  if (fp_) gzclose(fp_);
}

bool GzLineReader::fill() {
  if (eof_) return false;
  const int n = gzread(fp_, buf_.data(), static_cast<unsigned>(buf_.size()));
  if (n < 0) die("read error: %s", gzerror(fp_, nullptr));
  if (n == 0) {
    eof_ = true;
    return false;
  }
  pos_ = 0;
  len_ = static_cast<size_t>(n);
  return true;
}

int GzLineReader::peek() {
  if (pos_ >= len_ && !fill()) return -1;
  return static_cast<unsigned char>(buf_[pos_]);
}

bool GzLineReader::read_line(std::string& out) {
  out.clear();
  bool got = false;
  for (;;) {
    if (pos_ >= len_ && !fill()) break;
    got = true;
    const char* start = buf_.data() + pos_;
    const size_t avail = len_ - pos_;
    const char* nl = static_cast<const char*>(std::memchr(start, '\n', avail));
    if (nl) {
      out.append(start, static_cast<size_t>(nl - start));
      pos_ += static_cast<size_t>(nl - start) + 1;
      break;
    }
    out.append(start, avail);
    pos_ = len_;
  }
  if (!out.empty() && out.back() == '\r') out.pop_back();
  return got;
}

bool SeqReader::next(SeqRecord& r) {
  int c;
  while ((c = peek()) != -1 && c != '>' && c != '@') read_line(line_);  // skip blank/junk lines
  if (c == -1) return false;
  read_line(line_);
  const bool fastq = line_[0] == '@';
  const size_t ws = line_.find_first_of(" \t", 1);
  r.name.assign(line_, 1, ws == std::string::npos ? std::string::npos : ws - 1);
  r.seq.clear();
  r.qual.clear();
  while ((c = peek()) != -1 && c != '>' && c != '+' && c != '@') {
    read_line(line_);
    r.seq += line_;
  }
  if (fastq) {
    if (c != '+') die("malformed FASTQ record '%s'", r.name.c_str());
    read_line(line_);  // '+' separator
    while (r.qual.size() < r.seq.size() && read_line(line_)) r.qual += line_;
    if (r.qual.size() != r.seq.size()) die("FASTQ record '%s': quality length mismatch", r.name.c_str());
  }
  return true;
}

size_t SeqReader::next_batch(std::vector<SeqRecord>& out, size_t max) {
  if (out.size() < max) out.resize(max);
  size_t n = 0;
  while (n < max && next(out[n])) ++n;
  return n;
}

std::vector<SeqRecord> read_all_records(const std::string& path, size_t limit) {
  SeqReader rd(path);
  std::vector<SeqRecord> v;
  SeqRecord r;
  while (v.size() < limit && rd.next(r)) v.push_back(std::move(r));
  return v;
}

}  // namespace dnaln
