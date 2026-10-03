#include "evaluate.hpp"

#include <algorithm>
#include <cstdlib>
#include <cstring>

#include "common.hpp"
#include "seqio.hpp"

namespace dnaln {

TruthMap load_truth(const std::string& path) {
  TruthMap m;
  GzLineReader lr(path);
  std::string line;
  while (lr.read_line(line)) {
    if (line.empty() || line[0] == '#') continue;
    std::vector<std::string> f;
    size_t s = 0;
    for (;;) {
      const size_t t = line.find('\t', s);
      f.push_back(line.substr(s, t == std::string::npos ? std::string::npos : t - s));
      if (t == std::string::npos) break;
      s = t + 1;
    }
    if (f.size() < 6) die("malformed truth line: %s", line.c_str());
    TruthRec r;
    r.contig = f[1];
    r.start = std::atoll(f[2].c_str());
    r.end = std::atoll(f[3].c_str());
    r.rev = f[4] == "-";
    if (!parse_cigar(f[5].data(), f[5].size(), r.cigar)) die("bad CIGAR in truth: %s", f[5].c_str());
    m.emplace(f[0], std::move(r));
  }
  return m;
}

Evaluator::Evaluator(const TruthMap& truth, int tol) : truth_(truth), tol_(tol) {}

namespace {
// Reference position of every read base (reference orientation); -1 = inserted/clipped.
void base_positions(int64_t pos, const std::vector<uint32_t>& cigar, std::vector<int64_t>& out) {
  out.clear();
  for (uint32_t c : cigar) {
    const uint32_t op = c & 0xf, len = c >> 4;
    if (op == kCigM || op == kCigEq || op == kCigX) {
      for (uint32_t k = 0; k < len; ++k) out.push_back(pos++);
    } else if (op == kCigI || op == kCigS) {
      out.insert(out.end(), len, -1);
    } else if (op == kCigD || op == kCigN) {
      pos += len;
    }
  }
}
}  // namespace

void Evaluator::add(const std::string& name, bool mapped, const std::string& contig, int64_t pos0, bool rev, int mapq,
                    const std::vector<uint32_t>& cigar) {
  const auto it = truth_.find(name);
  if (it == truth_.end()) return;
  ++records_;
  if (!mapped || cigar.empty()) return;
  const TruthRec& t = it->second;
  ++mapped_;
  int64_t clip = 0;
  const uint32_t op0 = cigar.front() & 0xf;
  if (op0 == kCigS || op0 == kCigH) clip = cigar.front() >> 4;
  const bool ok = contig == t.contig && rev == t.rev && std::llabs(pos0 - clip - t.start) <= tol_;
  const int q = std::clamp(mapq, 0, 60);
  ++mq_mapped_[q];
  if (!ok) {
    ++mq_wrong_[q];
    return;
  }
  ++correct_;
  base_positions(t.start, t.cigar, tp_);
  base_positions(pos0, cigar, pp_);
  const size_t n = std::min(tp_.size(), pp_.size());
  for (size_t i = 0; i < n; ++i) {
    if (tp_[i] < 0) continue;
    ++bases_;
    bases_ok_ += pp_[i] == tp_[i];
  }
}

void Evaluator::add_sam_file(const std::string& path) {
  GzLineReader lr(path);
  std::string line, name, contig;
  std::vector<uint32_t> cigar;
  while (lr.read_line(line)) {
    if (line.empty() || line[0] == '@') continue;
    const char* f[6];
    size_t flen[6];
    const char* p = line.c_str();
    int nf = 0;
    while (nf < 6) {
      const char* tab = std::strchr(p, '\t');
      f[nf] = p;
      flen[nf] = tab ? static_cast<size_t>(tab - p) : std::strlen(p);
      ++nf;
      if (!tab) break;
      p = tab + 1;
    }
    if (nf < 6) die("malformed SAM line in %s", path.c_str());
    const int flag = std::atoi(f[1]);
    if (flag & 0x900) continue;  // secondary / supplementary
    name.assign(f[0], flen[0]);
    contig.assign(f[2], flen[2]);
    if (!parse_cigar(f[5], flen[5], cigar)) die("bad CIGAR in %s", path.c_str());
    add(name, !(flag & 4), contig, std::atoll(f[3]) - 1, flag & 16, std::atoi(f[4]), cigar);
  }
}

uint64_t Evaluator::mapped_at(int q) const {
  uint64_t s = 0;
  for (int i = std::max(q, 0); i <= 60; ++i) s += mq_mapped_[i];
  return s;
}
uint64_t Evaluator::wrong_at(int q) const {
  uint64_t s = 0;
  for (int i = std::max(q, 0); i <= 60; ++i) s += mq_wrong_[i];
  return s;
}

void print_eval_report(FILE* fp, const std::vector<std::string>& labels, const std::vector<const Evaluator*>& evs,
                       double max_gap_pp) {
  auto pct = [](uint64_t a, uint64_t b) { return b ? 100.0 * static_cast<double>(a) / b : 0.0; };
  std::fprintf(fp, "%-34s", "metric");
  for (const auto& l : labels) {
    const std::string s = l.size() > 22 ? "..." + l.substr(l.size() - 19) : l;
    std::fprintf(fp, " %22s", s.c_str());
  }
  std::fprintf(fp, "\n");
  auto row = [&](const char* label, auto fn) {
    std::fprintf(fp, "%-34s", label);
    for (const Evaluator* e : evs) std::fprintf(fp, " %22s", fn(*e).c_str());
    std::fprintf(fp, "\n");
  };
  char buf[64];
  row("reads in truth set", [&](const Evaluator& e) { return std::to_string(e.n_truth()); });
  row("primary records found", [&](const Evaluator& e) { return std::to_string(e.records()); });
  row("mapped", [&](const Evaluator& e) {
    std::snprintf(buf, sizeof buf, "%.3f%%", pct(e.mapped(), e.n_truth()));
    return std::string(buf);
  });
  row("ACCURACY (correct / all reads)", [&](const Evaluator& e) {
    std::snprintf(buf, sizeof buf, "%.3f%%", 100.0 * e.accuracy());
    return std::string(buf);
  });
  row("wrong (% of mapped)", [&](const Evaluator& e) {
    std::snprintf(buf, sizeof buf, "%.3f%%", pct(e.mapped() - e.correct(), e.mapped()));
    return std::string(buf);
  });
  for (int q : {60, 30, 10, 1}) {
    char label[64];
    std::snprintf(label, sizeof label, "MAPQ>=%-2d mapped | error rate", q);
    row(label, [&](const Evaluator& e) {
      std::snprintf(buf, sizeof buf, "%.2f%% | %.4f%%", pct(e.mapped_at(q), e.n_truth()),
                    pct(e.wrong_at(q), e.mapped_at(q)));
      return std::string(buf);
    });
  }
  row("base-level concordance", [&](const Evaluator& e) {
    std::snprintf(buf, sizeof buf, "%.3f%%", 100.0 * e.base_concordance());
    return std::string(buf);
  });
  if (evs.size() >= 2) {
    const Evaluator& ref = *evs.back();
    std::fprintf(fp, "\nAccuracy relative to %s:\n", labels.back().c_str());
    for (size_t i = 0; i + 1 < evs.size(); ++i) {
      const double gap = 100.0 * (evs[i]->accuracy() - ref.accuracy());
      std::fprintf(fp, "  %-28s %+.3f percentage points -> %s (target: within %.1f pp)\n", labels[i].c_str(), gap,
                   gap >= -max_gap_pp ? "PASS" : "FAIL", max_gap_pp);
    }
  }
}

}  // namespace dnaln
