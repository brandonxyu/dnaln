// evaluate.hpp — mapping accuracy against simulated ground truth (works on any SAM,
// so dnaln and minimap2 are scored by exactly the same code).
//
// A primary alignment is "correct" when it is on the true contig and strand and its
// unclipped start (POS minus leading soft/hard clip) is within `tol` bp of the true
// start. For correct reads we also measure base-level concordance: the fraction of
// read bases placed on exactly their true reference position.
#pragma once

#include <cstdint>
#include <cstdio>
#include <string>
#include <unordered_map>
#include <vector>

namespace dnaln {

struct TruthRec {
  std::string contig;
  int64_t start = 0, end = 0;
  bool rev = false;
  std::vector<uint32_t> cigar;  // reference orientation
};
using TruthMap = std::unordered_map<std::string, TruthRec>;

TruthMap load_truth(const std::string& path);

class Evaluator {
 public:
  explicit Evaluator(const TruthMap& truth, int tol = 10);
  void add(const std::string& name, bool mapped, const std::string& contig, int64_t pos0, bool rev, int mapq,
           const std::vector<uint32_t>& cigar);
  void add_sam_file(const std::string& path);

  uint64_t n_truth() const { return truth_.size(); }
  uint64_t records() const { return records_; }
  uint64_t mapped() const { return mapped_; }
  uint64_t correct() const { return correct_; }
  double accuracy() const { return truth_.empty() ? 0.0 : static_cast<double>(correct_) / truth_.size(); }
  double base_concordance() const { return bases_ ? static_cast<double>(bases_ok_) / bases_ : 0.0; }
  // Cumulative counts over alignments with MAPQ >= q.
  uint64_t mapped_at(int q) const;
  uint64_t wrong_at(int q) const;

 private:
  const TruthMap& truth_;
  int tol_;
  uint64_t records_ = 0, mapped_ = 0, correct_ = 0, bases_ = 0, bases_ok_ = 0;
  uint64_t mq_mapped_[61] = {}, mq_wrong_[61] = {};
  std::vector<int64_t> tp_, pp_;
};

// Side-by-side report; with >= 2 inputs also prints the accuracy delta vs the last one.
void print_eval_report(FILE* fp, const std::vector<std::string>& labels, const std::vector<const Evaluator*>& evs,
                       double max_gap_pp = 2.0);

}  // namespace dnaln
