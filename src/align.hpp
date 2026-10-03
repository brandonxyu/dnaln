// align.hpp — banded Smith–Waterman (local alignment, affine gaps).
//
// Band geometry. For query row i the band covers target columns j = i + d with
// d in [0, W), W = 2*bw + 1. The target window is placed so the expected diagonal
// sits at d = bw (window start = seed diagonal - bw), so the window has length
// qlen + 2*bw. In these "diagonal coordinates" the three DP predecessors are
//     diagonal (i-1, j-1) -> (i-1, d)
//     up       (i-1, j)   -> (i-1, d+1)     E: gap in the target (CIGAR I)
//     left     (i, j-1)   -> (i,   d-1)     F: gap in the query  (CIGAR D)
// so the diagonal and E terms come from the previous row and the rows can be
// updated in place; only F carries a dependency along the row.
//
//   H(i,d) = max(0, H(i-1,d) + s(q_i, t_{i+d}), E(i,d), F(i,d))
//   E(i,d) = max(H(i-1,d+1) - (go+ge), E(i-1,d+1) - ge)
//   F(i,d) = max(H(i,d-1)   - (go+ge), F(i,d-1)   - ge)
//
// Two implementations with bit-identical results:
//   ScalarAligner — one alignment at a time (the baseline).
//   SimdAligner   — inter-sequence vectorization: kLanes independent alignments, one
//                   per 16-bit SIMD lane, computed in lock-step (the approach of
//                   BWA-MEM2's banded SW). No intra-row shuffles are needed and every
//                   lane does exactly the scalar recurrence.
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "common.hpp"

namespace dnaln {

// A gap of length L costs gap_open + L * gap_ext.
struct ScoreParams {
  int match = 2, mismatch = 4, gap_open = 4, gap_ext = 2, ambig = 1;
  int bw = 16;  // band half-width

  int band_cols() const { return 2 * bw + 1; }
  int sub(uint8_t a, uint8_t b) const { return (a > 3 || b > 3) ? -ambig : (a == b ? match : -mismatch); }
};

struct AlignTask {
  const uint8_t* q = nullptr;  // query codes
  int32_t qlen = 0;
  const uint8_t* t = nullptr;  // target window codes; should hold qlen + 2*bw (shorter is N-padded)
  int32_t tlen = 0;
};

struct AlignResult {
  int32_t score = 0;
  int32_t qb = 0, qe = 0;  // aligned query range [qb, qe)
  int32_t tb = 0, te = 0;  // aligned target range [tb, te), window coordinates
  int32_t nm = 0;          // edit distance of the aligned part
  std::vector<uint32_t> cigar;  // M/I/D only

  void clear() {
    score = qb = qe = tb = te = nm = 0;
    cigar.clear();
  }
};

class Aligner {
 public:
  virtual ~Aligner() = default;
  // Aligns n independent tasks. With traceback=false only `score` is computed.
  virtual void align(const AlignTask* tasks, size_t n, bool traceback, AlignResult* out) = 0;
  virtual std::string name() const = 0;
};

class ScalarAligner final : public Aligner {
 public:
  explicit ScalarAligner(const ScoreParams& sp);
  void align(const AlignTask* tasks, size_t n, bool traceback, AlignResult* out) override;
  std::string name() const override { return "scalar"; }

 private:
  template <bool TB>
  void align_one(const AlignTask& tk, AlignResult& r);
  ScoreParams sp_;
  int mat_[5 * 8];
  AlignedBuffer<int32_t> h_, e_;
  AlignedBuffer<uint8_t> t_, dir_;
};

class SimdAligner final : public Aligner {
 public:
  explicit SimdAligner(const ScoreParams& sp);
  void align(const AlignTask* tasks, size_t n, bool traceback, AlignResult* out) override;
  std::string name() const override;
  static int lanes();
  static const char* isa();
  // Bytes of kernel scratch touched per lane-group for a given query length (cache analysis).
  size_t working_set_bytes(int qlen, bool traceback) const;

 private:
  template <bool TB>
  void align_group(const AlignTask* tasks, const uint32_t* idx, int n, AlignResult* out);
  ScoreParams sp_;
  ScalarAligner fallback_;  // reads too long for 16-bit scores
  AlignedBuffer<int16_t> h_, e_, dconst_;
  AlignedBuffer<uint8_t> q_, t_, dir_;
  int8_t sub16_[16];
  bool byte_scores_ok_;
};

// Direction byte written for every band cell when tracing back. No "H == 0" bit is
// needed: the traceback tracks the running score and stops when it reaches 0.
enum : uint8_t {
  kDirDiag = 0x80,  // H came from the diagonal (match/mismatch)
  kDirE = 0x40,     // otherwise: H came from E (insertion) if set, else from F (deletion)
  kDirExtE = 0x20,  // E at this cell extends E of (i-1, d+1)
  kDirExtF = 0x10,  // F at this cell extends F of (i, d-1)   (written as 0x1F)
};

// Walks a direction matrix stored as dir[(i * W + d) * stride + lane] back from the
// best cell (i_end, d_end), whose score is `score`.
void traceback(const uint8_t* dir, int stride, int lane, int W, int i_end, int d_end, int score,
               const uint8_t* q, const uint8_t* t, const ScoreParams& sp, AlignResult& r);

// Recomputes an alignment's score from its CIGAR (validation helper).
int score_cigar(const uint8_t* q, const uint8_t* t, int qb, int tb, const std::vector<uint32_t>& cigar,
                const ScoreParams& sp);

}  // namespace dnaln
