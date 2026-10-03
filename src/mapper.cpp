#include "mapper.hpp"

#include <algorithm>
#include <cstdlib>

namespace dnaln {

Mapper::Mapper(const Reference& ref, const MinimizerIndex& idx, const MapParams& p) : ref_(ref), idx_(idx), p_(p) {
  p_.ch.band = std::min(p_.ch.band, p_.sc.bw);  // never chain across more drift than the SW band covers
  if (p.simd) aligner_ = std::make_unique<SimdAligner>(p_.sc);
  else aligner_ = std::make_unique<ScalarAligner>(p_.sc);
}

void Mapper::seed_and_chain(const uint8_t* fwd, int len, uint32_t rid) {
  const int k = idx_.k();
  mins_.clear();
  sketch(fwd, static_cast<uint32_t>(len), k, idx_.w(), 0, mins_);

  // Two-stage software prefetch: all hash buckets first, then all occurrence lists,
  // so the independent cache misses of one read overlap instead of serializing.
  if (p_.prefetch)
    for (const Minimizer& m : mins_) idx_.prefetch(m.hash);
  hits_.clear();
  for (uint32_t i = 0; i < mins_.size(); ++i) {
    uint32_t n;
    const uint32_t* pos = idx_.lookup(mins_[i].hash, n);
    if (n == 0 || n > p_.max_occ) continue;
    if (p_.prefetch) __builtin_prefetch(pos);
    hits_.push_back({pos, n, i});
  }
  anchors_.clear();
  for (const Hit& h : hits_) {
    const Minimizer& m = mins_[h.mi];
    for (uint32_t j = 0; j < h.n; ++j) {
      const uint32_t v = h.pos[j];
      const uint32_t rev = (v & 1) ^ m.strand;
      const int32_t qpos = rev ? len - static_cast<int32_t>(m.pos) - k : static_cast<int32_t>(m.pos);
      anchors_.push_back({v >> 1, qpos, rev});
    }
  }
  stats_.anchors += anchors_.size();

  chainer_.run(anchors_, k, p_.ch, chains_);
  if (chains_.empty()) return;
  const int bw = p_.sc.bw;
  const double min_score = chains_[0].score * p_.pri_ratio;
  const size_t first = cands_.size();
  for (const Chain& c : chains_) {
    if (c.score < min_score || static_cast<int>(cands_.size() - first) >= p_.max_cand) break;
    const int64_t diag = c.dmin + (c.dmax - c.dmin) / 2;
    bool dup = false;
    for (size_t i = first; i < cands_.size() && !dup; ++i)
      dup = cands_[i].rev == c.strand && std::llabs(cands_[i].diag - diag) <= bw;
    if (dup) continue;
    Cand cd;
    cd.read = rid;
    cd.rev = c.strand;
    cd.contig = ref_.contig_of(c.rpos);
    cd.diag = diag;
    cd.wstart = diag - bw;
    cands_.push_back(cd);
  }
}

AlignTask Mapper::task_for(const Cand& c) const {
  AlignTask t;
  const int len = qlen_[c.read];
  t.q = codes_.data() + qoff_[c.read] + (c.rev ? len : 0);
  t.qlen = len;
  t.t = c.t;
  t.tlen = len + 2 * p_.sc.bw;
  return t;
}

void Mapper::map(const SeqRecord* reads, size_t n, std::vector<MapRecord>& out) {
  Timer t_seed;
  out.resize(n);
  qoff_.resize(n);
  qlen_.resize(n);
  size_t total = 0;
  for (size_t r = 0; r < n; ++r) total += reads[r].seq.size();
  codes_.resize(2 * total);
  cands_.clear();
  cbeg_.resize(n + 1);

  // 1-2. Encode both strands, seed, chain.
  size_t off = 0;
  for (size_t r = 0; r < n; ++r) {
    const std::string& s = reads[r].seq;
    const int len = static_cast<int>(s.size());
    qoff_[r] = off;
    qlen_[r] = len;
    uint8_t* fwd = codes_.data() + off;
    encode(s.data(), s.size(), fwd);
    revcomp_codes(fwd, s.size(), fwd + len);
    off += 2 * s.size();
    cbeg_[r] = static_cast<uint32_t>(cands_.size());
    seed_and_chain(fwd, len, static_cast<uint32_t>(r));
  }
  cbeg_[n] = static_cast<uint32_t>(cands_.size());
  stats_.candidates += cands_.size();

  // Reference windows: zero-copy pointers when inside the contig, else an N-padded copy.
  const int bw = p_.sc.bw, W = p_.sc.band_cols();
  size_t pad_bytes = 0;
  for (const Cand& c : cands_) {
    const int wl = qlen_[c.read] + 2 * bw;
    if (!ref_.inside(c.contig, c.wstart, wl)) pad_bytes += wl;
  }
  pad_.resize(pad_bytes);
  size_t po = 0;
  for (Cand& c : cands_) {
    const int wl = qlen_[c.read] + 2 * bw;
    if (ref_.inside(c.contig, c.wstart, wl)) {
      c.t = ref_.seq.data() + c.wstart;
    } else {
      ref_.fetch(c.contig, c.wstart, wl, pad_.data() + po);
      c.t = pad_.data() + po;
      po += wl;
    }
  }
  stats_.t_seed += t_seed.sec();

  Timer t_align;
  // 3. Rank the candidates of multi-hit reads with the score-only kernel.
  tasks_.clear();
  task_ref_.clear();
  for (size_t r = 0; r < n; ++r) {
    if (cbeg_[r + 1] - cbeg_[r] < 2) continue;
    for (uint32_t c = cbeg_[r]; c < cbeg_[r + 1]; ++c) {
      tasks_.push_back(task_for(cands_[c]));
      task_ref_.push_back(c);
    }
  }
  if (res_.size() < tasks_.size()) res_.resize(tasks_.size());
  aligner_->align(tasks_.data(), tasks_.size(), false, res_.data());
  for (size_t i = 0; i < tasks_.size(); ++i) {
    cands_[task_ref_[i]].sw = res_[i].score;
    stats_.cells += static_cast<uint64_t>(tasks_[i].qlen) * W;
  }

  chosen_.assign(n, -1);
  second_.assign(n, 0);
  for (size_t r = 0; r < n; ++r) {
    const uint32_t b = cbeg_[r], e = cbeg_[r + 1];
    if (b == e) continue;
    if (e - b == 1) {
      chosen_[r] = static_cast<int32_t>(b);
      continue;
    }
    int s1 = -1, ties = 0;
    for (uint32_t c = b; c < e; ++c) {
      if (cands_[c].sw > s1) {
        s1 = cands_[c].sw;
        ties = 1;
      } else if (cands_[c].sw == s1) {
        ++ties;
      }
    }
    int s2 = 0;
    if (ties > 1) s2 = s1;
    else
      for (uint32_t c = b; c < e; ++c)
        if (cands_[c].sw != s1) s2 = std::max(s2, cands_[c].sw);
    // Equal-best hits: pick one pseudo-randomly (but reproducibly) to avoid positional bias.
    uint32_t pick = fnv1a(reads[r].name) % static_cast<uint32_t>(ties);
    for (uint32_t c = b; c < e; ++c) {
      if (cands_[c].sw == s1 && pick-- == 0) {
        chosen_[r] = static_cast<int32_t>(c);
        break;
      }
    }
    second_[r] = s2;
  }

  // 4. Full alignment (with traceback) of each read's best candidate.
  tasks_.clear();
  task_ref_.clear();
  for (size_t r = 0; r < n; ++r) {
    if (chosen_[r] < 0) continue;
    tasks_.push_back(task_for(cands_[chosen_[r]]));
    task_ref_.push_back(static_cast<uint32_t>(r));
  }
  if (res_.size() < tasks_.size()) res_.resize(tasks_.size());
  aligner_->align(tasks_.data(), tasks_.size(), true, res_.data());
  for (const AlignTask& t : tasks_) stats_.cells += static_cast<uint64_t>(t.qlen) * W;
  stats_.t_align += t_align.sec();

  // 5. Records + MAPQ.
  for (size_t r = 0; r < n; ++r) {
    MapRecord& m = out[r];
    m.contig = -1;
    m.pos = 0;
    m.rev = m.mapq = 0;
    m.score = m.nm = 0;
    m.n_cand = static_cast<int32_t>(cbeg_[r + 1] - cbeg_[r]);
    m.cigar.clear();
  }
  for (size_t i = 0; i < tasks_.size(); ++i) {
    const uint32_t r = task_ref_[i];
    const AlignResult& a = res_[i];
    if (a.score < p_.min_aln_score) continue;
    const Cand& c = cands_[chosen_[r]];
    MapRecord& m = out[r];
    m.contig = c.contig;
    m.rev = static_cast<uint8_t>(c.rev);
    m.pos = c.wstart + a.tb - static_cast<int64_t>(ref_.contigs[c.contig].offset);
    m.score = a.score;
    m.nm = a.nm;
    if (a.qb > 0) m.cigar.push_back(static_cast<uint32_t>(a.qb) << 4 | kCigS);
    m.cigar.insert(m.cigar.end(), a.cigar.begin(), a.cigar.end());
    if (qlen_[r] > a.qe) m.cigar.push_back(static_cast<uint32_t>(qlen_[r] - a.qe) << 4 | kCigS);
    if (m.n_cand <= 1) {
      m.mapq = 60;
    } else {
      const int s2 = second_[r];
      const double q = s2 >= a.score ? 0.0 : 6.02 * (a.score - s2) / p_.sc.match + 0.499;
      m.mapq = static_cast<uint8_t>(std::min(60.0, q));
    }
    ++stats_.mapped;
  }
  stats_.reads += n;
}

void Mapper::write_sam_header(FILE* fp, const std::string& cmdline) const {
  std::fprintf(fp, "@HD\tVN:1.6\tSO:unsorted\n");
  for (const Contig& c : ref_.contigs) std::fprintf(fp, "@SQ\tSN:%s\tLN:%u\n", c.name.c_str(), c.len);
  std::fprintf(fp, "@PG\tID:dnaln\tPN:dnaln\tVN:%s\tCL:%s\n", kVersion, cmdline.c_str());
}

void Mapper::append_sam(const SeqRecord& r, const MapRecord& m, std::string& o) const {
  o += r.name;
  if (!m.mapped()) {
    o += "\t4\t*\t0\t0\t*\t*\t0\t0\t";
    o += r.seq.empty() ? "*" : r.seq;
    o += '\t';
    o += r.qual.empty() ? "*" : r.qual;
    o += '\n';
    return;
  }
  o += '\t';
  append_int(o, m.rev ? 16 : 0);
  o += '\t';
  o += ref_.contigs[m.contig].name;
  o += '\t';
  append_int(o, m.pos + 1);
  o += '\t';
  append_int(o, m.mapq);
  o += '\t';
  for (uint32_t c : m.cigar) {
    append_int(o, c >> 4);
    o += kCigOps[c & 0xf];
  }
  o += "\t*\t0\t0\t";
  if (m.rev) {
    revcomp_append(r.seq, o);
    o += '\t';
    if (r.qual.empty()) o += '*';
    else o.append(r.qual.rbegin(), r.qual.rend());
  } else {
    o += r.seq;
    o += '\t';
    o += r.qual.empty() ? "*" : r.qual;
  }
  o += "\tNM:i:";
  append_int(o, m.nm);
  o += "\tAS:i:";
  append_int(o, m.score);
  o += '\n';
}

}  // namespace dnaln
