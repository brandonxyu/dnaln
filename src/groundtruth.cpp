#include "groundtruth.hpp"

#include <algorithm>
#include <string>
#include <unordered_map>

#include "common.hpp"

namespace dnaln {

GtSet build_gt(const Reference& ref, const std::vector<SeqRecord>& reads, const TruthMap& truth, int bw,
               size_t limit) {
  GtSet g;
  std::unordered_map<std::string, int> cidx;
  for (size_t i = 0; i < ref.contigs.size(); ++i) cidx[ref.contigs[i].name] = static_cast<int>(i);
  std::vector<uint8_t> fwd;
  for (const SeqRecord& r : reads) {
    if (g.qlen.size() >= limit) break;
    const auto it = truth.find(r.name);
    if (it == truth.end()) continue;
    const auto ci = cidx.find(it->second.contig);
    if (ci == cidx.end()) continue;
    const int len = static_cast<int>(r.seq.size());
    fwd.resize(static_cast<size_t>(len));
    encode(r.seq.data(), r.seq.size(), fwd.data());
    const size_t qo = g.q.size();
    g.q.resize(qo + len);
    if (it->second.rev) revcomp_codes(fwd.data(), len, g.q.data() + qo);
    else std::copy(fwd.begin(), fwd.end(), g.q.begin() + static_cast<long>(qo));
    const int wl = len + 2 * bw;
    const size_t to = g.t.size();
    g.t.resize(to + wl);
    ref.fetch(ci->second, static_cast<int64_t>(ref.contigs[ci->second].offset) + it->second.start - bw, wl,
              g.t.data() + to);
    g.qoff.push_back(qo);
    g.toff.push_back(to);
    g.qlen.push_back(len);
    g.truth.push_back(&it->second);
  }
  g.tasks.resize(g.qlen.size());
  for (size_t i = 0; i < g.tasks.size(); ++i)
    g.tasks[i] = {g.q.data() + g.qoff[i], g.qlen[i], g.t.data() + g.toff[i], g.qlen[i] + 2 * bw};
  return g;
}

}  // namespace dnaln
