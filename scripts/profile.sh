#!/usr/bin/env bash
# Hardware-counter / cache profiling of the alignment kernel and the mapper.
#   Linux : perf stat (cycles, IPC, L1/LLC misses), else valgrind cachegrind
#   macOS : records an Instruments trace with xctrace (open it in Instruments.app)
# The built-in, portable cache experiment is section 4 of `dnaln bench`.
set -euo pipefail
cd "$(dirname "$0")/.."
REF="${REF:-data/ecoli.fa}"
[[ -s "$REF" ]] || REF=data/synthetic.fa
[[ -s "$REF" ]] || { echo "No reference: run scripts/get_ecoli.sh or SYNTHETIC=1 scripts/run_benchmark.sh first" >&2; exit 1; }
NAME="$(basename "${REF%.*}")"
SIM="data/${NAME}_sim150"
make -j4 >/dev/null
[[ -s "$SIM.fq" ]] || build/dnaln simulate "$REF" -n 100000 -o "$SIM"
mkdir -p results
CMD=(build/dnaln bench "$REF" "$SIM.fq" "$SIM.truth.tsv" --kernel-only --reps 1)

if command -v perf >/dev/null 2>&1; then
  perf stat -e cycles,instructions,L1-dcache-loads,L1-dcache-load-misses,LLC-loads,LLC-load-misses,branch-misses "${CMD[@]}"
elif command -v valgrind >/dev/null 2>&1; then
  valgrind --tool=cachegrind --cache-sim=yes --cachegrind-out-file=results/cachegrind.out "${CMD[@]}" --limit 20000
  cg_annotate results/cachegrind.out | head -60
elif [[ "$(uname)" == Darwin ]] && xcrun xctrace version >/dev/null 2>&1; then
  TEMPLATE="${TEMPLATE:-CPU Counters}"   # or "Time Profiler"
  rm -rf results/dnaln.trace
  xcrun xctrace record --template "$TEMPLATE" --output results/dnaln.trace --launch -- "${CMD[@]}"
  echo "Open results/dnaln.trace in Instruments to inspect counters / hot spots."
else
  echo "No profiler found (perf, valgrind, or Xcode's xctrace). Running the built-in cache sweep:"
  build/dnaln bench "$REF" "$SIM.fq" "$SIM.truth.tsv" --skip-map --limit 20000
fi
