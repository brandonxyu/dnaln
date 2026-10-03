#!/usr/bin/env bash
# End-to-end benchmark:
#   build -> simulate reads with ground truth -> benchmark harness (kernel speedup,
#   correctness, throughput, cache sweep) -> map with dnaln -> map with minimap2 (if
#   available) -> score both against the ground truth side by side.
#
# Environment overrides:
#   REF=path.fa      reference (default data/ecoli.fa; see scripts/get_ecoli.sh)
#   SYNTHETIC=1      use a generated E. coli-like genome if REF is missing (offline)
#   N_READS=100000   READ_LEN=150   SUB=0.01   INDEL=0.001   SEED=42
#   MM2=/path/to/minimap2
set -euo pipefail
cd "$(dirname "$0")/.."

N_READS="${N_READS:-100000}"
READ_LEN="${READ_LEN:-150}"
SUB="${SUB:-0.01}"
INDEL="${INDEL:-0.001}"
SEED="${SEED:-42}"
REF="${REF:-data/ecoli.fa}"
OUT=results
mkdir -p data "$OUT"

echo "==> Building (release)"
make -j"$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)" >/dev/null

if [[ ! -s "$REF" ]]; then
  if [[ "${SYNTHETIC:-0}" == 1 ]]; then
    REF=data/synthetic.fa
    [[ -s "$REF" ]] || build/dnaln randref -o "$REF"
  else
    echo "Reference '$REF' not found. Run ./scripts/get_ecoli.sh first," >&2
    echo "or run 'SYNTHETIC=1 $0' to use a generated E. coli-like genome." >&2
    exit 1
  fi
fi
NAME="$(basename "${REF%.*}")"
SIM="data/${NAME}_sim${READ_LEN}"

echo "==> Simulating $N_READS x ${READ_LEN} bp reads from $REF (sub=$SUB indel=$INDEL seed=$SEED)"
build/dnaln simulate "$REF" -n "$N_READS" -l "$READ_LEN" -e "$SUB" -i "$INDEL" -s "$SEED" -o "$SIM"

echo
echo "==> Benchmark harness"
build/dnaln bench "$REF" "$SIM.fq" "$SIM.truth.tsv" | tee "$OUT/bench_${NAME}.txt"

echo "==> dnaln map (end-to-end wall time, including index build and file I/O)"
build/dnaln map "$REF" "$SIM.fq" -o "$OUT/${NAME}.dnaln.sam"

MM2="${MM2:-}"
if [[ -z "$MM2" ]]; then
  if command -v minimap2 >/dev/null 2>&1; then MM2="$(command -v minimap2)"
  elif [[ -x third_party/minimap2/minimap2 ]]; then MM2=third_party/minimap2/minimap2; fi
fi

echo
if [[ -n "$MM2" ]]; then
  echo "==> minimap2 $("$MM2" --version) -ax sr -t1"
  "$MM2" -t1 -ax sr "$REF" "$SIM.fq" > "$OUT/${NAME}.minimap2.sam" 2> "$OUT/${NAME}.minimap2.log"
  grep -E "Real time" "$OUT/${NAME}.minimap2.log" || true
  echo
  echo "==> Accuracy vs ground truth (last column is the reference for the 2% check)"
  build/dnaln eval "$SIM.truth.tsv" "$OUT/${NAME}.dnaln.sam" "$OUT/${NAME}.minimap2.sam" | tee "$OUT/accuracy_${NAME}.txt"
else
  echo "minimap2 not found - run ./scripts/get_minimap2.sh (or set MM2=/path/to/minimap2)"
  echo "for the head-to-head accuracy comparison. Scoring dnaln alone:"
  echo
  build/dnaln eval "$SIM.truth.tsv" "$OUT/${NAME}.dnaln.sam" | tee "$OUT/accuracy_${NAME}.txt"
fi
echo
echo "Reports written to $OUT/"
