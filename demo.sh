#!/usr/bin/env bash
# One-command live demo for presenting dnaln.
#   ./demo.sh              step by step (press Enter to advance)
#   ./demo.sh --no-pause   run straight through
set -euo pipefail
cd "$(dirname "$0")"

echo "Preparing the demo..."
make -s -j8 >/dev/null

REF=data/ecoli.fa
NAME="E. coli K-12 MG1655"
if [[ ! -s "$REF" ]]; then
  ./scripts/get_ecoli.sh >/dev/null 2>&1 || true
fi
if [[ ! -s "$REF" ]]; then  # offline: fall back to a generated E. coli-like genome
  REF=data/synthetic.fa
  NAME="synthetic E. coli-like genome"
  [[ -s "$REF" ]] || build/dnaln randref -o "$REF" 2>/dev/null
fi

READS="data/demo_$(basename "${REF%.*}")"
[[ -s "$READS.fq" ]] || build/dnaln simulate "$REF" -n 100000 -l 150 -o "$READS" 2>/dev/null

mkdir -p results
ARGS=(--genome-name "$NAME" --sam results/demo.dnaln.sam)
MM2=""
if command -v minimap2 >/dev/null 2>&1; then MM2="$(command -v minimap2)"
elif [[ -x third_party/minimap2/minimap2 ]]; then MM2=third_party/minimap2/minimap2; fi
if [[ -n "$MM2" ]]; then
  "$MM2" -t1 -ax sr "$REF" "$READS.fq" > results/demo.minimap2.sam 2> results/demo.minimap2.log
  MM2_TIME="$(sed -n 's/.*Real time: \([0-9.]*\) sec.*/\1/p' results/demo.minimap2.log)"
  ARGS+=(--mm2-sam results/demo.minimap2.sam)
  [[ -n "$MM2_TIME" ]] && ARGS+=(--mm2-time "$MM2_TIME")
fi

exec build/dnaln demo "${ARGS[@]}" "$@" "$REF" "$READS.fq" "$READS.truth.tsv"
