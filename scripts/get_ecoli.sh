#!/usr/bin/env bash
# Downloads the E. coli K-12 substr. MG1655 reference genome (RefSeq NC_000913.3,
# 4,641,652 bp, ~4.7 MB uncompressed) from NCBI into data/ecoli.fa.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p data
OUT=data/ecoli.fa
if [[ -s "$OUT" ]]; then echo "$OUT already exists"; exit 0; fi

EFETCH="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000913.3&rettype=fasta&retmode=text"
FTP="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz"

echo "Downloading E. coli K-12 MG1655 (NC_000913.3) from NCBI..."
if curl -fsSL --retry 3 "$EFETCH" -o "$OUT.tmp" && head -c1 "$OUT.tmp" | grep -q '>'; then
  mv "$OUT.tmp" "$OUT"
else
  echo "E-utilities download failed; trying the NCBI FTP mirror..."
  curl -fsSL --retry 3 "$FTP" | gunzip -c > "$OUT.tmp"
  mv "$OUT.tmp" "$OUT"
fi
echo "Saved $OUT ($(grep -v '>' "$OUT" | tr -d '\n' | wc -c | tr -d ' ') bp)"
