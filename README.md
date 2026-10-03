# dnaln

**A fast DNA short-read aligner.** `dnaln` takes sequencing reads (FASTQ) and a reference
genome (FASTA), works out where in the genome each read came from, and writes the
alignments in the standard SAM format used by samtools, IGV and most other genomics tools.

* **Fast:** about 330,000 reads per second on a single CPU core (E. coli, 150 bp reads)
* **Accurate:** places 98.70% of simulated reads correctly, vs 98.69% for minimap2
* **Simple:** one self-contained program, one command

## Install

**macOS or Linux, prebuilt (no compiler needed):**

```bash
curl -fsSL https://raw.githubusercontent.com/brandonxyu/dnaln/main/install.sh | sh
```

This puts `dnaln` in `~/.local/bin` and tells you if that folder needs adding to your PATH.

**From source** (needs a C++17 compiler, `make` and zlib; on macOS run
`xcode-select --install` first, on Ubuntu/Debian `sudo apt install build-essential zlib1g-dev`):

```bash
git clone https://github.com/brandonxyu/dnaln.git
cd dnaln
make
make install PREFIX=~/.local     # or: sudo make install   (to /usr/local/bin)
```

Check it worked with `dnaln --version`.

## Usage

```bash
dnaln map genome.fa reads.fq > alignments.sam
```

* **Genome:** FASTA, one or many sequences, plain or gzipped (`.fa.gz`).
* **Reads:** FASTQ or FASTA, plain or gzipped (`.fq.gz`).
* **Output:** SAM on standard output (or `-o file.sam`). Progress and timing go to the
  terminal.

To get a sorted, indexed BAM file you can open in IGV (requires
[samtools](https://www.htslib.org/)):

```bash
dnaln map genome.fa reads.fq.gz | samtools sort -o reads.bam
samtools index reads.bam
```

Common options (all have sensible defaults; see `dnaln map --help`):

| Option | Meaning | Default |
|---|---|---|
| `-o FILE` | write SAM to FILE instead of standard output | stdout |
| `--bw INT` | band width: the largest insertion/deletion aligned through | 16 |
| `-k INT`, `-w INT` | seed (minimizer) length and window size | 15, 10 |
| `-A`, `-B`, `-O`, `-E` | match score, mismatch / gap-open / gap-extend penalties | 2, 4, 4, 2 |
| `--min-score INT` | report reads scoring below this as unmapped | 40 |

### Try it without your own data

dnaln can generate a test genome and simulated reads with known true positions, then grade
its own answers:

```bash
dnaln randref -o genome.fa                   # 4.6 Mbp E. coli-like test genome
dnaln simulate genome.fa -n 100000 -o sim    # sim.fq + sim.truth.tsv (the true answers)
dnaln map genome.fa sim.fq -o sim.sam
dnaln eval sim.truth.tsv sim.sam             # % of reads placed correctly
```

## Is it right for your data?

**Good fit:** short reads (roughly 50–1,000 bp, e.g. Illumina) against genomes up to about
2 billion bases: bacteria, viruses, yeast, worm, fly, Arabidopsis and so on.

**Not yet supported:**
* Genomes larger than 2.1 billion bases, such as human (positions are stored in 32 bits).
* Long reads (Oxford Nanopore, PacBio): use [minimap2](https://github.com/lh3/minimap2).
* Paired-end awareness: you can map each file of a pair separately, but mate information is
  not used or reported.
* Multiple CPU threads: dnaln uses one core.
* Secondary and supplementary alignments: only the best placement of each read is reported.

## Performance

Measured on an Apple M5 Pro, one thread, 100,000 simulated 150 bp reads against
E. coli K-12 (1% substitutions, 0.1% indels), compared with minimap2 2.31 (`-ax sr -t1`):

| | dnaln | minimap2 |
|---|---|---|
| Reads placed correctly | 98.70% | 98.69% |
| Wrong answers among high-confidence (MAPQ 60) reads | 0 | 1 |
| Total time, including index build | 0.38 s | 0.60 s |

The remaining 1.3% are reads from E. coli's identical repeated regions, which no aligner can
place uniquely; both tools mark them as ambiguous (MAPQ 0).

## How it works

dnaln finds short exact matches (*minimizer seeds*) between each read and the genome using a
hash index, links consistent seeds into candidate locations (*chaining*), and lines the read
up base by base at the best candidates with a banded Smith–Waterman alignment. The alignment
step uses SIMD vector instructions (NEON on ARM, SSE4.1/AVX2 on x86) to align 16 reads at
once, 6.7–9.4x faster than plain code with bit-identical results.

The [developer guide](docs/DEVELOPMENT.md) covers the design in depth, plus building,
testing, benchmarking, the live demo and publishing releases.

## License

MIT. See [LICENSE](LICENSE).
