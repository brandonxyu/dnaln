# dnaln

[![CI](https://github.com/brandonxyu/dnaln/actions/workflows/ci.yml/badge.svg)](https://github.com/brandonxyu/dnaln/actions/workflows/ci.yml)

**A fast DNA short-read aligner.** `dnaln` takes sequencing reads (FASTQ) and a reference
genome (FASTA), works out where in the genome each read came from, and writes the
alignments in the standard SAM format used by samtools, IGV and most other genomics tools.

* **Fast:** about 330,000 reads per second on one CPU core of an Apple M5 Pro
  (150 bp reads, E. coli genome)
* **Accurate:** places 98.70% of simulated reads correctly, vs 98.69% for minimap2
* **Simple:** one self-contained program, one command

## Supported platforms

| Platform | Prebuilt download | Build from source |
|---|---|---|
| macOS, Apple Silicon or Intel | ✓ | ✓ |
| Linux x86-64 | ✓ | ✓ |
| Linux ARM64 | ✓ | ✓ |

* **macOS:** one universal download for both chip types, built for macOS 11 or later. CI
  tests it natively on an Apple Silicon Mac and on an Intel Mac (macOS 15).
* **Linux:** the downloads are statically linked, so they don't depend on your system's
  glibc or other shared libraries at run time. CI runs them in minimal and older Linux
  userlands (BusyBox, Alpine, CentOS 7 and Ubuntu 18.04 containers) on x86-64 and ARM64.
* **From source:** every change is built and tested on Linux x86-64, Linux ARM64 and Apple
  Silicon Macs.
* **Windows:** not supported natively.

## Install

**Prebuilt (no compiler needed):**

```bash
curl -fsSL https://raw.githubusercontent.com/brandonxyu/dnaln/main/install.sh | sh
```

This downloads the latest release, verifies its SHA-256 checksum and puts `dnaln` in
`~/.local/bin`. If that folder isn't on your PATH yet, it prints the line to add.

**From source** needs a C++17 compiler, `make` and zlib:

* macOS: `xcode-select --install`
* Ubuntu/Debian: `sudo apt install build-essential zlib1g-dev`
* Fedora/RHEL: `sudo dnf install gcc-c++ make zlib-devel`

```bash
git clone https://github.com/brandonxyu/dnaln.git
cd dnaln
make
make install PREFIX="$HOME/.local"   # or: sudo make install   (installs to /usr/local/bin)
```

`$HOME/.local/bin` is often not on your PATH (it isn't by default on macOS). If
`dnaln --version` says "command not found", add it once and open a new terminal:

```bash
echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.zshrc   # bash users: ~/.bashrc
```

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

### Example output

Mapping 100,000 reads to the E. coli genome (Apple M5 Pro):

```
$ dnaln map ecoli.fa reads.fq -o reads.sam
[dnaln] reference: 1 contig(s), 4641652 bp (loaded in 0.00 s)
[dnaln] index: k=15 w=10, 838542 distinct minimizers, 867198 positions, 35.3 MB, built in 0.07 s
[dnaln] alignment kernel: NEON x16 lanes
[dnaln] mapped 100000 / 100000 reads (100.00%) in 0.311 s -> 321291 reads/s single-threaded (excl. index build)
[dnaln] time split: seed+chain 0.163 s, alignment 0.115 s (4.95 GCUPS), parse+output 0.033 s; total wall 0.38 s
```

Each line of `reads.sam` describes one read (sequence and qualities shortened here):

```
sim_10  16  NC_000913.3  4370989  60  144M1D6M  *  0  0  ACGACACTGAACATACGAAT...  IIIIIIIIII...  NM:i:2  AS:i:288
```

That is: read `sim_10` maps to the reverse strand (`16`) of `NC_000913.3` at position
4,370,989 with mapping quality 60 (confident). It aligns as 144 matching positions, one
deleted base and 6 more (`144M1D6M`), with 2 differences from the genome (`NM:i:2`).

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

**Good fit:** short reads (roughly 50–1,000 bp, e.g. Illumina) against bacterial, viral and
other small genomes. It has been benchmarked on E. coli (4.6 Mbp).

**Larger genomes** up to 2.1 billion bases are supported but not yet benchmarked. Mapping to
E. coli peaks at 95 MB of memory, and most of that grows with genome size, at roughly
10–15 bytes per base. That estimate gives about 1.2–1.8 GB for a 120 Mbp genome such as
fruit fly or Arabidopsis.

**Not yet supported:**
* Genomes larger than 2.1 billion bases, such as human (positions are stored in 32 bits).
* Long reads (Oxford Nanopore, PacBio): use [minimap2](https://github.com/lh3/minimap2).
* Paired-end awareness: you can map each file of a pair separately, but mate information is
  not used or reported.
* Multiple CPU threads: dnaln uses one core.
* Secondary and supplementary alignments: only the best placement of each read is reported.

## Performance

Measured on one core of an Apple M5 Pro: 100,000 simulated 150 bp reads (1% substitutions,
0.1% indels) against E. coli K-12, compared with minimap2 2.31 (`-ax sr -t1`):

| | dnaln | minimap2 |
|---|---|---|
| Reads placed correctly | 98.70% | 98.69% |
| Wrong answers among high-confidence (MAPQ 60) reads | 0 of 97,821 | 1 of 97,611 |
| Total time, including index build | 0.38 s | 0.60 s |

The remaining 1.3% are reads from E. coli's identical repeated regions, which no aligner can
place uniquely; both tools mark them as ambiguous (MAPQ 0).

These are single runs on one machine. On GitHub's shared CI machines, the 20,000-read
end-to-end runs in the test jobs on Linux x86-64, Linux ARM64 and Apple Silicon measure
about 130,000–220,000 reads per second. The shorter smoke tests that check the release
binaries are for correctness only; their throughput is not a comparable benchmark. The
comparison with minimap2 covers single-end short reads on one core only:
minimap2 is a much more general tool (long reads, paired-end, spliced alignment,
multithreading). To reproduce, run `./scripts/run_benchmark.sh` (see the
[developer guide](docs/DEVELOPMENT.md)).

## How it works

dnaln finds short exact matches (*minimizer seeds*) between each read and the genome using a
hash index, links consistent seeds into candidate locations (*chaining*), and lines the read
up base by base at the best candidates with a banded Smith–Waterman alignment. The alignment
step uses SIMD vector instructions (NEON on ARM, SSE4.1/AVX2 on x86) to align 16 reads at
once: 6.7–9.4x faster than plain code on an Apple M5 Pro, with bit-identical results.

The [developer guide](docs/DEVELOPMENT.md) covers the design in depth, plus building,
testing, benchmarking, the live demo and publishing releases.

## License

MIT. See [LICENSE](LICENSE).
