# dnaln developer guide

Everything for working on dnaln itself: building, testing, benchmarking, the live demo,
how the aligner works internally, and how to publish a release. For using dnaln, see the
[README](../README.md).

## Building and testing in VS Code

1. **Prerequisites (macOS):** the Xcode Command Line Tools, which provide `clang++`, `make`,
   zlib and `git`. Install them with `xcode-select --install` if needed.
   On Linux: `g++` or `clang++`, `make` and `zlib1g-dev`.
2. **Open the folder** in VS Code (*File → Open Folder…*). Accept the recommended
   extensions: *C/C++* (IntelliSense) and *CodeLLDB* (debugging).
3. **Build:** press **⇧⌘B**, which runs the default task *build (release)* (`make`).
4. **Run tasks** with *⌘⇧P → Tasks: Run Task*:
   * `run unit tests`: 11 tests, about 75k checks
   * `download E. coli reference`: NC_000913.3 from NCBI (~4.7 MB)
   * `build minimap2 (for comparison)`: clones and builds minimap2 into `third_party/`
   * `full benchmark (E. coli)`: simulate, benchmark, map, then compare with minimap2
   * `full benchmark (offline synthetic genome)`: needs no downloads
5. **Debug:** pick *Debug unit tests* or *Debug dnaln map* in the Run & Debug panel (F5).
   This builds with AddressSanitizer and UndefinedBehaviorSanitizer first. *Debug dnaln map*
   reads `data/ecoli.fa` and `data/ecoli_sim150.fq`, which `./scripts/run_benchmark.sh`
   creates.

## Building and testing from the terminal

```bash
make -j8 && make test                    # build and run unit tests
make debug                               # ASan + UBSan build in build-debug/
./scripts/get_ecoli.sh                   # data/ecoli.fa (E. coli K-12 MG1655)
./scripts/get_minimap2.sh                # optional: minimap2 for the accuracy comparison
./scripts/run_benchmark.sh               # everything; reports land in results/
SYNTHETIC=1 ./scripts/run_benchmark.sh   # offline variant, no downloads
```

Individual commands:

```bash
build/dnaln simulate data/ecoli.fa -n 100000 -l 150 -o data/sim       # sim.fq + sim.truth.tsv
build/dnaln map data/ecoli.fa data/sim.fq -o results/dnaln.sam        # SAM output
build/dnaln eval data/sim.truth.tsv results/dnaln.sam results/mm2.sam # side-by-side accuracy
build/dnaln bench data/ecoli.fa data/sim.fq data/sim.truth.tsv        # benchmark harness
build/dnaln map --help                                                # all options
```

## Live demo (for showing it to someone)

One command, which works from the VS Code terminal (*Terminal → New Terminal*):

```bash
./demo.sh
```

Or in VS Code: *⌘⇧P → Tasks: Run Task → ▶ Run demo*. It walks through six steps in
plain language, pausing for **Enter** between each so you can talk over it (about 1 minute):

1. Load the E. coli genome and build the search index (0.07 s).
2. Introduce 100,000 simulated DNA reads whose true origin is known.
3. Place one read and show it lined up letter by letter against the genome, with
   sequencing errors and a missing letter highlighted.
4. Race plain C++ against SIMD on 500,000 alignments, with live bars (6.7x faster, same answers).
5. Map all 100,000 reads (about 320,000 reads/s), compared with minimap2's time.
6. Grade every answer against the truth: 98.70% correct vs minimap2's 98.69%, and 0 errors
   among confident calls.

`./demo.sh --no-pause` runs straight through. The demo downloads E. coli if it's missing
(or falls back to a synthetic genome offline) and includes minimap2 when it's installed.

## Benchmark results

Measured with `scripts/run_benchmark.sh` on an Apple M5 Pro (NEON, one thread) against the
real **E. coli K-12 MG1655** genome (NC_000913.3, 4,641,652 bp). Workload: 100,000 simulated
150 bp reads (1% substitutions, 0.1% indels), compared with **minimap2 2.31** (`-ax sr -t1`).

| Target | Result |
|---|---|
| ≥ 20,000 reads/s, single thread | **337,000 reads/s** (in-memory); 321,000 reads/s with file I/O |
| ≥ 6x SIMD speedup over scalar | **9.43x** score-only, **6.74x** with traceback |
| Correct on 100,000 ground-truth alignments | 100% bit-identical to scalar; 100% provably optimal |
| Within 2% of minimap2's accuracy | dnaln **98.70%** vs minimap2 **98.69%** correctly placed (+0.005 pp) |

Both tools misplace the same ~1.3% of reads. These come from E. coli's long exact repeats
(rRNA operons, IS elements), where a read has several equally good locations and both tools
correctly report MAPQ 0. Total wall time, including index build, was 0.38 s for dnaln and
0.60 s for minimap2. minimap2 places individual bases slightly better (99.95% vs 99.87%
base-level concordance) because it extends alignments to the read ends, where local
Smith–Waterman soft-clips a terminal mismatch.

```
metric                                    ecoli.dnaln.sam     ecoli.minimap2.sam
mapped                                           100.000%                99.991%
ACCURACY (correct / all reads)                    98.699%                98.694%
MAPQ>=60 mapped | error rate             97.82% | 0.0000%       97.61% | 0.0010%
MAPQ>=30 mapped | error rate             98.03% | 0.0000%       98.07% | 0.0020%
base-level concordance                            99.866%                99.945%

== 1. Banded Smith-Waterman kernel: scalar baseline vs SIMD ==
   mode        kernel              time (ms)     GCUPS   speedup    identical to scalar
   score-only  scalar                  557.0      0.89     1.00x                      -
   score-only  NEON x16 lanes           59.1      8.38     9.43x      100,000 / 100,000
   traceback   scalar                  682.8      0.72     1.00x                      -
   traceback   NEON x16 lanes          101.4      4.88     6.74x      100,000 / 100,000

== 2. Correctness on 100,000 ground-truth simulated alignments ==
   SIMD == scalar (score, coordinates, CIGAR, NM; both modes)     100,000 / 100,000  PASS
   CIGAR re-scores to the reported SW score                       100,000 / 100,000  PASS
   Optimality: SW score >= score of the true alignment path       100,000 / 100,000  PASS
   Read bases aligned to their true reference position                      99.866%
```

Accuracy holds across read lengths and error rates (20k reads each, on the synthetic
E. coli-like genome from `dnaln randref`):

| Reads | Accuracy | Error rate at MAPQ ≥ 60 | Throughput |
|---|---|---|---|
| 100 bp, 1% sub, 0.1% indel | 99.38% | 0% | 523k reads/s |
| 250 bp, 1% sub, 0.1% indel | 99.62% | 0% | 199k reads/s |
| 150 bp, 3% sub, 0.5% indel | 99.32% | 0% | 367k reads/s |
| 150 bp, 5% sub, 1% indel | 97.77% | 0% | 384k reads/s |

## How it works

```
reads ─► minimizers ─► index lookup ─► anchors ─► chaining ─► candidate loci
           (w,k)        (prefetched)              (DP)          (strand, diagonal)
                                                                     │
SAM ◄── MAPQ ◄── traceback SW on best ◄── score-only SW on all ◄─────┘
                 (CIGAR, NM)               candidates (multi-hit reads only)
```

**Minimizer index** (`minimizer.cpp`, `index.cpp`). Every k-mer (default k = 15) is reduced
to its canonical strand and hashed with an invertible integer hash. In each window of w = 10
consecutive k-mers, the k-mer with the smallest hash is the minimizer (about 2/(w+1) of
positions). The sliding-window minimum is a streamed **van Herk/Gil-Werman** block
prefix/suffix scan with no data-dependent branches. That was 2x faster than the classic
monotone deque, whose pops mispredict constantly. The index is an open-addressing hash
table (16-byte buckets, load factor ≤ 0.5) pointing into one contiguous occurrence array,
so a lookup touches about one cache line plus a sequential scan.

**Chaining** (`chain.cpp`). Seed hits become anchors in strand-oriented read coordinates.
A minimap2-style DP links anchors with a gap cost and extracts disjoint chains. Chains
scoring at least 50% of the best become candidate loci, deduplicated by diagonal.

**Banded Smith–Waterman** (`align.cpp`, `align_simd.cpp`). The kernel does local alignment
with affine gaps (match 2, mismatch −4, gap 4 + 2·len). The band is stored in diagonal
coordinates: row i covers columns j = i + d, d ∈ [0, 2·bw]. Then the diagonal and vertical
(E) predecessors both come from the previous row, the rows update in place, and only the
horizontal (F) term carries a dependency along the row.

**SIMD: inter-sequence vectorization.** Instead of vectorizing inside one alignment, which
needs shuffles for the F dependency, 16 alignments run in lock-step, one per 16-bit lane
(the approach of BWA-MEM2's banded SW). Each lane computes exactly the scalar recurrence,
which is why the results can be bit-identical. Key optimizations:

* **Byte-table scoring.** Sequences are transposed into lane-interleaved byte rows: query
  codes are pre-shifted, N is encoded as `0x80`, and `idx = q | t`. One NEON `TBX` (x86:
  `PSHUFB`) then scores 16 lanes, including the N penalty, replacing about five
  compare/select ops per 8 lanes.
* **Two registers per logical vector** (`kUnroll = 2` on NEON). This gives the out-of-order
  core two independent dependency chains through the latency-bound F→H path.
* **1 byte per cell of traceback.** Four direction flags are packed with NEON
  shift-insert (`vsri` + `vshrn`). There is no "H = 0" flag: the traceback re-derives the
  running score and stops when it reaches 0, which removes two ops per cell.
* **Two-phase alignment.** A score-only pass ranks the candidates of multi-hit reads, and
  only the winner gets the traceback pass.

**Cache design and profiling.** Per 16-lane group the hot working set is one band row of H/E
(about 2 KB), the transposed sequences (about 5 KB) and the traceback matrix (86 KB for 150
bp reads at bw = 16). All of it fits in the 128 KB L1 data cache of Apple M-series
performance cores; on CPUs with a 32–48 KB L1 the traceback matrix lives in L2. Scratch
buffers are 64-byte aligned and reused. Section 4 of `dnaln bench` profiles this: it sweeps
the band width so the working set grows past L1 into L2 and reports throughput at each
point. Index lookups use two-stage software prefetching (buckets, then occurrence lists).
That gains only about 2% on E. coli, because the 35 MB index mostly stays cache-resident; it
matters more for larger genomes. `scripts/profile.sh` collects hardware counters with `perf`
or cachegrind on Linux, or records an Instruments trace with `xctrace` on macOS.

**MAPQ.** A read with a single candidate gets 60. Otherwise MAPQ ≈ 6·(S1 − S2)/match,
capped at 60, where S1 and S2 are the best and second-best scores. Equal-best hits are
broken pseudo-randomly with a read-name hash, as minimap2 does, to avoid positional bias.

## Validation methodology

* **Simulator** (`dnaln simulate`). Samples reads uniformly from both strands with
  substitutions and short indels, and records the true CIGAR in `*.truth.tsv`.
* **Kernel correctness** (`dnaln bench`, section 2). One banded alignment per read at its
  true locus (100,000 alignments). The harness checks four things:
  1. SIMD and scalar outputs are identical (score, coordinates, CIGAR, NM).
  2. Every CIGAR re-scores to the reported score.
  3. **Optimality:** the SW score is at least the score of the true alignment path. A local
     optimum must beat any in-band path, including the truth.
  4. Base-level agreement with the true alignment.
* **Mapping accuracy** (`dnaln eval`). A read is correct when it is on the true strand and
  its unclipped start is within 10 bp of the truth. The report also gives error rate by
  MAPQ threshold and base-level concordance. It parses any SAM, so dnaln and minimap2 are
  scored by identical code.
* **Unit tests** (`make test`):
  * minimizers vs a brute-force reference
  * the scalar kernel vs a naive full-matrix DP
  * SIMD vs scalar on randomized tasks: four scoring schemes, N bases, ragged lane groups,
    padded targets
  * traceback invariants
  * an end-to-end mapping test

  CI runs them on every push on Linux x86-64 (GCC, AVX2), Linux ARM64 (GCC, NEON) and
  macOS (Clang, NEON). They also pass with the SSE4.1 kernel (on Apple Silicon,
  `make ARCH="-arch x86_64 -msse4.1"` builds an x86 version that runs under Rosetta) and
  with the portable fallback (`make ARCH=-DDNALN_NO_SIMD`), and run clean under ASan and
  UBSan (`make debug`).

## Releasing a new version

1. Update `kVersion` in `src/common.hpp`, commit and push to `main`.
2. Wait for CI to pass. Every push already builds the release binaries once
   (`make PORTABLE=1`) and then tests those exact files:
   * Linux x86-64 and ARM64: statically linked, so there is no runtime dependency on the
     system's glibc. Each runs in BusyBox, Alpine, CentOS 7 and Ubuntu 18.04 containers.
   * macOS: one universal binary (Apple Silicon + Intel, built for macOS 11+), run natively
     on an Apple Silicon runner and on an Intel runner (`macos-15-intel`). Each must report
     the expected SIMD kernel (NEON or SSE4.1).

   They use baseline CPU features, so SSE4.1 rather than AVX2 on x86. You can download them
   from the run's *Artifacts* section.
3. Tag and push the tag: `git tag v0.1.0 && git push origin v0.1.0`. CI checks that the tag
   matches `dnaln --version`, then the `publish` job (the only job with write access)
   creates the GitHub Release with the binaries and a `SHA256SUMS` file.
4. `install.sh` always downloads the latest release and verifies its checksum, so users
   get the new version automatically.

## Project layout

```
src/
  common.hpp/.cpp      encoding, CIGAR helpers, RNG, aligned buffers
  seqio.*              FASTA/FASTQ reader (plain or gzip)
  reference.*          concatenated 2-bit-coded reference + contig table
  minimizer.*          (w,k)-minimizer sketching (branchless sliding-window min)
  index.*              minimizer hash index
  chain.*              colinear chaining DP
  simd.hpp             NEON / SSE4.1 / AVX2 / portable SIMD abstraction
  align.hpp/.cpp       banded SW API, scalar baseline kernel, traceback
  align_simd.cpp       inter-sequence SIMD kernel
  mapper.*             mapping pipeline, MAPQ, SAM output
  simulate.*           read simulator + synthetic genome
  evaluate.*           accuracy evaluation against ground truth
  bench.cpp            benchmark harness
  demo.cpp             narrated live demo for presentations (run via ./demo.sh)
  groundtruth.*        alignment tasks at each simulated read's true locus
  main.cpp, cli.hpp    command-line interface
tests/test_main.cpp    unit + integration tests
scripts/               get_ecoli.sh, get_minimap2.sh, run_benchmark.sh, profile.sh
demo.sh                one-command narrated demo
install.sh             installer for prebuilt release binaries
.github/workflows/     CI: tests and release builds on every push; GitHub Release on version tags
docs/                  this guide
.vscode/               build/test/benchmark tasks, debug configs, IntelliSense settings
Makefile               primary build (CMakeLists.txt provided as an alternative)
```
