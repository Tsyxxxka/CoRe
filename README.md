# CoRe: lossless graph summarization without negative corrections

This is the source code of our paper "Processing-Efficient Lossless Graph Summarization". All experiments are conducted on a Linux server with an AMD EPYC 7443 24-Core processor and 1,024 GB of RAM. The summarization and query-processing algorithms are implemented in single-threaded C++. The C++ implementations are compiled with g++ 8.5 using -O3. For query evaluation, we perform one warm-up run followed by multiple measured runs and report the average runtime. For BFS and DFS, we use the same randomly selected source vertex across methods in each run.

## Files

- `src/core.cpp`: CoRe positive-only summary constructor; candidate generation,
  initial lazy merging, positive-gain coarsening, bounded local repartitioning,
  vertex relocation, and independent summary encoding.
- `src/common.hpp`: graph loader, summary data structures, and encoder.
- `scripts/preprocess.py`: converts an edge-list dataset to the common undirected
  simple graph `.bin` format (requires NumPy).
- `querys/`: four query executables (BFS/DFS/PageRank/Coloring).

## Build

From this directory:

```bash
bash scripts/build.sh
```

Requires GCC or Clang supporting C++20 and the standard library (`std::span`,
`std::bit_width`). By default, produces `bin/core`.

## Input graph format

Require preprocess to remove self-loops, generating undirected simple graphs.

```bash
python3 scripts/preprocess.py --help
```

## Build one NCF summary

```bash
mkdir -p results/final/summaries
bin/core data/caida.bin positive results/final/summaries/caida_core
```

The output prefix creates four raw binary files:

```text
caida_core.partition.bin    # n int32 super-node identifiers
caida_core.se.bin           # (super-node A, super-node B) int32 pairs
caida_core.cp.bin           # (vertex u, vertex v) positive corrections
caida_core.cm.bin           # empty for CoRe (negative corrections)
```