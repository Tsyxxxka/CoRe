# Four representative graph queries

This is a **single-threaded** query package including PageRank, BFS, DFS and Graph Coloring.

## Build

```bash
bash querys/scripts/build.sh
```

## Run one query on one representation

```bash
querys/bin/bfs CoRe data/caida.bin results/final/summaries/caida_core \
    querys/output/caida_bfs_core.jsonl 3 1 -1

querys/bin/bfs Mags data/caida.bin /path/to/mags/caida_mags \
    querys/output/caida_bfs_mags.jsonl 3 1 -1

querys/bin/dfs CSR data/caida.bin data/caida.bin \
    querys/output/caida_dfs_csr.jsonl 3 1 -1

querys/bin/pagerank CoRe data/caida.bin results/final/summaries/caida_core \
    querys/output/caida_pr_core.jsonl 3 1 -1

querys/bin/coloring CoRe data/caida.bin results/final/summaries/caida_core \
    querys/output/caida_coloring_core.jsonl 3 1 -1
```

`MODE` is `CSR`, `Mags` or `CoRe`; subsequent arguments are graph file,
summary prefix (ignored for CSR), output JSONL, timed repetitions (default 3),
warmups (default 1), BFS/DFS source ID (default `-1`: **seeded random vertex** in the
**original graph**, shared by CSR, Mags, and CoRe). 
