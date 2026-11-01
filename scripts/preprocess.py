#!/usr/bin/env python3
"""Make a Mags-compatible undirected simple graph.

Input:
  whitespace-delimited edge list; plain text or .gz.
  Lines beginning with '#' are comments.
  If a record has more than two columns (e.g. CAIDA relationship/weight),
  only the first two endpoint columns are used.

Cleaning:
  1) remove self-loops;
  2) ignore edge direction by sorting the two endpoints;
  3) remove duplicate/reverse-duplicate edges;
  4) retain only vertices incident to at least one cleaned edge.

Binary:
  little-endian int64 n,m; then m pairs of little-endian int32 dense IDs.
Dense IDs follow sorted surviving original IDs, identically for every algorithm.
"""
import argparse, gzip, hashlib, json, time
from pathlib import Path
import numpy as np


def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def preprocess(src, prefix):
    t = time.perf_counter()
    src = Path(src)
    prefix = Path(prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)

    op = gzip.open if src.suffix == ".gz" else open
    with op(src, "rt") as f:
        # Important for CAIDA: the file can contain a third relationship/weight
        # column. Mags uses only graph topology, so read the two endpoint columns.
        a = np.loadtxt(
            f,
            dtype=np.int64,
            comments="#",
            ndmin=2,
            usecols=(0, 1),
        )

    if a.size == 0:
        a = np.empty((0, 2), dtype=np.int64)

    raw = len(a)
    loops = int(np.sum(a[:, 0] == a[:, 1]))

    # Remove self-loops first.
    b = a[a[:, 0] != a[:, 1]]
    directed_unique = np.unique(b, axis=0)
    same_direction = len(b) - len(directed_unique)
    flips = int(np.sum(b[:, 0] > b[:, 1]))

    # Drop direction, then deduplicate.
    b.sort(axis=1)
    clean = np.unique(b, axis=0)
    m = len(clean)

    # Mags Table 2 counts the graph after cleaning. Vertices that occur only in
    # removed self-loops therefore disappear. This is why ids must be computed
    # from `clean`, not from the raw input.
    ids = np.unique(clean)
    n = len(ids)

    if n > np.iinfo(np.int32).max:
        raise ValueError("number of surviving vertices exceeds int32")

    dense = np.searchsorted(ids, clean).astype("<i4")

    with open(str(prefix) + ".bin", "wb") as f:
        np.asarray([n, m], dtype="<i8").tofile(f)
        dense.tofile(f)

    ids.astype("<i8").tofile(str(prefix) + ".ids.bin")

    # Text is used only for portability and independent read-back.
    with gzip.GzipFile(str(prefix) + ".clean.txt.gz", "wb", mtime=0) as f:
        import io
        with io.TextIOWrapper(f) as s:
            s.write("# Undirected, unique, no self-loops; original vertex IDs\n")
            np.savetxt(s, clean, fmt="%d")

    with open(str(prefix) + ".bin", "rb") as f:
        nh, mh = np.fromfile(f, "<i8", 2)
        ed = np.fromfile(f, "<i4").reshape(-1, 2)

    assert nh == n and mh == m and len(ed) == m
    assert not m or np.all(ed[:, 0] < ed[:, 1])
    assert np.array_equal(ed, np.unique(ed, axis=0))
    assert np.array_equal(ids[ed], clean)

    out = dict(
        dataset=prefix.name,
        input_filename=src.name,
        input_sha256=digest(src),
        n=int(n),
        input_records=raw,
        self_loops_removed=loops,
        same_direction_duplicates_removed=same_direction,
        reverse_duplicates_removed=int(len(directed_unique) - m),
        records_reoriented=flips,
        m=int(m),
        remaining_self_loops=0,
        remaining_duplicates=0,
        id_order="sorted surviving original IDs",
        edges_sha256=hashlib.sha256(dense.tobytes()).hexdigest(),
        graph_file_sha256=digest(str(prefix) + ".bin"),
        preprocessing_seconds=time.perf_counter() - t,
        readback_verified=True,
    )
    Path(str(prefix) + ".preprocessing.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out), flush=True)
    return out


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("input")
    p.add_argument("output_prefix")
    a = p.parse_args()
    preprocess(a.input, a.output_prefix)
