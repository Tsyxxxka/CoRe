#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
mkdir -p "$ROOT/bin"
CXX="${CXX:-g++}"
read -ra FLAGS <<< "${CXXFLAGS:--O3 -DNDEBUG -std=c++17 -march=native}"
for task in bfs dfs pagerank coloring; do
  "$CXX" "${FLAGS[@]}" -I"$ROOT/src" "$ROOT/src/$task.cpp" -o "$ROOT/bin/$task"
  echo "Built $ROOT/bin/$task"
done
