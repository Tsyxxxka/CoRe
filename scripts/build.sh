#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
mkdir -p "$ROOT/bin"
CXX="${CXX:-g++}"
read -ra FLAGS <<< "${CXXFLAGS:--O3 -DNDEBUG -std=c++20 -march=native}"
"$CXX" "${FLAGS[@]}" -I"$ROOT/src" "$ROOT/src/core.cpp" -o "$ROOT/bin/core"
echo "Built $ROOT/bin/core"
