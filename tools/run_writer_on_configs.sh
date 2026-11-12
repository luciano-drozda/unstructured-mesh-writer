#!/bin/bash
# Generate the reference configs and run writer.jl on them.
# Usage: tools/run_writer_on_configs.sh WRITER_JL OUTDIR
# OUTDIR gets the groups structured/, triangle/ and tetra/, each with one
# .yaml, .h5 and .xdmf file per config, and one log file per group.
# Use OUTDIR as WRITER_REF_DIR for test_writer.jl, or compare two such
# folders with tools/compare_h5.jl.
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
W=$(realpath "$1"); OUT=$2
rm -rf "$OUT"; mkdir -p "$OUT"
python3 "$HERE/gen_configs.py" "$OUT"
for g in structured triangle tetra; do
    julia "$W" "$OUT/$g" > "$OUT/$g.log" 2>&1 || echo "$g: exit $?"
done
