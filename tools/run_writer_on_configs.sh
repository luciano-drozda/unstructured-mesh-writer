#!/bin/bash
# Generate the reference configs and run writer.jl on them.
# Usage: tools/run_writer_on_configs.sh WRITER_JL OUTDIR
# OUTDIR gets the groups structured/, triangle/, tetra/, hierarchy/ and step/, each
# with one .yaml file per config, the .h5 and .xdmf files, and one log file
# per group. A config with a hierarchy gives one .h5 and one .xdmf file
# for each level (name_L0, name_L1, and so on).
# Use OUTDIR as WRITER_REF_DIR for test_writer.jl, or compare two such
# folders with tools/compare_h5.jl.
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
W=$(realpath "$1"); OUT=$2
rm -rf "$OUT"; mkdir -p "$OUT"
python3 "$HERE/gen_configs.py" "$OUT"
for g in structured triangle tetra hierarchy step; do
    julia --project="$HERE/.." "$W" "$OUT/$g" > "$OUT/$g.log" 2>&1 || echo "$g: exit $?"
done
