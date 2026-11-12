# Tools

- `gen_configs.py OUTDIR`: writes the reference config set (groups `structured`, `triangle`, `tetra`).
- `run_writer_on_configs.sh WRITER_JL OUTDIR`: generates the configs and runs `writer.jl` on them.
  The folder is a valid `WRITER_REF_DIR` for `test_writer.jl`.
- `compare_h5.jl REF NEW [--ref-column-major]`: bit-by-bit comparison of the `.h5` files of two folders.
  Use the flag only for a reference written before schema version 1 (column-major matrices).
- `compare_schema_py.py JULIA_DIR PYTHON_DIR`: compares every dataset and attribute of the Julia and the
  Python output of the same configs, up to edge order. It needs h5py. Known differences, both by
  convention: the edge order (Julia: first appearance, Python: sorted) and the sign of the 3-D cell
  orientation (Julia positive, Python negative).
- `check_style.py`, `dup_scan.py`, `check_ste.py`: checks of `writer.jl` (house style of the kernels,
  duplicate names, long sentences in comments). Usage: `python3 tools/check_style.py writer.jl`.
