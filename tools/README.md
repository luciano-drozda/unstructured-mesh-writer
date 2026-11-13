# Tools

- `gen_configs.py OUTDIR`: writes the reference config set (groups `structured`, `triangle`, `tetra`, `hierarchy`).
  The group `hierarchy` has 8 configs that cover every mesh type, with 3 or 4 levels. Each one gives one `.h5` file for each level.
- `run_writer_on_configs.sh WRITER_JL OUTDIR`: generates the configs and runs `writer.jl` on them.
  The folder is a valid `WRITER_REF_DIR` for `test_writer.jl`.
- `compare_h5.jl REF NEW [--ref-column-major]`: bit-by-bit comparison of the `.h5` files of two folders.
  Use the flag only for a reference written before schema version 1 (column-major matrices).
- `compare_schema_py.py JULIA_DIR PYTHON_DIR`: compares every dataset and attribute of the Julia and the
  Python output of the same configs, up to edge order. It needs h5py. The one known difference, by
  convention, is the edge order (Julia: first appearance, Python: sorted). The structured meshes agree
  in every other dataset and attribute.
- `check_style.py`, `dup_scan.py`, `check_ste.py`: checks of `writer.jl` (house style of the kernels,
  duplicate names, long sentences in comments). Usage: `python3 tools/check_style.py writer.jl`.
- `quality_stats.py FOLDER`: dihedral statistics of every tetrahedral `.h5` file in a folder. It works for the files of both writers.
- `hierarchy_multigrid_check.jl`: reads the level files like a solver and checks the hierarchy with a small multigrid code. It builds the prolongation matrix from `fine_vertex_parents`. It checks that the matrix copies linear functions exactly and that the Galerkin product equals the stiffness matrix of the coarse mesh. Then it runs V(2,2)-cycles for a Poisson problem and prints the convergence factor. Usage: `julia --project=. tools/hierarchy_multigrid_check.jl [NAME ...]`. Without a name, the script runs all nine (about 4 minutes). A periodic direction is allowed: the vertices that are periodic images of each other are one unknown. The names are `triangle_circle`, `gridsplit_stretched`, `kuhn_box`, `tetra_sphere`, `triangle_periodic`, `gridsplit_periodic_x`, `kuhn_periodic`, `kuhn_fully_periodic` and `tetra_periodic`.
- `hierarchy_quality.jl`: compares the diagonal rules of the 3-D hierarchy refinement (max-min and the 3 fixed diagonals) on some meshes. For each level it prints the dihedral angles, the radius-edge ratio and the mean-ratio shape. Usage: `julia --project=. tools/hierarchy_quality.jl`. A compare script for `writer.py` is not needed: `writer.py` has no hierarchy, so `compare_schema_py.py` skips the level files.
