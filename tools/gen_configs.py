# Generate the reference config set. Usage: gen_configs.py outdir
import itertools, os, sys
out = sys.argv[1]
def w(group, name, text):
    d = os.path.join(out, group); os.makedirs(d, exist_ok=True)
    open(os.path.join(d, name + ".yaml"), "w").write(text)
def yb(b): return "[" + ", ".join("true" if x else "false" for x in b) + "]"
# structured 2-D gridsplit
for diag, val in (("main", "backslash"), ("flip", "slash"), ("alternate", "alternate")):
    w("structured", f"gs_{diag}", f"mesh_type: gridsplit\nbox_length: [2.0, 1.0]\nn_divisions: [7, 5]\ndiagonal: {val}\n")
w("structured", "gs_stretch", "mesh_type: gridsplit\nbox_length: [2.0, 1.0]\nn_divisions: [7, 5]\nstretching_ratio: [1.2, 0.9]\n")
w("structured", "gs_sym", "mesh_type: gridsplit\nbox_length: [2.0, 1.0]\nn_divisions: [8, 7]\nstretching_ratio: [1.3, 1.1]\nstretching_symmetric: [true, true]\n")
for p in itertools.product([False, True], repeat=2):
    w("structured", "gs_per_" + "".join("t" if x else "f" for x in p), f"mesh_type: gridsplit\nbox_length: [2.0, 1.0]\nn_divisions: [6, 4]\nperiodic: {yb(p)}\n")
# default mesh_type
w("structured", "gs_default", "box_length: [1.0, 1.0]\nn_divisions: [4, 4]\n")
# equilateral
for ny in (6, 7):
    for py in (False, True):
        w("structured", f"eq_ny{ny}_{'p' if py else 'n'}", f"mesh_type: equilateral\nbox_length: [2.0, 1.0]\nn_divisions: [6, {ny}]\nperiodic: [false, {'true' if py else 'false'}]\n")
w("structured", "eq_perx", "mesh_type: equilateral\nbox_length: [2.0, 1.0]\nn_divisions: [6, 6]\nperiodic: [true, false]\n")
# kuhn
for diag in ("main", "flip_x", "flip_y", "flip_z", "alternate_layers"):
    w("structured", f"kuhn_{diag}", f"mesh_type: kuhn\nbox_length: [1.0, 1.0, 1.0]\nn_divisions: [3, 4, 3]\ndiagonal: {diag}\n")
w("structured", "kuhn_stretch", "mesh_type: kuhn\nbox_length: [1.0, 2.0, 1.0]\nn_divisions: [3, 4, 3]\nstretching_ratio: [1.2, 0.9, 1.1]\nstretching_symmetric: [false, true, false]\n")
for p in itertools.product([False, True], repeat=3):
    w("structured", "kuhn_per_" + "".join("t" if x else "f" for x in p), f"mesh_type: kuhn\nbox_length: [1.0, 1.0, 1.0]\nn_divisions: [3, 3, 3]\nperiodic: {yb(p)}\n")
# triangle
opt = "mesh_options:\n  mesh_size: 0.2\n  mesh_size_obstacle: 0.08\n  grading: 0.3\n"
circ = "  - shape: circle\n    center: [1.0, 0.5]\n    radius: 0.2\n"
rect = "  - shape: rectangle\n    center: [2.0, 0.5]\n    half_extents: [0.25, 0.15]\n"
w("triangle", "tri_none", "mesh_type: triangle\nbox_length: [2.0, 1.0]\nmesh_options:\n  mesh_size: 0.2\n")
w("triangle", "tri_circle", f"mesh_type: triangle\nbox_length: [2.0, 1.0]\nobstacles:\n{circ}{opt}")
w("triangle", "tri_rect", f"mesh_type: triangle\nbox_length: [3.0, 1.0]\nobstacles:\n{rect}{opt}")
w("triangle", "tri_several", f"mesh_type: triangle\nbox_length: [3.0, 1.0]\nobstacles:\n{circ}{rect}{opt}")
for p in itertools.product([False, True], repeat=2):
    w("triangle", "tri_per_" + "".join("t" if x else "f" for x in p), f"mesh_type: triangle\nbox_length: [3.0, 1.0]\nperiodic: {yb(p)}\nobstacles:\n{circ}{rect}{opt}")
# tetra
topt = "mesh_options:\n  mesh_size: 0.4\n  mesh_size_obstacle: 0.2\n  grading: 0.3\n"
sph = "  - shape: sphere\n    center: [0.5, 0.5, 0.5]\n    radius: 0.2\n"
bx = "  - shape: box\n    center: [1.4, 0.5, 0.5]\n    half_extents: [0.15, 0.2, 0.15]\n"
w("tetra", "tet_none", "mesh_type: tetra\nbox_length: [1.0, 1.0, 1.0]\nmesh_options:\n  mesh_size: 0.34\n")
w("tetra", "tet_box", f"mesh_type: tetra\nbox_length: [2.0, 1.0, 1.0]\nobstacles:\n{bx}{topt}")
w("tetra", "tet_sphere", f"mesh_type: tetra\nbox_length: [1.0, 1.0, 1.0]\nobstacles:\n{sph}{topt}")
w("tetra", "tet_several", f"mesh_type: tetra\nbox_length: [2.0, 1.0, 1.0]\nobstacles:\n{sph}{bx}{topt}")
for p in [(False,False,False),(True,False,False),(False,True,False),(False,False,True),(True,True,True),(False,True,True)]:
    w("tetra", "tet_per_" + "".join("t" if x else "f" for x in p), f"mesh_type: tetra\nbox_length: [2.0, 1.0, 1.0]\nperiodic: {yb(p)}\nobstacles:\n{sph}{bx}{topt}")
# hierarchy: nested levels of every mesh type. The writer makes one file pair per level.
hh = "hierarchy:\n  levels: 3\n"
w("hierarchy", "hier_gs_alt_per", f"mesh_type: gridsplit\nbox_length: [1.0, 1.0]\nn_divisions: [3, 2]\nperiodic: [true, false]\ndiagonal: alternate\n{hh}")
w("hierarchy", "hier_gs_stretch", "mesh_type: gridsplit\nbox_length: [1.0, 1.0]\nn_divisions: [4, 3]\ndiagonal: slash\nstretching_ratio: [1.3, 1.5]\nstretching_symmetric: [false, true]\nhierarchy:\n  levels: 4\n")
w("hierarchy", "hier_eq_per", f"mesh_type: equilateral\nbox_length: [1.0, 0.866]\nn_divisions: [4, 4]\nperiodic: [false, true]\n{hh}")
w("hierarchy", "hier_tri_several", f"mesh_type: triangle\nbox_length: [3.0, 1.0]\nperiodic: [true, true]\nobstacles:\n{circ}{rect}mesh_options:\n  mesh_size: 0.4\n  mesh_size_obstacle: 0.2\n{hh}")
w("hierarchy", "hier_kuhn_per", f"mesh_type: kuhn\nbox_length: [1.0, 1.0, 1.0]\nn_divisions: [2, 2, 2]\nperiodic: [true, false, true]\ndiagonal: main\n{hh}")
w("hierarchy", "hier_kuhn_alt_stretch", f"mesh_type: kuhn\nbox_length: [1.0, 1.0, 1.0]\nn_divisions: [2, 3, 2]\ndiagonal: alternate_layers\nstretching_ratio: [1.3, 1.0, 1.5]\nstretching_symmetric: [false, false, true]\n{hh}")
w("hierarchy", "hier_tet_none", f"mesh_type: tetra\nbox_length: [1.0, 1.0, 1.0]\nmesh_options:\n  mesh_size: 0.5\n{hh}")
w("hierarchy", "hier_tet_several", f"mesh_type: tetra\nbox_length: [2.0, 1.0, 1.0]\nperiodic: [false, true, true]\nobstacles:\n{sph}{bx}mesh_options:\n  mesh_size: 0.5\n  mesh_size_obstacle: 0.3\n{hh}")
