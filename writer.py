#!/usr/bin/env python3
"""
writer.py
==============
Generate and store 2-D / 3-D unstructured box meshes, optionally with
obstacles carved out as holes.

Five mesh types are supported:

  triangle    Delaunay-refined triangulation        (MeshPy / Triangle)
  tetra       Delaunay-refined tetrahedralization   (MeshPy / TetGen)
  gridsplit   Cartesian grid split into right triangles  (no external lib)
  equilateral Offset-row triangular lattice              (no external lib)
  kuhn        Kuhn / Freudenthal 6-tet-per-cube          (no external lib)

Usage
-----
    python writer.py config.yaml      # process a single config
    python writer.py configs_dir/     # process every '*.yaml' directly
                                       # inside configs_dir (non-recursive);
                                       # each gets its own .h5/.xdmf pair,
                                       # named and located per 'Output
                                       # files' below. A failure in one
                                       # config is reported and does not
                                       # stop the remaining ones from being
                                       # processed.

YAML keys
---------
  mesh_type   : triangle | tetra | gridsplit | equilateral | kuhn
  box_length  : [Lx, Ly]  or  [Lx, Ly, Lz]
  n_divisions : [nx, ny]  or  [nx, ny, nz]
      MeshPy types   — boundary vertices per edge; interior filled automatically.
      Structured types — number of cells per direction.
  periodic    : [bool, bool]  or  [bool, bool, bool]
  quality     : (MeshPy types only)
      min_angle  : float  (Triangle, degrees, default 20)
      max_area   : float  (Triangle; auto if omitted)
      max_volume : float  (TetGen;   auto if omitted)
  diagonal    : (gridsplit / kuhn only, see notes below)
  stretching_ratio     : (gridsplit / kuhn only, see 'Grid stretching' below)
  stretching_symmetric : (gridsplit / kuhn only, see 'Grid stretching' below)
  obstacles   : (triangle / tetra only — see 'Obstacles' below)

Output files
------------
  The .h5 and .xdmf files are written next to the config
  file, sharing its base name. E.g. running
      python writer.py meshes/triangle.yaml
  writes
      meshes/triangle.h5
      meshes/triangle.xdmf

  Provenance
  ----------
  After the .xdmf file above is written, three more blocks are appended to
  it as XML comments, after the closing </Xdmf> tag:
    1) the YAML config used for this run, verbatim
    2) the terminal output produced while building this mesh (every
       print() from load_config() through write_xdmf() -- i.e. everything
       a plain 'python writer.py config.yaml' would show)
    3) the full source of writer.py itself, as it existed at run time
  This makes every .xdmf file a self-contained record of exactly how it
  was produced. Literal '--' is not legal inside an XML comment, so each
  block is passed through _xml_comment_safe() first.

Notes on structured types
--------------------------
  gridsplit
      All triangles are congruent right triangles with legs dx and dy.
      Setting dx == dy gives isosceles right triangles.
      Optional 'diagonal' key selects which diagonal splits each cell:
          backslash (default) : cut along (i+1,j)-(i,j+1), i.e. "\"
          slash               : cut along (i,j)-(i+1,j+1),   i.e. "/"
          alternate           : checkerboard mix of the two (still all
                                congruent triangles, up to reflection)

      Optional 'stretching_ratio'/'stretching_symmetric' (per-axis lists,
      default [1.0, 1.0] / [false, false]) apply geometric grid stretching
      along x and/or y independently — see 'Grid stretching' below. Note:
      enabling stretching on either axis makes dx and/or dy vary across the
      grid, so cells are NO LONGER congruent; the structured tensor-product
      topology (and hence periodicity/diagonal handling) is unaffected.

  equilateral
      Even rows: nx+1 points at x = 0, dx, ..., Lx
      Odd  rows: nx+2 points at x = 0, dx/2, ..., Lx  (boundary points added)
      Interior triangles are equilateral when dy = dx * sqrt(3)/2, i.e.
          Ly / ny  =  (Lx / nx) * sqrt(3)/2  ≈  0.866 * Lx/nx
      Left/right boundary triangles are right-angle (unavoidable on a box).
      For y-periodic meshes choose ny even so top/bottom rows are both even.
      No stretching support (the offset-row layout does not generalize to
      a simple per-axis stretched grid).

  kuhn
      Each axis-aligned cube → 6 tetrahedra via all permutations of the three
      coordinate-step directions, walking from one cube corner to the
      opposite corner (the cube's internal "body diagonal").
      With dx=dy=dz=h all 6 tets per cube are congruent Sommerville type-1
      tets (edges h,h,h,h√2,h√2,h√3). For a non-cubic box (dx≠dy≠dz) the 6
      permutations within a single cube are NOT all congruent to each other
      — this is inherent to the Kuhn decomposition itself, independent of
      the 'diagonal' option below.

      Optional 'diagonal' key selects which body diagonal each cube splits
      along. A cube has 4 distinct body diagonals; all are uniformly
      conforming choices when applied to every cube:
          main             (default) : (0,0,0)-(1,1,1) — 'fixed' is a
                                        backward-compatible alias for this
          flip_x                     : (1,0,0)-(0,1,1) — main, mirrored along x
          flip_y                     : (0,1,0)-(1,0,1) — main, mirrored along y
          flip_z                     : (0,0,1)-(1,1,0) — main, mirrored along z
          alternate_layers            : the diagonal's starting corner varies
                                        independently per layer index,
                                        (i%2, j%2, k%2) — the family (up to
                                        relabeling) that removes a single
                                        global diagonal direction while
                                        remaining conforming across every
                                        shared cube face.
      For cubic cells (dx=dy=dz) the four fixed variants are related by the
      cube's rotational symmetry (same local tet shapes, different global
      grain direction); for non-cubic cells they are genuinely different
      tet-shape families.

      Optional 'stretching_ratio'/'stretching_symmetric' (per-axis lists,
      default [1.0, 1.0, 1.0] / [false, false, false]) apply geometric grid
      stretching along x, y and/or z independently — see 'Grid stretching'
      below. As with gridsplit, enabling stretching breaks the 'congruent
      Sommerville tet' guarantee described above; the 6-tets-per-cube
      decomposition and its conformity/periodicity properties are otherwise
      unaffected.

Grid stretching (gridsplit / kuhn only)
------------------------------------------
  stretching_ratio     : per-axis list, e.g. [1.0, 1.0] (2-D) or
                          [1.0, 1.0, 1.0] (3-D). 1.0 = uniform spacing
                          (default). Any other positive value applies a
                          geometric progression to the cell widths along
                          that axis.
  stretching_symmetric : per-axis list of bool, default all false.

  Non-symmetric (symmetric: false), ratio r != 1:
      Cell widths w_i = w0 * r**i for i = 0..n-1, chosen so they sum to the
      axis length. r > 1 clusters cells near coord=0 (small cells growing
      towards coord=L); r < 1 clusters cells near coord=L.

  Symmetric (symmetric: true), ratio r != 1:
      The axis is split at its midpoint; each half independently follows
      the same one-sided progression, mirrored about the midpoint.
      r > 1 clusters cells at BOTH boundaries (sparse at the center) — the
      classic wall-clustered grid used e.g. for boundary-layer resolution
      in channel/duct flows. r < 1 clusters cells at the CENTER (sparse at
      both boundaries).
      If n_divisions along that axis is odd, the two mirrored halves have
      slightly different cell counts, producing a small kink at the exact
      center; use an even n_divisions for a perfectly symmetric grid (a
      [warn] is printed when this occurs).

  A stretching ratio has no visible effect on an axis with only 1 cell
  (one-sided) or fewer than 2 cells per half (symmetric with n_divisions
  < 4), since a lone cell always spans its full span regardless of ratio.

  Stretching only redistributes vertex positions along each axis; it does
  not change the topology (cell count, diagonal pattern, periodicity
  matching) established by n_divisions/diagonal — only the physical
  spacing dx/dy/dz becomes non-uniform, so cells are no longer congruent
  once stretching is enabled on an axis.

Obstacles (triangle / tetra only)
----------------------------------
  Only the Delaunay-refined types support obstacles: they are implemented as
  PSLG/PLC holes (a closed boundary loop/surface + one interior point telling
  Triangle/TetGen "do not mesh here"). The structured types (gridsplit,
  equilateral, kuhn) tile the box uniformly with no facility for excluding
  regions, so 'obstacles' is rejected for those mesh_types.

  obstacles: a list of dicts, each with:
    2-D (mesh_type: triangle):
      shape: circle
        center: [cx, cy]
        radius: r
        n_segments: int (default 32)     # polygon approximation resolution
      shape: rectangle
        center: [cx, cy]
        half_extents: [hx, hy]
    3-D (mesh_type: tetra):
      shape: sphere
        center: [cx, cy, cz]
        radius: r
        subdivisions: int (default 2)    # icosphere refinement level
      shape: box
        center: [cx, cy, cz]
        half_extents: [hx, hy, hz]

  Each obstacle gets its own boundary marker: 10+idx for 2-D obstacles,
  100+idx for 3-D obstacles (idx = 0-based position in the 'obstacles' list),
  distinct from the box's own 1-4 / 1-6 face markers. Unlike the box (whose
  6 faces get individually distinguishable markers), a single obstacle's
  entire surface shares ONE marker — it is meant to tag "wall of obstacle
  #idx" as a whole, not distinguish parts of it.

  Sphere obstacles are polyhedral approximations (icosphere), not exact
  spheres; increase 'subdivisions' for a smoother surface at the cost of
  more surface triangles. There is no local mesh-size refinement near
  obstacles (unlike Gmsh's size fields) — refinement is governed globally
  by 'quality.max_volume' / 'max_area'.

  Placement is sanity-checked with simple bounding-region heuristics
  (printed as [warn], not hard failures): an obstacle reaching the domain
  boundary, or two obstacles whose bounding regions appear to overlap.
  Final geometric validity is enforced by Triangle/TetGen themselves.

HDF5 layout
-----------
  Attributes: mesh_type, dimension, box_length, n_vertices, n_cells, n_edges

  /points              float64  (n_pts,  3)     z=0 for 2-D
  /cell_to_vertex      int64    (n_cells, d)    d=3 (tri) | 4 (tet)
  /cell_to_edge        int64    (n_cells, e)    e=3 (tri) | 6 (tet); each
                                                entry indexes a row of
                                                /vertex_to_vertex. Local edge
                                                order follows
                                                itertools.combinations of the
                                                cell's local vertex indices:
                                                  tri  (0,1),(0,2),(1,2)
                                                  tet  (0,1),(0,2),(0,3),
                                                       (1,2),(1,3),(2,3)
  /cell_vol            float64  (n_cells,)      triangle area (2-D) |
                                                tetra volume (3-D)
  /vertex_to_vertex    int64    (n_edges, 2)    all unique sorted edges
  /vertex_to_cell/
      offsets          int64    (n_pts+1,)      CSR row pointers
      indices          int64    (Σ degree,)     CSR cell-index data
  /vertex_vol          float64  (n_pts,)        barycentric dual volume:
                                                sum, over cells incident to
                                                the vertex (per
                                                /vertex_to_cell), of
                                                cell_vol/(d+1) — i.e. each
                                                cell's area/volume split
                                                equally among its vertices
                                                (1/3 for triangles, 1/4 for
                                                tetrahedra). Periodicity IS
                                                enforced: vertices linked via
                                                /periodicity (x_pairs,
                                                y_pairs, z_pairs) — including
                                                transitively across more
                                                than one direction, e.g. a
                                                periodic-box corner — share
                                                ONE combined value, the sum
                                                of their individually
                                                computed contributions.
  /boundary/
      point_markers    int32    (n_pts,)        0=interior; see convention
      edge_markers     int32    (n_edges,)      2-D only
  /periodicity/
      Attrs: periodic_x, periodic_y [, periodic_z]   bool
      x_pairs          int64    (n, 2)          [v_lo, v_hi]
      y_pairs / z_pairs                         idem
  /obstacles/                                   present only if obstacles given
      Attrs: count
      obstacle_<idx>/  Attrs: shape, marker, center, radius | half_extents

Boundary marker convention
--------------------------
  2-D:  1=bottom(y=0)  2=right(x=Lx)  3=top(y=Ly)  4=left(x=0)
        10+idx = obstacle idx's surface
  3-D:  1=xmin  2=xmax  3=ymin  4=ymax  5=zmin  6=zmax   0=interior
        100+idx = obstacle idx's surface
  Corner/edge vertices receive the highest-priority face marker (lowest number).
"""

import io
import sys
from itertools import combinations, permutations as _iperms
from pathlib import Path

import h5py
import numpy as np
import yaml

# ── Constants ──────────────────────────────────────────────────────────────────
_MK2  = {"bottom": 1, "right": 2, "top": 3, "left": 4}
_MK3  = {"xmin": 1, "xmax": 2, "ymin": 3, "ymax": 4, "zmin": 5, "zmax": 6}
_DIRS = ["x", "y", "z"]

_2D_TYPES = {"triangle", "gridsplit", "equilateral"}
_3D_TYPES = {"tetra", "kuhn"}
_ALL_TYPES = _2D_TYPES | _3D_TYPES
_OBSTACLE_TYPES = {"triangle", "tetra"}
_OBSTACLE_SHAPES_2D = {"circle", "rectangle"}
_OBSTACLE_SHAPES_3D = {"sphere", "box"}
_STRETCH_TYPES = {"gridsplit", "kuhn"}

# Pre-computed Kuhn permutations (used in _build_kuhn)
_KUHN_PERMS = list(_iperms([0, 1, 2]))   # 6 permutations of {x, y, z}


# ── stdout capture (provenance) ─────────────────────────────────────────────────

class _Tee:
    """Mirror writes to multiple streams at once (used to capture the run's
    terminal output for embedding in the .xdmf file, while still printing
    normally to the real terminal)."""
    def __init__(self, *streams):
        self.streams = streams
    def write(self, data):
        for s in self.streams:
            s.write(data)
    def flush(self):
        for s in self.streams:
            s.flush()


# ── Configuration ──────────────────────────────────────────────────────────────

def load_config(path: str) -> dict:
    with open(path) as fh:
        cfg = yaml.safe_load(fh)
    mt = cfg.get("mesh_type", "triangle")
    if mt not in _ALL_TYPES:
        raise ValueError(f"mesh_type must be one of {sorted(_ALL_TYPES)}, got '{mt}'")
    dim = 2 if mt in _2D_TYPES else 3
    cfg.setdefault("periodic",  [False] * dim)
    cfg.setdefault("quality",   {})
    cfg.setdefault("obstacles", [])
    cfg.setdefault("stretching_ratio",     [1.0] * dim)
    cfg.setdefault("stretching_symmetric", [False] * dim)

    for key in ("box_length", "n_divisions", "periodic",
               "stretching_ratio", "stretching_symmetric"):
        if len(cfg[key]) != dim:
            raise ValueError(f"'{key}' needs {dim} entries for mesh_type='{mt}'")
    return cfg


# ── Unique-point registry (MeshPy builders only) ───────────────────────────────

class _Reg:
    _PREC = 10
    def __init__(self): self._pts, self._idx = [], {}
    def add(self, *c) -> int:
        key = tuple(round(v, self._PREC) for v in c)
        if key not in self._idx:
            self._idx[key] = len(self._pts); self._pts.append(c)
        return self._idx[key]
    def array(self) -> np.ndarray:
        return np.asarray(self._pts, dtype=np.float64)


# ── Obstacle validation ────────────────────────────────────────────────────────

def _validate_obstacles(obstacles: list, box: list, dim: int, mesh_type: str) -> None:
    """Raise on malformed obstacle specs; print [warn] for placement heuristics."""
    if not obstacles:
        return
    if mesh_type not in _OBSTACLE_TYPES:
        raise ValueError(
            f"'obstacles' is only supported for mesh_type in {sorted(_OBSTACLE_TYPES)} "
            f"(structured types tile the box uniformly with no facility for "
            f"excluding regions); got mesh_type='{mesh_type}'")

    valid_shapes = _OBSTACLE_SHAPES_2D if dim == 2 else _OBSTACLE_SHAPES_3D
    reach = []   # (center array, scalar reach) per obstacle, for overlap heuristic

    for i, obs in enumerate(obstacles):
        shape = obs.get("shape")
        if shape not in valid_shapes:
            raise ValueError(f"obstacle {i}: shape must be one of {sorted(valid_shapes)} "
                             f"for a {dim}-D mesh, got '{shape}'")
        c = obs.get("center")
        if c is None or len(c) != dim:
            raise ValueError(f"obstacle {i}: 'center' must have {dim} entries")
        if shape in ("circle", "sphere"):
            r = obs.get("radius")
            if r is None or r <= 0:
                raise ValueError(f"obstacle {i}: 'radius' must be a positive number")
            half = [r] * dim
        else:
            he = obs.get("half_extents")
            if he is None or len(he) != dim:
                raise ValueError(f"obstacle {i}: 'half_extents' must have {dim} entries")
            half = list(he)

        for d in range(dim):
            if c[d] - half[d] <= 0 or c[d] + half[d] >= box[d]:
                print(f"  [warn] obstacle {i}: extends to/past the domain boundary "
                      f"along axis {_DIRS[d]}")
        reach.append((np.array(c, dtype=float), max(half)))

    for i in range(len(reach)):
        for j in range(i + 1, len(reach)):
            ci, hi = reach[i]; cj, hj = reach[j]
            if np.linalg.norm(ci - cj) < (hi + hj):
                print(f"  [warn] obstacles {i} and {j}: bounding regions appear to "
                      f"overlap — check placement")


# ── Obstacle geometry helpers ──────────────────────────────────────────────────

def _circle_points(cx, cy, r, n):
    """n points evenly spaced around a circle of radius r centered at (cx,cy)."""
    ang = np.linspace(0., 2*np.pi, n, endpoint=False)
    return [(cx + r*np.cos(a), cy + r*np.sin(a)) for a in ang]


def _rectangle_points(cx, cy, hx, hy):
    """4 corners of an axis-aligned rectangle, CCW, centered at (cx,cy)."""
    return [(cx-hx, cy-hy), (cx+hx, cy-hy), (cx+hx, cy+hy), (cx-hx, cy+hy)]


def _icosphere(cx, cy, cz, r, subdivisions=2):
    """
    Closed, watertight triangulated sphere via recursive icosahedron
    subdivision (outward-pointing normals, i.e. away from the center).

    Returns (vertices (n,3), faces (m,3) local indices, surf_tol).
    surf_tol is the max deviation of a flat facet from the true sphere —
    needed because Steiner points TetGen inserts on this polyhedral
    surface satisfy dist(center) ∈ [r - surf_tol, r], not exactly r.
    """
    t = (1.0 + np.sqrt(5.0)) / 2.0
    base = [
        (-1,  t,  0), ( 1,  t,  0), (-1, -t,  0), ( 1, -t,  0),
        ( 0, -1,  t), ( 0,  1,  t), ( 0, -1, -t), ( 0,  1, -t),
        ( t,  0, -1), ( t,  0,  1), (-t,  0, -1), (-t,  0,  1),
    ]
    vlist = [np.array(v, dtype=np.float64) for v in base]
    for i in range(len(vlist)):
        vlist[i] = vlist[i] / np.linalg.norm(vlist[i])

    faces = [
        (0,11,5), (0,5,1), (0,1,7), (0,7,10), (0,10,11),
        (1,5,9), (5,11,4), (11,10,2), (10,7,6), (7,1,8),
        (3,9,4), (3,4,2), (3,2,6), (3,6,8), (3,8,9),
        (4,9,5), (2,4,11), (6,2,10), (8,6,7), (9,8,1),
    ]

    for _ in range(subdivisions):
        midcache = {}
        def midpoint(i1, i2):
            key = (min(i1, i2), max(i1, i2))
            if key in midcache:
                return midcache[key]
            m = (vlist[i1] + vlist[i2]) / 2.0
            m = m / np.linalg.norm(m)
            vlist.append(m)
            midcache[key] = len(vlist) - 1
            return midcache[key]

        new_faces = []
        for a, b, c in faces:
            ab, bc, ca = midpoint(a, b), midpoint(b, c), midpoint(c, a)
            new_faces += [(a, ab, ca), (b, bc, ab), (c, ca, bc), (ab, bc, ca)]
        faces = new_faces

    verts_unit = np.array(vlist, dtype=np.float64)
    centroids  = verts_unit[np.array(faces)].mean(axis=1)
    surf_tol_unit = 1.0 - float(np.min(np.linalg.norm(centroids, axis=1)))

    verts = verts_unit * r + np.array([cx, cy, cz])
    return verts, np.array(faces, dtype=np.int64), surf_tol_unit * r


def _box_obstacle_surface(reg, cx, cy, cz, hx, hy, hz, marker):
    """Closed axis-aligned box surface (12 triangles, outward normals)."""
    x0, x1 = cx-hx, cx+hx
    y0, y1 = cy-hy, cy+hy
    z0, z1 = cz-hz, cz+hz
    facets, fmk = [], []
    def quad(a, b, c, d):
        facets.append([a, b, c]); facets.append([a, c, d]); fmk.extend([marker, marker])
    quad(reg.add(x0,y0,z0), reg.add(x0,y0,z1), reg.add(x0,y1,z1), reg.add(x0,y1,z0))  # xmin
    quad(reg.add(x1,y0,z0), reg.add(x1,y1,z0), reg.add(x1,y1,z1), reg.add(x1,y0,z1))  # xmax
    quad(reg.add(x0,y0,z0), reg.add(x1,y0,z0), reg.add(x1,y0,z1), reg.add(x0,y0,z1))  # ymin
    quad(reg.add(x0,y1,z0), reg.add(x0,y1,z1), reg.add(x1,y1,z1), reg.add(x1,y1,z0))  # ymax
    quad(reg.add(x0,y0,z0), reg.add(x0,y1,z0), reg.add(x1,y1,z0), reg.add(x1,y0,z0))  # zmin
    quad(reg.add(x0,y0,z1), reg.add(x1,y0,z1), reg.add(x1,y1,z1), reg.add(x0,y1,z1))  # zmax
    return facets, fmk


# ── MeshPy-based builders ──────────────────────────────────────────────────────

def _build_2d(cfg: dict):
    """Delaunay-refined triangulation via Triangle / MeshPy, optionally with
    circular/rectangular obstacles carved out as holes (see 'obstacles' key)."""
    import meshpy.triangle as triangle

    Lx, Ly = cfg["box_length"]; nx, ny = cfg["n_divisions"]; q = cfg["quality"]
    obstacles = cfg.get("obstacles", [])
    dx, dy  = Lx / max(nx-1, 1), Ly / max(ny-1, 1)
    max_area = q.get("max_area") or (0.5 * dx * dy)
    min_ang  = q.get("min_angle", 20.0)
    reg = _Reg()
    xv, yv = np.linspace(0., Lx, nx), np.linspace(0., Ly, ny)
    segs, smk = [], []

    def strip(seq, mk):
        for a, b in zip(seq[:-1], seq[1:]): segs.append((a, b)); smk.append(mk)

    def ring(seq, mk):
        for a, b in zip(seq, seq[1:] + seq[:1]): segs.append((a, b)); smk.append(mk)

    strip([reg.add(x, 0.) for x in xv],        _MK2["bottom"])
    strip([reg.add(Lx, y) for y in yv],         _MK2["right"])
    strip([reg.add(x, Ly) for x in xv[::-1]],  _MK2["top"])
    strip([reg.add(0., y) for y in yv[::-1]],  _MK2["left"])

    holes = []
    obstacle_meta = []
    for idx, obs in enumerate(obstacles):
        marker = 10 + idx
        shape  = obs["shape"]
        cx, cy = obs["center"]
        if shape == "circle":
            r = obs["radius"]; n = obs.get("n_segments", 32)
            pts_ring = [reg.add(x, y) for x, y in _circle_points(cx, cy, r, n)]
        else:   # rectangle
            hx, hy = obs["half_extents"]
            pts_ring = [reg.add(x, y) for x, y in _rectangle_points(cx, cy, hx, hy)]
        ring(pts_ring, marker)
        holes.append((cx, cy))
        obstacle_meta.append({**obs, "marker": marker})

    info = triangle.MeshInfo()
    info.set_points(reg.array().tolist())
    info.set_facets(segs, facet_markers=smk)
    if holes:
        info.set_holes(holes)

    # MeshPy ≥ 2024 renamed generate_edges → generate_faces
    # and mesh.edges / mesh.edge_markers → mesh.faces / mesh.face_markers.
    # The try/except keeps the script working on both old and new installs.
    try:
        mesh    = triangle.build(info, max_volume=max_area, min_angle=min_ang,
                                 generate_faces=True)
        edges   = np.asarray(mesh.faces,        dtype=np.int64)
        edge_mk = np.asarray(mesh.face_markers, dtype=np.int32)
    except TypeError:           # older MeshPy: generate_faces not recognised
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mesh = triangle.build(info, max_volume=max_area, min_angle=min_ang,
                                  generate_edges=True)
        edges   = np.asarray(mesh.edges,        dtype=np.int64)
        edge_mk = np.asarray(mesh.edge_markers, dtype=np.int32)

    return (np.asarray(mesh.points,        dtype=np.float64),
            np.asarray(mesh.elements,      dtype=np.int64),
            edges, np.asarray(mesh.point_markers, dtype=np.int32), edge_mk,
            obstacle_meta)


def _build_3d(cfg: dict):
    """Delaunay-refined tetrahedralization via TetGen / MeshPy, optionally with
    spherical/box obstacles carved out as holes (see 'obstacles' key)."""
    import meshpy.tet as tet

    Lx, Ly, Lz = cfg["box_length"]; nx, ny, nz = cfg["n_divisions"]; q = cfg["quality"]
    obstacles = cfg.get("obstacles", [])
    dx = Lx/max(nx-1,1); dy = Ly/max(ny-1,1); dz = Lz/max(nz-1,1)
    max_vol = q.get("max_volume") or (dx*dy*dz/6.)
    reg = _Reg()
    xv = np.linspace(0., Lx, nx); yv = np.linspace(0., Ly, ny); zv = np.linspace(0., Lz, nz)
    facets, fmk = [], []
    def quad(a, b, c, d, mk):
        # New pybind11 MeshPy: each facet is a flat vertex list, not [[v0,v1,v2]]
        facets.append([a,b,c]); facets.append([a,c,d]); fmk.extend([mk, mk])
    for j in range(ny-1):
        for k in range(nz-1):
            y0,y1,z0,z1 = yv[j],yv[j+1],zv[k],zv[k+1]
            quad(reg.add(0.,y0,z0),reg.add(0.,y0,z1),reg.add(0.,y1,z1),reg.add(0.,y1,z0),_MK3["xmin"])
            quad(reg.add(Lx,y0,z0),reg.add(Lx,y1,z0),reg.add(Lx,y1,z1),reg.add(Lx,y0,z1),_MK3["xmax"])
    for i in range(nx-1):
        for k in range(nz-1):
            x0,x1,z0,z1 = xv[i],xv[i+1],zv[k],zv[k+1]
            quad(reg.add(x0,0.,z0),reg.add(x1,0.,z0),reg.add(x1,0.,z1),reg.add(x0,0.,z1),_MK3["ymin"])
            quad(reg.add(x0,Ly,z0),reg.add(x0,Ly,z1),reg.add(x1,Ly,z1),reg.add(x1,Ly,z0),_MK3["ymax"])
    for i in range(nx-1):
        for j in range(ny-1):
            x0,x1,y0,y1 = xv[i],xv[i+1],yv[j],yv[j+1]
            quad(reg.add(x0,y0,0.),reg.add(x0,y1,0.),reg.add(x1,y1,0.),reg.add(x1,y0,0.),_MK3["zmin"])
            quad(reg.add(x0,y0,Lz),reg.add(x1,y0,Lz),reg.add(x1,y1,Lz),reg.add(x0,y1,Lz),_MK3["zmax"])

    # ── Obstacles: closed surfaces carved out as TetGen holes ──────────────
    holes = []
    obstacle_meta = []      # kept for the point-marker fallback below
    for idx, obs in enumerate(obstacles):
        marker = 100 + idx
        shape  = obs["shape"]
        cx, cy, cz = obs["center"]
        if shape == "sphere":
            r = obs["radius"]; subdiv = obs.get("subdivisions", 2)
            verts, faces, sagitta = _icosphere(cx, cy, cz, r, subdiv)
            idx_map = [reg.add(*v) for v in verts]
            for a, b, c in faces:
                facets.append([idx_map[a], idx_map[b], idx_map[c]]); fmk.append(marker)
            obstacle_meta.append({**obs, "marker": marker, "surf_tol": max(sagitta, 1e-9)})
        else:   # box
            hx, hy, hz = obs["half_extents"]
            ob_facets, ob_fmk = _box_obstacle_surface(reg, cx, cy, cz, hx, hy, hz, marker)
            facets += ob_facets; fmk += ob_fmk
            obstacle_meta.append({**obs, "marker": marker})
        holes.append((cx, cy, cz))

    info = tet.MeshInfo()
    info.set_points(reg.array().tolist()); info.set_facets(facets, markers=fmk)
    if holes:
        info.set_holes(holes)

    mesh = tet.build(info, options=tet.Options("pqe"), max_volume=max_vol)
    pts   = np.asarray(mesh.points,        dtype=np.float64)
    cells = np.asarray(mesh.elements,      dtype=np.int64)

    # Try TetGen's native point_markers first (cheap: no extra pass over pts).
    # Known failure mode: with no obstacles it can come back entirely
    # unallocated. Defensively also check the *values* are all within the
    # set of markers we actually assigned to input facets (0 = interior,
    # 1-6 = box faces, obstacle markers) — if TetGen ever produced stray
    # values outside that set we fall back rather than trust them silently.
    valid_markers = {0, 1, 2, 3, 4, 5, 6} | {obs["marker"] for obs in obstacle_meta}
    try:
        pt_mk = np.asarray(mesh.point_markers, dtype=np.int32)
        if pt_mk.size == 0 or not np.all(np.isin(pt_mk, list(valid_markers))):
            raise RuntimeError("unallocated or contains values outside the "
                               "expected marker set")
    except RuntimeError as exc:
        print(f"  [info] TetGen point_markers {exc}; deriving from coordinates")
        pt_mk = _pt_markers_3d(pts, Lx, Ly, Lz, obstacle_meta)

    try:
        edges = np.asarray(mesh.edges, dtype=np.int64)
        if edges.ndim != 2 or edges.shape[1] != 2 or len(edges) == 0: raise ValueError
    except (AttributeError, ValueError) as exc:
        print(f"  [info] TetGen edge output unavailable ({exc}); deriving from cells")
        edges = _edges_from_cells(cells)

    return pts, cells, edges, pt_mk, obstacle_meta


# ── Grid stretching (gridsplit / kuhn only) ────────────────────────────────────

def _stretched_coords(L: float, n: int, ratio: float = 1.0,
                      symmetric: bool = False) -> np.ndarray:
    """
    Compute n+1 vertex coordinates along [0, L] with optional geometric grid
    stretching. See the module docstring's 'Grid stretching' section for the
    full semantics of 'ratio' and 'symmetric'.

    ratio == 1.0 -> uniform spacing, regardless of 'symmetric'.
    symmetric and n < 2 -> falls back to the one-sided formula (a lone cell,
    or a domain that cannot be split into at least 1 cell per half, always
    just spans its own length; there is nothing to mirror).
    """
    if abs(ratio - 1.0) < 1e-12:
        return np.linspace(0., L, n + 1)

    def _one_sided(length, n_cells):
        w0 = length * (ratio - 1.0) / (ratio**n_cells - 1.0)
        return w0 * ratio**np.arange(n_cells)

    if not symmetric or n < 2:
        widths = _one_sided(L, n)
    else:
        n1, n2 = n // 2, n - n // 2
        half = L / 2.0
        w_left  = _one_sided(half, n1) if n1 else np.array([])
        w_right = (_one_sided(half, n2) if n2 else np.array([]))[::-1]
        widths  = np.concatenate([w_left, w_right])

    coords = np.concatenate([[0.], np.cumsum(widths)])
    coords[-1] = L   # guard against floating-point drift
    return coords


def _validate_stretching(ratio: list, symmetric: list, n_divisions: list,
                         dim: int, mesh_type: str) -> None:
    """Raise on malformed stretching specs; print [warn] for odd-n symmetric axes."""
    is_default = all(abs(r - 1.0) < 1e-12 for r in ratio) and not any(symmetric)
    if mesh_type not in _STRETCH_TYPES and not is_default:
        raise ValueError(
            f"'stretching_ratio'/'stretching_symmetric' are only supported for "
            f"mesh_type in {sorted(_STRETCH_TYPES)} (grid-based structured "
            f"types); got mesh_type='{mesh_type}'")
    for d, r in enumerate(ratio):
        if r <= 0:
            raise ValueError(f"stretching_ratio[{d}] must be > 0, got {r}")
    if mesh_type in _STRETCH_TYPES:
        for d, (sym, n) in enumerate(zip(symmetric, n_divisions)):
            if sym and n % 2:
                print(f"  [warn] stretching_symmetric on axis {_DIRS[d]}: "
                      f"n_divisions={n} is odd -- the two mirrored halves "
                      f"will have slightly different cell counts, giving a "
                      f"small kink at the center; use an even value for a "
                      f"perfectly symmetric grid")


# ── Structured builders ────────────────────────────────────────────────────────

def _build_gridsplit(cfg: dict):
    """
    Cartesian (nx+1)×(ny+1) grid, each rectangle split into two triangles.
    Produces 2·nx·ny congruent right triangles with legs dx=Lx/nx, dy=Ly/ny
    (congruent up to reflection when diagonal='alternate'), UNLESS grid
    stretching is enabled on an axis, in which case dx and/or dy vary and
    cells are no longer congruent (see 'stretching_ratio'/'stretching_symmetric'
    below and the module docstring's 'Grid stretching' section).

    diagonal:
        'backslash' (default) — split along (i+1,j)-(i,j+1)   "\\"
        'slash'                — split along (i,j)-(i+1,j+1)   "/"
        'alternate'            — checkerboard mix of the two, i.e. the
                                  diagonal direction flips on alternating
                                  cells (i+j odd vs even). This removes the
                                  directional bias a single fixed diagonal
                                  introduces (e.g. for isotropic-looking
                                  gradients), while every triangle is still
                                  congruent to every other up to reflection.

    stretching_ratio / stretching_symmetric:
        Per-axis [rx, ry] / [sym_x, sym_y], default [1.0, 1.0] / [false, false].
        Replaces the uniform np.linspace grid along each axis with a
        geometrically stretched one; see _stretched_coords.
    """
    Lx, Ly = cfg["box_length"]; nx, ny = cfg["n_divisions"]
    diagonal = cfg.get("diagonal", "backslash")
    valid = {"backslash", "slash", "alternate"}
    if diagonal not in valid:
        raise ValueError(f"'diagonal' must be one of {sorted(valid)}, got '{diagonal}'")

    rx, ry = cfg.get("stretching_ratio", [1.0, 1.0])
    sx, sy = cfg.get("stretching_symmetric", [False, False])
    xs = _stretched_coords(Lx, nx, rx, sx)
    ys = _stretched_coords(Ly, ny, ry, sy)
    XX, YY = np.meshgrid(xs, ys)
    pts = np.column_stack([XX.ravel(), YY.ravel()])   # row-major: j*(nx+1)+i

    def vid(i, j): return j*(nx+1)+i

    cells = []
    for j in range(ny):
        for i in range(nx):
            bl, br, tr, tl = vid(i,j), vid(i+1,j), vid(i+1,j+1), vid(i,j+1)
            use_backslash = (diagonal == "backslash" or
                             (diagonal == "alternate" and (i + j) % 2 == 0))
            if use_backslash:               # diagonal br-tl  "\"
                cells.append([bl, br, tl])
                cells.append([br, tr, tl])
            else:                            # diagonal bl-tr  "/"
                cells.append([bl, br, tr])
                cells.append([bl, tr, tl])

    return pts, np.array(cells, dtype=np.int64)


def _build_equilateral(cfg: dict):
    """
    Offset-row triangular lattice on [0,Lx]×[0,Ly].

    Row layout
    ----------
      Even rows (j=0,2,…):  nx+1 points at x = 0, dx, …, Lx
      Odd  rows (j=1,3,…):  nx+2 points at x = 0, dx/2, …, Lx-dx/2, Lx
        (x=0 and x=Lx are added as boundary anchors)

    Triangle count per strip:  2·nx+1  (nx interior + nx-1 interior + 2 boundary)
    Total triangles:  ny·(2·nx+1)

    Equilateral condition
    ---------------------
      Interior triangles are equilateral when dy = dx·√3/2.
      Left/right boundary triangles are right-angle (unavoidable on a box).
      For y-periodic meshes choose ny even (both bounding rows are even-type).
    """
    Lx, Ly = cfg["box_length"]; nx, ny = cfg["n_divisions"]
    dx, dy = Lx/nx, Ly/ny
    dy_eq  = dx*np.sqrt(3)/2
    if abs(dy - dy_eq)/dy_eq > 0.02:
        print(f"  [info] equilateral: interior tris equilateral at "
              f"dy={dy_eq:.5f} (Ly/ny={dy:.5f}, {100*(dy/dy_eq-1):+.1f}%)\n"
              f"         Suggestion: set Ly = {nx*dy_eq*ny/nx:.5f} or ny = "
              f"{max(1,round(Ly/dy_eq))}")
    if ny % 2 and cfg.get("periodic", [False,False])[1]:
        print("  [warn] equilateral: ny is odd; top row is odd-type — "
              "y-periodic pairing may be incomplete")

    pts = []; row_start = []
    for j in range(ny+1):
        row_start.append(len(pts)); y = j*dy
        if j % 2 == 0:                          # even: nx+1 pts
            for i in range(nx+1): pts.append([i*dx, y])
        else:                                   # odd:  nx+2 pts
            pts.append([0., y])
            for i in range(nx): pts.append([(i+.5)*dx, y])
            pts.append([Lx, y])

    pts   = np.array(pts, dtype=np.float64)
    cells = []
    for j in range(ny):
        s0, s1 = row_start[j], row_start[j+1]
        if j % 2 == 0:
            # even (nx+1) → odd (nx+2)
            # E: s0+0 … s0+nx   |   O: s1+0(x=0), s1+1…s1+nx(int), s1+nx+1(x=Lx)
            cells.append([s0,    s1,    s1+1])              # left  boundary
            for i in range(nx):
                cells.append([s0+i, s0+i+1, s1+i+1])       # upward   ▲
                if i < nx-1:
                    cells.append([s0+i+1, s1+i+1, s1+i+2]) # downward ▽
            cells.append([s0+nx, s1+nx, s1+nx+1])           # right boundary
        else:
            # odd (nx+2) → even (nx+1)
            # O: s0+0(x=0), s0+1…s0+nx(int), s0+nx+1(x=Lx)   |   E: s1+0…s1+nx
            cells.append([s0,      s1,    s0+1])             # left  boundary
            for i in range(nx):
                cells.append([s0+i+1, s1+i, s1+i+1])        # downward ▽
                if i < nx-1:
                    cells.append([s0+i+1, s0+i+2, s1+i+1])  # upward   ▲
            cells.append([s0+nx+1, s0+nx, s1+nx])            # right boundary

    return pts, np.array(cells, dtype=np.int64)


def _build_kuhn(cfg: dict):
    """
    Kuhn / Freudenthal decomposition: 6 tetrahedra per cube.

    For a chosen starting corner offset o=(ox,oy,oz) ∈ {0,1}³ within cube
    (i,j,k), and for each permutation σ=(σ₁,σ₂,σ₃) of {x,y,z}, the walk
    starts at corner o and flips one coordinate at a time (in axis order σ)
    until it reaches the opposite corner (1-ox, 1-oy, 1-oz):
        v₀ = (i+ox,     j+oy,     k+oz)
        v₁ = v₀ with coordinate σ₁ flipped
        v₂ = v₁ with coordinate σ₂ flipped
        v₃ = v₂ with coordinate σ₃ flipped  = (i+1-ox, j+1-oy, k+1-oz)
    All 6 permutations share the same body-diagonal edge v₀-v₃.

    diagonal:
        'main'             (default) — o=(0,0,0) for every cube.
                                        'fixed' is a backward-compatible alias.
        'flip_x'                     — o=(1,0,0) for every cube.
        'flip_y'                     — o=(0,1,0) for every cube.
        'flip_z'                     — o=(0,0,1) for every cube.
        'alternate_layers'            — o=(i%2, j%2, k%2): the diagonal
                                        flips independently per layer index.

    All five modes are conforming (adjacent cubes agree on shared-face
    triangulations). The four fixed variants ('main', 'flip_x', 'flip_y',
    'flip_z') apply the same constant offset to every cube, so a face
    shared by two cubes trivially carries the same diagonal on both sides.
    'alternate_layers' is the general conforming family derived from the
    weaker constraint that the offset component normal to a face direction
    may only depend on the layer index along that direction.

    stretching_ratio / stretching_symmetric:
        Per-axis [rx, ry, rz] / [sym_x, sym_y, sym_z], default
        [1.0, 1.0, 1.0] / [false, false, false]. Replaces the uniform
        np.linspace grid along each axis with a geometrically stretched
        one; see _stretched_coords. Breaks the 'congruent Sommerville tet'
        guarantee on any axis where it is enabled.
    """
    Lx, Ly, Lz = cfg["box_length"]; nx, ny, nz = cfg["n_divisions"]

    raw = cfg.get("diagonal", "main")
    diagonal = "main" if raw == "fixed" else raw   # 'fixed' = backward-compat alias

    fixed_offsets = {
        "main":   (0, 0, 0),
        "flip_x": (1, 0, 0),
        "flip_y": (0, 1, 0),
        "flip_z": (0, 0, 1),
    }
    if diagonal not in fixed_offsets and diagonal != "alternate_layers":
        valid_display = sorted(fixed_offsets) + ["alternate_layers", "fixed"]
        raise ValueError(f"'diagonal' must be one of {valid_display}, got '{raw}'")

    rx, ry, rz = cfg.get("stretching_ratio", [1.0, 1.0, 1.0])
    sx, sy, sz = cfg.get("stretching_symmetric", [False, False, False])
    xs = _stretched_coords(Lx, nx, rx, sx)
    ys = _stretched_coords(Ly, ny, ry, sy)
    zs = _stretched_coords(Lz, nz, rz, sz)

    def vid(i, j, k): return k*(ny+1)*(nx+1) + j*(nx+1) + i

    pts = np.array([[xs[i], ys[j], zs[k]]
                    for k in range(nz+1)
                    for j in range(ny+1)
                    for i in range(nx+1)], dtype=np.float64)

    if diagonal == "alternate_layers":
        def start_offset(i, j, k): return (i % 2, j % 2, k % 2)
    else:
        _o = fixed_offsets[diagonal]
        def start_offset(i, j, k): return _o

    cells = []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                o = list(start_offset(i, j, k))
                for perm in _KUHN_PERMS:
                    cur = list(o)
                    tet = [vid(i+cur[0], j+cur[1], k+cur[2])]
                    for axis in perm:
                        cur[axis] = 1 - cur[axis]
                        tet.append(vid(i+cur[0], j+cur[1], k+cur[2]))
                    cells.append(tet)

    return pts, np.array(cells, dtype=np.int64)


# ── Boundary markers for structured meshes ─────────────────────────────────────

def _pt_markers_2d(pts: np.ndarray, Lx: float, Ly: float,
                   tol: float = 1e-9) -> np.ndarray:
    """Assign boundary markers; lower number wins at corners (bottom=1 has priority)."""
    mk = np.zeros(len(pts), dtype=np.int32)
    x, y = pts[:, 0], pts[:, 1]
    mk[np.abs(x)        < tol] = _MK2["left"]
    mk[np.abs(y - Ly)   < tol] = _MK2["top"]
    mk[np.abs(x - Lx)   < tol] = _MK2["right"]
    mk[np.abs(y)        < tol] = _MK2["bottom"]   # highest priority
    return mk


def _edge_markers_2d(pts: np.ndarray, edges: np.ndarray,
                     Lx: float, Ly: float, tol: float = 1e-9) -> np.ndarray:
    """An edge is a boundary edge iff both endpoints lie on the same boundary face."""
    a, b = edges[:, 0], edges[:, 1]
    xa, ya = pts[a, 0], pts[a, 1]
    xb, yb = pts[b, 0], pts[b, 1]
    mk = np.zeros(len(edges), dtype=np.int32)
    mk[(np.abs(xa)      < tol) & (np.abs(xb)      < tol)] = _MK2["left"]
    mk[(np.abs(ya-Ly)   < tol) & (np.abs(yb-Ly)   < tol)] = _MK2["top"]
    mk[(np.abs(xa-Lx)   < tol) & (np.abs(xb-Lx)   < tol)] = _MK2["right"]
    mk[(np.abs(ya)      < tol) & (np.abs(yb)       < tol)] = _MK2["bottom"]
    return mk


def _pt_markers_3d(pts: np.ndarray, Lx: float, Ly: float, Lz: float,
                   obstacles=None, tol: float = 1e-9) -> np.ndarray:
    """
    Assign boundary markers from coordinates: box faces (1-6, lower number
    wins at corners/edges) plus, if given, obstacle surfaces (100+idx).

    For sphere obstacles the surface is a polyhedral (icosphere)
    approximation, so points on it satisfy dist(center) ∈ [r - surf_tol, r]
    rather than exactly r; 'surf_tol' (computed in _icosphere) accounts
    for that. Box obstacles are exact planes, checked directly.
    """
    mk = np.zeros(len(pts), dtype=np.int32)
    x, y, z = pts[:, 0], pts[:, 1], pts[:, 2]
    mk[np.abs(z - Lz) < tol] = _MK3["zmax"]
    mk[np.abs(z)      < tol] = _MK3["zmin"]
    mk[np.abs(y - Ly) < tol] = _MK3["ymax"]
    mk[np.abs(y)      < tol] = _MK3["ymin"]
    mk[np.abs(x - Lx) < tol] = _MK3["xmax"]
    mk[np.abs(x)      < tol] = _MK3["xmin"]   # highest priority among box faces

    for obs in (obstacles or []):
        marker = obs["marker"]
        interior = (mk == 0)   # never overwrite a box-face marker
        if obs["shape"] == "sphere":
            cx, cy, cz = obs["center"]; r = obs["radius"]
            band = max(obs.get("surf_tol", 1e-9), tol)
            dist = np.sqrt((x-cx)**2 + (y-cy)**2 + (z-cz)**2)
            on_surf = (dist <= r + tol) & (dist >= r - band)
        else:  # box
            cx, cy, cz = obs["center"]; hx, hy, hz = obs["half_extents"]
            x0,x1,y0,y1,z0,z1 = cx-hx,cx+hx,cy-hy,cy+hy,cz-hz,cz+hz
            near_x = (np.abs(x-x0)<tol) | (np.abs(x-x1)<tol)
            near_y = (np.abs(y-y0)<tol) | (np.abs(y-y1)<tol)
            near_z = (np.abs(z-z0)<tol) | (np.abs(z-z1)<tol)
            in_x = (x>=x0-tol) & (x<=x1+tol)
            in_y = (y>=y0-tol) & (y<=y1+tol)
            in_z = (z>=z0-tol) & (z<=z1+tol)
            on_surf = (near_x&in_y&in_z) | (near_y&in_x&in_z) | (near_z&in_x&in_y)
        mk[interior & on_surf] = marker

    return mk


# ── Edge extraction ────────────────────────────────────────────────────────────

def _edges_from_cells(cells: np.ndarray) -> np.ndarray:
    """Extract all unique sorted edges from cell-to-vertex connectivity."""
    n = cells.shape[1]
    pairs = np.concatenate(
        [np.sort(cells[:, [i, j]], axis=1) for i, j in combinations(range(n), 2)]
    )
    return np.unique(pairs, axis=0).astype(np.int64)


# ── Cell-to-edge ───────────────────────────────────────────────────────────────

def _cell_to_edge(cells: np.ndarray, edges: np.ndarray, n_pts: int) -> np.ndarray:
    """
    cell_to_edge[c, e] = row index into 'edges' (vertex_to_vertex) of the
    e-th local edge of cell c. Local edge order follows
    itertools.combinations(range(d), 2) over the cell's d local vertex slots
    (matches the order _edges_from_cells uses to enumerate edges):
        triangle (d=3): (0,1), (0,2), (1,2)
        tetra    (d=4): (0,1), (0,2), (0,3), (1,2), (1,3), (2,3)

    Implementation: encode each undirected vertex pair (a,b) as a single
    integer key = min(a,b)*n_pts + max(a,b), sort the encoded 'edges' once,
    then binary-search each cell's local-edge keys against that sorted array
    -- O((n_cells*n_local + n_edges) log n_edges) overall, no per-cell dict.
    """
    d = cells.shape[1]
    local_pairs = list(combinations(range(d), 2))

    e_sorted   = np.sort(edges, axis=1).astype(np.int64)
    key_edges  = e_sorted[:, 0] * n_pts + e_sorted[:, 1]
    order      = np.argsort(key_edges)          # order[i] = original edge index
    key_sorted = key_edges[order]

    cell_to_edge = np.empty((len(cells), len(local_pairs)), dtype=np.int64)
    for li, (a, b) in enumerate(local_pairs):
        va, vb = cells[:, a].astype(np.int64), cells[:, b].astype(np.int64)
        key = np.minimum(va, vb) * n_pts + np.maximum(va, vb)
        pos = np.searchsorted(key_sorted, key)
        cell_to_edge[:, li] = order[pos]

    return cell_to_edge


# ── Cell volumes / areas ────────────────────────────────────────────────────────

def _cell_volumes(pts: np.ndarray, cells: np.ndarray, dim: int) -> np.ndarray:
    """
    Per-cell measure: triangle area (2-D) or tetrahedron volume (3-D).

    2-D (shoelace formula on the x,y components; z is 0 for these meshes):
        area = 0.5 * |x0(y1-y2) + x1(y2-y0) + x2(y0-y1)|

    3-D (scalar triple product):
        volume = |(p1-p0) · [(p2-p0) × (p3-p0)]| / 6
    """
    if dim == 2:
        p0, p1, p2 = pts[cells[:, 0], :2], pts[cells[:, 1], :2], pts[cells[:, 2], :2]
        return 0.5 * np.abs(
            p0[:, 0] * (p1[:, 1] - p2[:, 1]) +
            p1[:, 0] * (p2[:, 1] - p0[:, 1]) +
            p2[:, 0] * (p0[:, 1] - p1[:, 1])
        )
    else:
        p0 = pts[cells[:, 0]]; p1 = pts[cells[:, 1]]
        p2 = pts[cells[:, 2]]; p3 = pts[cells[:, 3]]
        cross = np.cross(p2 - p0, p3 - p0)
        return np.abs(np.einsum('ij,ij->i', p1 - p0, cross)) / 6.0


# ── Vertex volumes (barycentric dual, periodicity-enforced) ────────────────────

def _vertex_volumes(cells: np.ndarray, cell_vol: np.ndarray, n_pts: int) -> np.ndarray:
    """
    Per-vertex barycentric dual volume: each cell's area/volume is split
    equally among its own vertices via the barycenter -- 1/3 per vertex for
    triangles, 1/4 for tetrahedra -- then summed over every cell incident to
    a given vertex. Does NOT account for periodicity by itself; see
    _enforce_periodicity_vertex_vol, applied afterwards in main().
    """
    n_per = cells.shape[1]                       # 3 (triangle) | 4 (tetra)
    share = np.repeat(cell_vol, n_per) / n_per    # same order as cells.ravel()
    return np.bincount(cells.ravel(), weights=share, minlength=n_pts)


def _enforce_periodicity_vertex_vol(vertex_vol: np.ndarray, per_pairs: dict) -> np.ndarray:
    """
    Merge periodic vertex volumes so that vertices identified with each other
    via /periodicity (x_pairs, y_pairs, z_pairs) share a single combined dual
    volume -- the sum of their individually-computed /vertex_vol entries --
    rather than each carrying only its own local share.

    Uses a union-find over the vertices that actually appear in any pair
    (most vertices are untouched, since only mesh boundary vertices can be
    periodic). Pairs from different directions are merged transitively, so
    a vertex periodic in more than one direction at once (e.g. a domain
    corner under doubly/triply periodic boundaries) is correctly grouped
    with all of its images, not just the ones sharing a single axis.
    """
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        root = x
        while parent[root] != root:
            root = parent[root]
        while parent[x] != root:          # path compression
            parent[x], x = root, parent[x]
        return root

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for arr in per_pairs.values():
        for v_lo, v_hi in arr:
            union(int(v_lo), int(v_hi))

    if not parent:          # nothing periodic -- return unchanged
        return vertex_vol

    merged = vertex_vol.copy()
    groups: dict = {}
    for v in parent:
        groups.setdefault(find(v), []).append(v)

    for members in groups.values():
        total = vertex_vol[members].sum()
        merged[members] = total

    return merged


# ── Vertex-to-cell (CSR) ───────────────────────────────────────────────────────

def _vertex_to_cell(n_pts: int, cells: np.ndarray):
    """
    CSR vertex-to-cell map.
    cells incident to vertex v:  indices[ offsets[v] : offsets[v+1] ]
    """
    n_per = cells.shape[1]
    v_ids = cells.ravel()
    c_ids = np.repeat(np.arange(len(cells), dtype=np.int64), n_per)
    order = np.argsort(v_ids, kind="stable")
    c_ids, v_ids = c_ids[order], v_ids[order]
    counts  = np.bincount(v_ids, minlength=n_pts).astype(np.int64)
    offsets = np.zeros(n_pts+1, dtype=np.int64); offsets[1:] = np.cumsum(counts)
    return offsets, c_ids


# ── Periodic pair detection ────────────────────────────────────────────────────

def _periodic_pairs(pts: np.ndarray, box: list, periodic: list,
                    bnd_tol: float = 1e-9, match_tol: float = 1e-7) -> dict:
    """
    For each periodic direction d match vertices on face coord_d≈0 to
    vertices on face coord_d≈L_d by nearest-neighbour on remaining coords.
    """
    pairs: dict = {}
    d_pts = pts.shape[1]
    for d in range(len(periodic)):
        if not periodic[d]: continue
        L, dname = box[d], _DIRS[d]
        other = [i for i in range(d_pts) if i != d]
        lo = np.where(np.abs(pts[:, d])     < bnd_tol)[0]
        hi = np.where(np.abs(pts[:, d] - L) < bnd_tol)[0]
        if lo.size == 0 or hi.size == 0:
            print(f"  [warn] periodic_{dname}: no boundary vertices found")
            pairs[dname] = np.empty((0,2), dtype=np.int64); continue
        if lo.size != hi.size:
            print(f"  [warn] periodic_{dname}: {lo.size} lo vs {hi.size} hi vertices")
        lo_co, hi_co = pts[lo][:, other], pts[hi][:, other]
        matched, n_miss = [], 0
        for i, vi in enumerate(lo):
            dist = np.linalg.norm(hi_co - lo_co[i], axis=1); j = int(np.argmin(dist))
            if dist[j] < match_tol: matched.append([vi, hi[j]])
            else: n_miss += 1
        pairs[dname] = (np.array(matched, dtype=np.int64) if matched
                        else np.empty((0,2), dtype=np.int64))
        status = f"{len(matched):,} pairs" + (f", {n_miss} unmatched" if n_miss else "")
        print(f"  periodic_{dname}: {status}")
    return pairs


# ── HDF5 writer ────────────────────────────────────────────────────────────────

def write_hdf5(h5_path: str,
               pts3d: np.ndarray, cells: np.ndarray, edges: np.ndarray,
               c2e: np.ndarray, cell_vol: np.ndarray, vertex_vol: np.ndarray,
               v2c_off: np.ndarray, v2c_idx: np.ndarray,
               pt_mk: np.ndarray, edge_mk,
               periodic: list, per_pairs: dict,
               box: list, dim: int, obstacles=None) -> None:
    obstacles = obstacles or []
    kw = dict(compression="gzip", compression_opts=4)
    with h5py.File(h5_path, "w") as f:
        f.attrs.update({"mesh_type" : "triangle" if dim==2 else "tetra",
                        "dimension" : dim, "box_length": list(box),
                        "n_vertices": int(len(pts3d)), "n_cells": int(len(cells)),
                        "n_edges"   : int(len(edges))})
        f.create_dataset("points",           data=pts3d, **kw)
        f.create_dataset("cell_to_vertex",   data=cells, **kw)
        d_c2e = f.create_dataset("cell_to_edge", data=c2e, **kw)
        d_c2e.attrs["description"] = ("cell_to_edge[c,e] indexes a row of "
                                      "/vertex_to_vertex; local edge order = "
                                      "itertools.combinations of the cell's "
                                      "local vertex indices")
        d_vol = f.create_dataset("cell_vol", data=cell_vol, **kw)
        d_vol.attrs["description"] = ("triangle area (2-D) or tetrahedron "
                                      "volume (3-D) of each cell, same order "
                                      "as /cell_to_vertex")
        f.create_dataset("vertex_to_vertex", data=edges, **kw)
        g = f.create_group("vertex_to_cell")
        g.attrs["format"] = "CSR"
        g.attrs["description"] = "indices[offsets[v]:offsets[v+1]] = cells of vertex v"
        g.create_dataset("offsets", data=v2c_off, **kw)
        g.create_dataset("indices", data=v2c_idx, **kw)
        d_vv = f.create_dataset("vertex_vol", data=vertex_vol, **kw)
        d_vv.attrs["description"] = ("barycentric dual volume per vertex: sum "
                                     "over incident cells of cell_vol/(d+1) "
                                     "(1/3 tri, 1/4 tet); vertices linked via "
                                     "/periodicity are merged to share one "
                                     "combined value")
        b = f.create_group("boundary")
        b.create_dataset("point_markers", data=pt_mk, **kw)
        if edge_mk is not None:
            b.create_dataset("edge_markers", data=edge_mk, **kw)
        p = f.create_group("periodicity")
        for d_name, is_per in zip(_DIRS[:dim], periodic):
            p.attrs[f"periodic_{d_name}"] = bool(is_per)
        for d_name, arr in per_pairs.items():
            p.attrs[f"periodic_{d_name}"] = True
            if len(arr): p.create_dataset(f"{d_name}_pairs", data=arr, **kw)
        if obstacles:
            o = f.create_group("obstacles")
            o.attrs["count"] = len(obstacles)
            for i, obs in enumerate(obstacles):
                og = o.create_group(f"obstacle_{i}")
                og.attrs["shape"]  = obs["shape"]
                og.attrs["marker"] = obs["marker"]
                og.attrs["center"] = list(obs["center"])
                if obs["shape"] in ("circle", "sphere"):
                    og.attrs["radius"] = obs["radius"]
                else:
                    og.attrs["half_extents"] = list(obs["half_extents"])
    print(f"\n  HDF5 written → {h5_path}")
    print(f"    vertices  : {len(pts3d):>10,}")
    print(f"    cells     : {len(cells):>10,}")
    print(f"    edges     : {len(edges):>10,}")
    print(f"    cell_vol  : min={cell_vol.min():.6g}  max={cell_vol.max():.6g}")
    print(f"    vertex_vol: min={vertex_vol.min():.6g}  max={vertex_vol.max():.6g}"
          f"  sum={vertex_vol.sum():.6g}")
    if obstacles:
        print(f"    obstacles : {len(obstacles):>10,}")


# ── XDMF writer ────────────────────────────────────────────────────────────────

def write_xdmf(xdmf_path: str, h5_path: str,
               n_pts: int, n_cells: int, n_edges: int, dim: int) -> None:
    """
    Two grids in a Spatial Collection:
      'cells' — volume mesh (Triangle / Tetrahedron topology), with
                BoundaryMarker + VertexVolume (point data) and CellVolume
                (cell data)
      'edges' — wire skeleton (Polyline topology); commented out by default
                (see below) to avoid interior mesh edges rendering as clutter
                on top of the 'cells' surface -- uncomment the XML block in
                the .xdmf file if you want to inspect connectivity directly.
    Toggle visibility per grid in ParaView's Pipeline Browser.
    Obstacle holes need no special handling here — they simply have no
    cells, so ParaView renders the cavity automatically; obstacle boundary
    markers are visible via the BoundaryMarker point attribute (10+/100+).
    XDMF and HDF5 must reside in the same directory.
    """
    h5  = Path(h5_path).name
    tt  = "Triangle" if dim==2 else "Tetrahedron"
    nn  = 3 if dim==2 else 4
    geo = f"{n_pts} 3"
    def di(dims, dset, dtype="Int", prec=None):
        pr = f' Precision="{prec}"' if prec else ""
        return (f'<DataItem Format="HDF" DataType="{dtype}"{pr} '
                f'Dimensions="{dims}">{h5}:{dset}</DataItem>')
    cell_grid = (
        f'    <Grid Name="cells" GridType="Uniform">\n'
        f'      <Topology TopologyType="{tt}" NumberOfElements="{n_cells}">\n'
        f'        {di(f"{n_cells} {nn}", "/cell_to_vertex")}\n'
        f'      </Topology>\n'
        f'      <Geometry GeometryType="XYZ">\n'
        f'        {di(geo, "/points", "Float", 8)}\n'
        f'      </Geometry>\n'
        f'      <Attribute Name="BoundaryMarker" AttributeType="Scalar" Center="Node">\n'
        f'        {di(f"{n_pts}", "/boundary/point_markers")}\n'
        f'      </Attribute>\n'
        f'      <Attribute Name="VertexVolume" AttributeType="Scalar" Center="Node">\n'
        f'        {di(f"{n_pts}", "/vertex_vol", "Float", 8)}\n'
        f'      </Attribute>\n'
        f'      <Attribute Name="CellVolume" AttributeType="Scalar" Center="Cell">\n'
        f'        {di(f"{n_cells}", "/cell_vol", "Float", 8)}\n'
        f'      </Attribute>\n'
        f'    </Grid>')
    edge_grid = (
        f'    <Grid Name="edges" GridType="Uniform">\n'
        f'      <Topology TopologyType="Polyline" NumberOfElements="{n_edges}"\n'
        f'                NodesPerElement="2">\n'
        f'        {di(f"{n_edges} 2", "/vertex_to_vertex")}\n'
        f'      </Topology>\n'
        f'      <Geometry GeometryType="XYZ">\n'
        f'        {di(geo, "/points", "Float", 8)}\n'
        f'      </Geometry>\n'
        f'    </Grid>')
    xdmf = ('<?xml version="1.0" ?>\n'
            '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n'
            '<Xdmf Version="3.0">\n  <Domain>\n'
            '    <!-- Spatial collection: "cells" and "edges" share /points. -->\n'
            '    <!-- Toggle visibility per grid in ParaView Pipeline Browser. -->\n'
            '    <Grid Name="MeshCollection" GridType="Collection"\n'
            '          CollectionType="Spatial">\n'
            f'{cell_grid}\n'
            f'<!-- {edge_grid} -->\n'
            '    </Grid>\n  </Domain>\n</Xdmf>\n')
    Path(xdmf_path).write_text(xdmf)
    print(f"  XDMF written → {xdmf_path}")


# ── Provenance (appended to the .xdmf as XML comments) ──────────────────────────

def _xml_comment_safe(text: str) -> str:
    """
    Make arbitrary text safe to embed inside an XML comment <!-- ... -->.
    XML forbids '--' anywhere inside a comment, and forbids the comment's
    content from ending in '-' (which would otherwise yield '--->'). Loop
    until no '--' remains, since replacing one occurrence can create a new
    adjacent pair when the source has runs of 3+ hyphens.
    """
    safe = text
    while "--" in safe:
        safe = safe.replace("--", "- -")
    if safe.endswith("-"):
        safe += " "
    return safe


def _append_provenance(xdmf_path: str, yaml_path: str, yaml_text: str,
                       terminal_text: str, script_text: str) -> None:
    """
    Append three XML comment blocks to the end of the .xdmf file, after the
    closing </Xdmf> tag, so every .xdmf is a self-contained record of
    exactly how it was produced:
      1) the YAML config used for this run
      2) the terminal output produced while building this mesh
      3) the full source of writer.py itself
    """
    def block(label: str, content: str) -> str:
        return (f"\n\n<!-- ===== {label} =====\n"
                f"{_xml_comment_safe(content)}\n"
                f"===== end {label} ===== -->\n")
    with open(xdmf_path, "a") as f:
        f.write("\n")
        f.write(block(f"YAML CONFIG: {Path(yaml_path).name}", yaml_text))
        f.write(block("TERMINAL OUTPUT", terminal_text))
        f.write(block("WRITER.PY SOURCE", script_text))
    print(f"  Provenance appended → {xdmf_path}")


# ── Entry point ────────────────────────────────────────────────────────────────

def main(yaml_path: str) -> None:
    # Capture everything printed during this run (mirrored to the real
    # terminal via _Tee) so it can be embedded in the .xdmf file afterwards.
    _stdout_buf  = io.StringIO()
    _real_stdout = sys.stdout
    sys.stdout   = _Tee(_real_stdout, _stdout_buf)
    try:
        cfg = load_config(yaml_path)
        mt  = cfg["mesh_type"]
        per = cfg["periodic"]
        box = cfg["box_length"]
        dim = 2 if mt in _2D_TYPES else 3

        _validate_obstacles(cfg.get("obstacles", []), box, dim, mt)
        _validate_stretching(cfg["stretching_ratio"], cfg["stretching_symmetric"],
                             cfg["n_divisions"], dim, mt)

        print(f"Building '{mt}' mesh …")

        edge_mk = None          # only set for 2-D types that provide it
        obstacle_meta = []      # only set for 'triangle' / 'tetra'

        # ── dispatch to builder ──────────────────────────────────────────────
        if mt == "triangle":
            pts, cells, edges, pt_mk, edge_mk, obstacle_meta = _build_2d(cfg)
        elif mt == "tetra":
            pts, cells, edges, pt_mk, obstacle_meta = _build_3d(cfg)
        elif mt == "gridsplit":
            pts, cells = _build_gridsplit(cfg)
            edges  = _edges_from_cells(cells)
            pt_mk  = _pt_markers_2d(pts, *box)
            edge_mk = _edge_markers_2d(pts, edges, *box)
        elif mt == "equilateral":
            pts, cells = _build_equilateral(cfg)
            edges  = _edges_from_cells(cells)
            pt_mk  = _pt_markers_2d(pts, *box)
            edge_mk = _edge_markers_2d(pts, edges, *box)
        elif mt == "kuhn":
            pts, cells = _build_kuhn(cfg)
            edges  = _edges_from_cells(cells)
            pt_mk  = _pt_markers_3d(pts, *box)

        # ── pad 2-D coordinates to 3-D (z = 0) for uniform XDMF geometry ────
        if dim == 2:
            pts3d      = np.column_stack([pts, np.zeros(len(pts))])
            pts_for_per = pts              # use native 2-D coords for pairing
        else:
            pts3d = pts_for_per = pts

        print("Computing cell-to-edge connectivity …")
        c2e = _cell_to_edge(cells, edges, len(pts3d))

        print("Computing cell volumes …")
        cell_vol = _cell_volumes(pts3d, cells, dim)

        print("Computing vertex volumes …")
        vertex_vol = _vertex_volumes(cells, cell_vol, len(pts3d))

        print("Computing vertex-to-cell connectivity …")
        v2c_off, v2c_idx = _vertex_to_cell(len(pts3d), cells)

        print("Detecting periodic vertex pairs …")
        per_pairs = _periodic_pairs(pts_for_per, box, per)

        print("Enforcing periodicity on vertex volumes …")
        vertex_vol = _enforce_periodicity_vertex_vol(vertex_vol, per_pairs)

        hdf5_path = str(Path(yaml_path).with_suffix(".h5"))
        xdmf_path = str(Path(yaml_path).with_suffix(".xdmf"))

        write_hdf5(hdf5_path,
                   pts3d, cells, edges, c2e, cell_vol, vertex_vol,
                   v2c_off, v2c_idx, pt_mk, edge_mk,
                   per, per_pairs, box, dim, obstacle_meta)
        write_xdmf(xdmf_path, hdf5_path,
                   len(pts3d), len(cells), len(edges), dim)

        terminal_text = _stdout_buf.getvalue()
    finally:
        sys.stdout = _real_stdout   # always restore, even on error

    # ── append provenance: YAML config + terminal output + writer.py source ──
    yaml_text   = Path(yaml_path).read_text()
    script_text = Path(__file__).read_text()
    _append_provenance(xdmf_path, yaml_path, yaml_text, terminal_text, script_text)


def _process_yaml_or_dir(path_str: str) -> None:
    """
    CLI entry dispatcher: accepts either a single YAML config file or a
    directory. For a directory, every '*.yaml' file directly inside it
    (non-recursive) is processed in turn via main(), each producing its own
    '.h5'/'.xdmf' pair (see 'Output files' in the module docstring). A
    failure in one config is reported and does not stop the remaining
    configs in the directory from being processed; the process exits with
    a nonzero status if any config failed.
    """
    path = Path(path_str)

    if path.is_dir():
        yaml_files = sorted(path.glob("*.yaml"))
        if not yaml_files:
            sys.exit(f"No .yaml files found in directory: {path}")
        print(f"Found {len(yaml_files)} YAML config(s) in {path}\n")

        failures = []
        for i, yf in enumerate(yaml_files, 1):
            print(f"[{i}/{len(yaml_files)}] {yf.name}")
            try:
                main(str(yf))
            except Exception as exc:
                print(f"  [error] {yf.name}: {exc}")
                failures.append(yf.name)
            print()

        if failures:
            sys.exit(f"Completed with {len(failures)} failure(s): "
                     f"{', '.join(failures)}")
        return

    if not path.is_file():
        sys.exit(f"Path not found: {path}")
    main(str(path))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(f"Usage: {sys.argv[0]} config.yaml | config_dir/")
    _process_yaml_or_dir(sys.argv[1])