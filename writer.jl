#!/usr/bin/env julia
#=
Build a mesh of a box, in 2-D or in 3-D. This file holds five mesh
types. It needs no other file of this project, and it does not use the
Gmsh package. It replaces the three files writer.jl (the old version),
writer_delaunay.jl, and writer_delaunay_3d.jl. The file writer_gmsh.jl
stays separate. The YAML schema, the HDF5 layout, and the XDMF layout
match the layouts of writer_gmsh.jl.

Usage:
    julia writer.jl config.yaml
    julia writer.jl configs_dir/

The key mesh_type selects the method. The default is gridsplit.

    mesh_type    Dim  Method
    gridsplit    2    Grid of rectangles. Each rectangle has 2 triangles.
    equilateral  2    Offset-row lattice of near-equilateral triangles.
    kuhn         3    Grid of boxes. Each box has 6 tetrahedra.
    triangle     2    Delaunay refinement (Ruppert), with obstacles.
    tetra        3    Delaunay refinement (Shewchuk), with obstacles.

Keys in the YAML file:

    Key                 Structured types        triangle and tetra
    box_length          Required. dim entries.  Required. dim entries.
    periodic            Optional. dim entries.  Optional. dim entries.
                        Default: all false.     Default: all false.
    n_divisions         Required. dim entries.  Error.
    diagonal            gridsplit and kuhn.     Error.
                        Default: backslash (gridsplit), main (kuhn).
    stretching_ratio    gridsplit and kuhn.     Error.
    stretching_symmetric  gridsplit and kuhn.   Error.
    obstacles           Error if not empty.     Optional.
    mesh_options        Error.                  Required.
    hierarchy           Optional.               Optional.
    gmsh_options        Error.                  Error.

A config that has mesh_options and no mesh_type is an error. An error
for a key that does not fit the mesh type names the key and the type.

    diagonal for gridsplit: backslash, slash, alternate.
    diagonal for kuhn: main, flip_x, flip_y, flip_z, alternate_layers.

An obstacle record has the keys shape, center, and one more key.
    triangle: shape circle with radius, or shape rectangle with
              half_extents (2 entries).
    tetra:    shape sphere with radius, or shape box with
              half_extents (3 entries).
An obstacle that crosses or touches the box is an error. Two obstacles
with overlapping bounding boxes are an error.

Keys in mesh_options:
    mesh_size           Target edge length far from obstacles. Required.
    mesh_size_obstacle  Target edge length on obstacle boundaries.
                        Default: mesh_size.
    grading             Growth of the target size per unit of distance
                        from an obstacle. Default: 0.3.
    min_angle           triangle only. Smallest cell angle in degrees
                        that refinement enforces. The value must be in
                        (0, 20.7]. Default: 20.0.
    radius_edge_bound   tetra only. Largest allowed ratio of the
                        circumradius to the shortest edge. The value
                        must be 2.0 or more. Default: 2.0.
    min_dihedral        tetra only. Goal for the smallest dihedral angle,
                        in degrees. Flips of faces and edges try to
                        remove the cells below the goal. The goal is not
                        a guarantee. The value must be in [0, 60]. The
                        value 0 turns the flips off. Default: 10.0.
    max_dihedral        tetra only. Goal for the largest dihedral angle,
                        in degrees. The flips use it in the same way as
                        min_dihedral. The value must be 0 or in [90, 180).
                        The value 0 turns the goal off. Default: 0.
    smoothing_sweeps    Number of Jacobi smoothing sweeps. Default: 3.
    max_vertices        Safety limit on the vertex count.
                        Default: 2000000 for triangle, 1000000 for tetra.

Keys in hierarchy (a mesh hierarchy for geometric multigrid):
    levels              Number of levels, coarsest included. A whole
                        number, 1 or more. Required.
    boundary            frozen. The only value for now. The new boundary
                        vertices are edge midpoints, so an obstacle keeps
                        the polygon (polyhedron) of level 0.
    max_cells           Limit on the cells of the finest level. The
                        writer stops with an error before it refines
                        when the estimate is above the limit.
                        Default: 2000000.
The mesh that the other keys describe is level 0, the coarsest mesh.
Level k comes from k red refinements. Each triangle gives 4 triangles.
Each tetrahedron gives 8 tetrahedra.
The files are stem_L0.h5, stem_L1.h5, and so on, with one XDMF file
each. Every file has the group hierarchy. A file after level 0 has the
datasets fine_vertex_parents (n_vertices by 2) and cell_parent
(n_cells), 0-based. The fine vertices 1 to nV are the coarse vertices.
The midpoint of coarse edge e is fine vertex nV + e. The children of
coarse cell k are the rows c (k - 1) + 1 to c k. The number c is 4 in
2-D and 8 in 3-D. The writer makes no transfer operator. A tetrahedron
gives 8 tetrahedra (Bey). The writer picks the octahedron diagonal of
each cell with the max-min rule. It keeps the diagonal whose worst
interior child has the best mean-ratio shape.

Markers: the 2-D box sides are 1 bottom, 2 right, 3 top, 4 left. The 3-D
box faces are 1 xmin, 2 xmax, 3 ymin, 4 ymax, 5 zmin, 6 zmax. Obstacle
number i, counted from 0, has the marker 10 + i in 2-D and 100 + i in
3-D. An interior vertex has the marker 0. The lowest marker wins at a
shared vertex.

Output: one HDF5 file and one XDMF file next to each YAML file. The
attribute mesh_type of the HDF5 file names the cell type: triangle for
every 2-D mesh and tetra for every 3-D mesh, as in writer.py.
The attribute schema_version is 1, as in writer.py: row-major (N, k)
datasets and 0-based connectivity. It increases when a dataset layout or
a dataset meaning changes.
The dataset boundary/edge_markers exists for 2-D meshes only. The group
obstacles exists only when the config has obstacles.

Method of the Delaunay types (triangle and tetra):
    1. Build the boundary vertices and the boundary segments (and, in
       3-D, the facets: each sphere is an icosphere).
    2. Build a Delaunay triangulation of these vertices.
    3. Split segments (and subfacets) until each one is an edge (or a
       face) of the mesh and no vertex encroaches on it.
    4. Set a flag on each cell. A flood fill that stops at the boundary
       sets the flags.
    5. Insert circumcenters of bad cells until every cell in the domain
       has the required shape and size.
    6. Remove the cells outside the domain and renumber.
    7. Smooth the interior vertices with Jacobi sweeps.
Refinement does not remove sliver cells in 3-D. Smoothing improves
them, but it cannot move a boundary vertex.

Index convention: every Julia array here stays 1-based, the normal
Julia way. A stored connectivity value is 0-based in the HDF5 file, to
match the XDMF file. A row of cell_to_vertex, an entry in
vertex_to_vertex, and a CSR index are such values. The conversion
happens in main.

Names: a function of the 2-D Delaunay mesher ends with 2. A function of
the 3-D Delaunay mesher ends with 3. A kernel that works for one
dimension only has 2d or 3d in its name. Code that serves all
dimensions has no suffix.

The file has these parts, in this order.
    1. The computing kernels. They follow the house style in
       skill-stade: no keyword argument, no allocation, positional
       scalar-first indexing, only the allowed math functions, and a
       plain comment header instead of a docstring.
    2. The shared driver code: exact predicates, the queue, and the
       edge and vertex-to-cell helpers.
    3. The 2-D mesher core and the 3-D mesher core.
    4. The builders. Each builder gives one result record to main.
    5. The config check and the periodicity helpers.
    6. The input and output code, main, and the entry point.
Parts 1 to 5 lie between the two "package-free part" markers. A test
script can evaluate that text without the YAML and HDF5 packages.
=#

using YAML
using HDF5
using LinearAlgebra
using Random

# =====================================================================
# package-free part: begins
# =====================================================================

# Names of the axes.
_dirs() = ("x", "y", "z")

# The mesh types, grouped by how the mesh is built.
_structured_types() = ("gridsplit", "equilateral", "kuhn")
_delaunay_types() = ("triangle", "tetra")
_mesh_types() = ("gridsplit", "equilateral", "kuhn", "triangle", "tetra")

# The obstacle shapes that each dimension accepts.
_obstacle_shapes_2d() = ("circle", "rectangle")
_obstacle_shapes_3d() = ("sphere", "box")

# The vertex pairs of the edges of one triangle and of one tetrahedron,
# as local vertex numbers.
_tri_local_edges() = ((1, 2), (1, 3), (2, 3))
_tet_local_edges() = ((1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4))

# The 6 permutations of the 3 axes, 1-based, in the order a cube walks
# corner to corner for the Kuhn decomposition, and each one's parity.
_kuhn_perm_axis() = [1 2 3; 1 3 2; 2 1 3; 2 3 1; 3 1 2; 3 2 1]
_kuhn_perm_parity() = [1, -1, -1, 1, 1, -1]

# Print one line. Keep a copy in logBuf for the provenance record.
function _log(logBuf::IOBuffer, msg::AbstractString)
    println(msg)
    println(logBuf, msg)
end

# =====================================================================
# computing kernels -- house style from skill-stade
# =====================================================================

# ── Grid kernels (structured meshes) ──────────────────────────────────

# stretch_positions_kernel(n, ratio, symmetric, pos)
#
# Fill pos with n+1 breakpoints in [0,1], from a geometric-ratio
# grading, for later scaling by a box length.
#
# n: number of cells
# ratio: ratio between consecutive cell widths (1.0 gives uniform
#   spacing)
# symmetric: 0 grades from one side to the other. 1 grades from both
#   ends toward the middle. Each side has its own one-sided grading,
#   and the two sides mirror each other. With an odd n the right half
#   has one cell more. With n = 1 the grading is one-sided.
# pos: length n+1 output array; pos[i] is the i-th breakpoint, 1-based
function stretch_positions_kernel(n, ratio, symmetric, pos)
    if symmetric == 1 && n >= 2
        # With an odd n the right half gets the extra cell, as in writer.py.
        nLeft = div(n, 2)
        nRight = n - div(n, 2)
        wSumLeft = 0.0
        w = 1.0
        # total raw (unnormalized) left-side width -- carried sum
        for i = 1:nLeft
            wSumLeft += w
            w *= ratio
        end
        pos[1] = 0.0
        w = 1.0
        acc = 0.0
        for i = 1:nLeft
            acc += w
            pos[i + 1] = 0.5 * acc / wSumLeft
            w *= ratio
        end
        wSumRight = 0.0
        w = 1.0
        for i = 1:nRight
            wSumRight += w
            w *= ratio
        end
        pos[n + 1] = 1.0
        w = 1.0
        acc = 0.0
        # The right side starts at its own end and moves inward, toward
        # the midpoint that the left loop set. The last term would set
        # that midpoint again, so the loop leaves it out.
        for i = 1:nRight - 1
            acc += w
            pos[n - i + 1] = 1.0 - 0.5 * acc / wSumRight
            w *= ratio
        end
    else
        wSum = 0.0
        w = 1.0
        for i = 1:n
            wSum += w
            w *= ratio
        end
        pos[1] = 0.0
        w = 1.0
        acc = 0.0
        for i = 1:n
            acc += w
            pos[i + 1] = acc / wSum
            w *= ratio
        end
    end
    return nothing
end

# grid_points_2d_kernel(nx, ny, px, py, pts)
#
# Fill pts with every point of a tensor-product grid.
#
# nx: number of cells along x
# ny: number of cells along y
# px: length nx+1 breakpoint array along x
# py: length ny+1 breakpoint array along y
# pts: (nx+1)*(ny+1)-by-2 output array; row order is x fastest
function grid_points_2d_kernel(nx, ny, px, py, pts)
    nCols = nx + 1
    for idx = 1:(ny + 1) * nCols
        j = div(idx - 1, nCols)
        i = mod(idx - 1, nCols)
        pts[idx, 1] = px[i + 1]
        pts[idx, 2] = py[j + 1]
    end
    return nothing
end

# grid_triangles_2d_kernel(nx, ny, diagonalMode, cells)
#
# Fill cells with 2 triangles per grid cell, split along one
# diagonal, over the point grid grid_points_2d_kernel fills for the
# same nx, ny.
#
# nx: number of cells along x
# ny: number of cells along y
# diagonalMode: 0 splits every cell from bottom-right to top-left
#   (backslash). 1 splits every cell from bottom-left to top-right
#   (slash). 2 alternates by checkerboard parity. It uses backslash
#   where i + j is even and slash where i + j is odd.
# cells: (2*nx*ny)-by-3 output array, point rows 1-based
function grid_triangles_2d_kernel(nx, ny, diagonalMode, cells)
    nCols = nx + 1
    for idx = 1:nx * ny
        j = div(idx - 1, nx)
        i = mod(idx - 1, nx)
        bl = j * nCols + i + 1
        br = j * nCols + i + 2
        tl = (j + 1) * nCols + i + 1
        tr = (j + 1) * nCols + i + 2
        cellParity = mod(i + j, 2)
        useSlash = diagonalMode
        if diagonalMode == 2
            useSlash = cellParity
        end
        if useSlash == 1
            cells[2 * idx - 1, 1] = bl
            cells[2 * idx - 1, 2] = br
            cells[2 * idx - 1, 3] = tr
            cells[2 * idx, 1] = bl
            cells[2 * idx, 2] = tr
            cells[2 * idx, 3] = tl
        else
            cells[2 * idx - 1, 1] = bl
            cells[2 * idx - 1, 2] = br
            cells[2 * idx - 1, 3] = tl
            cells[2 * idx, 1] = br
            cells[2 * idx, 2] = tr
            cells[2 * idx, 3] = tl
        end
    end
    return nothing
end

# equilateral_points_kernel(nx, ny, Lx, Ly, pts)
#
# Fill pts with the points of an offset-row triangular lattice. An even
# row (0, 2, ...) has nx+1 points at x = 0, dx, ..., Lx. An odd row has
# nx+2 points: x = 0, then the midpoints dx/2, 3dx/2, ..., then x = Lx.
# The two extra points of an odd row lie on the left and right box
# sides. Row j lies at y = j*dy. The rows follow each other in pts.
#
# nx: number of columns of cells
# ny: number of rows of cells
# Lx: box length along x
# Ly: box length along y
# pts: output array with x in column 1 and y in column 2. It has
#   (ny div 2 + 1)*(nx+1) + ((ny+1) div 2)*(nx+2) rows.
function equilateral_points_kernel(nx, ny, Lx, Ly, pts)
    nEven = div(ny, 2) + 1
    nOdd = div(ny + 1, 2)
    # Even rows: j = 2r. Rows before row j hold j*(nx+1) + div(j, 2) points.
    for idx = 1:nEven * (nx + 1)
        r = div(idx - 1, nx + 1)
        i = mod(idx - 1, nx + 1)
        row = idx + r * (nx + 2)
        pts[row, 1] = Lx * i / nx
        pts[row, 2] = Ly * (2 * r) / ny
    end
    # Odd rows: j = 2r + 1. The row starts after (2r+1)*(nx+1) + r points.
    # Left boundary point of each odd row.
    for r = 0:nOdd - 1
        row = (2 * r + 1) * (nx + 1) + r + 1
        pts[row, 1] = 0.0
        pts[row, 2] = Ly * (2 * r + 1) / ny
    end
    # Interior (midpoint) points of each odd row.
    for idx = 1:nOdd * nx
        r = div(idx - 1, nx)
        i = mod(idx - 1, nx)
        row = (2 * r + 1) * (nx + 1) + r + 2 + i
        pts[row, 1] = Lx * (i + 0.5) / nx
        pts[row, 2] = Ly * (2 * r + 1) / ny
    end
    # Right boundary point of each odd row.
    for r = 0:nOdd - 1
        row = (2 * r + 1) * (nx + 1) + r + nx + 2
        pts[row, 1] = Lx
        pts[row, 2] = Ly * (2 * r + 1) / ny
    end
    return nothing
end

# equilateral_triangles_kernel(nx, ny, cells)
#
# Fill cells with the triangles of the lattice that
# equilateral_points_kernel fills. Each band of cells between two
# rows holds 2*nx+1 triangles. One right-angled triangle lies at each
# box side. nx triangles point one way and nx-1 point the other way.
# Band j starts at cell row j*(2*nx+1). The vertex order is not yet
# oriented, fix_orientation_kernel does that.
#
# nx: number of columns of cells
# ny: number of rows of cells
# cells: (ny*(2*nx+1))-by-3 output array, point rows 1-based
function equilateral_triangles_kernel(nx, ny, cells)
    per = 2 * nx + 1
    nEvenBands = div(ny + 1, 2)
    nOddBands = div(ny, 2)
    # Even band j = 2r. Its lower row is even and its upper row is odd.
    # The variable s0 holds the 0-based number of the first point of
    # the lower row. The variable s1 holds it for the upper row.
    for r = 0:nEvenBands - 1
        s0 = 2 * r * (nx + 1) + r
        s1 = (2 * r + 1) * (nx + 1) + r
        base = 2 * r * per
        cells[base + 1, 1] = s0 + 1
        cells[base + 1, 2] = s1 + 1
        cells[base + 1, 3] = s1 + 2
        cells[base + per, 1] = s0 + nx + 1
        cells[base + per, 2] = s1 + nx + 1
        cells[base + per, 3] = s1 + nx + 2
    end
    for idx = 1:nEvenBands * nx
        r = div(idx - 1, nx)
        i = mod(idx - 1, nx)
        s0 = 2 * r * (nx + 1) + r
        s1 = (2 * r + 1) * (nx + 1) + r
        base = 2 * r * per
        # triangle that points up
        cells[base + 2 + 2 * i, 1] = s0 + i + 1
        cells[base + 2 + 2 * i, 2] = s0 + i + 2
        cells[base + 2 + 2 * i, 3] = s1 + i + 2
    end
    for idx = 1:nEvenBands * (nx - 1)
        r = div(idx - 1, nx - 1)
        i = mod(idx - 1, nx - 1)
        s0 = 2 * r * (nx + 1) + r
        s1 = (2 * r + 1) * (nx + 1) + r
        base = 2 * r * per
        # triangle that points down
        cells[base + 3 + 2 * i, 1] = s0 + i + 2
        cells[base + 3 + 2 * i, 2] = s1 + i + 2
        cells[base + 3 + 2 * i, 3] = s1 + i + 3
    end
    # Odd band j = 2r+1, from an odd row to an even row (nx+1 points).
    for r = 0:nOddBands - 1
        s0 = (2 * r + 1) * (nx + 1) + r
        s1 = (2 * r + 2) * (nx + 1) + r + 1
        base = (2 * r + 1) * per
        cells[base + 1, 1] = s0 + 1
        cells[base + 1, 2] = s1 + 1
        cells[base + 1, 3] = s0 + 2
        cells[base + per, 1] = s0 + nx + 2
        cells[base + per, 2] = s0 + nx + 1
        cells[base + per, 3] = s1 + nx + 1
    end
    for idx = 1:nOddBands * nx
        r = div(idx - 1, nx)
        i = mod(idx - 1, nx)
        s0 = (2 * r + 1) * (nx + 1) + r
        s1 = (2 * r + 2) * (nx + 1) + r + 1
        base = (2 * r + 1) * per
        # triangle that points down
        cells[base + 2 + 2 * i, 1] = s0 + i + 2
        cells[base + 2 + 2 * i, 2] = s1 + i + 1
        cells[base + 2 + 2 * i, 3] = s1 + i + 2
    end
    for idx = 1:nOddBands * (nx - 1)
        r = div(idx - 1, nx - 1)
        i = mod(idx - 1, nx - 1)
        s0 = (2 * r + 1) * (nx + 1) + r
        s1 = (2 * r + 2) * (nx + 1) + r + 1
        base = (2 * r + 1) * per
        # triangle that points up
        cells[base + 3 + 2 * i, 1] = s0 + i + 2
        cells[base + 3 + 2 * i, 2] = s0 + i + 3
        cells[base + 3 + 2 * i, 3] = s1 + i + 2
    end
    return nothing
end

# grid_points_3d_kernel(nx, ny, nz, px, py, pz, pts)
#
# Fill pts with every point of a tensor-product grid.
#
# nx, ny, nz: number of cells along x, y, z
# px, py, pz: breakpoint arrays, length nx+1, ny+1, nz+1
# pts: (nx+1)*(ny+1)*(nz+1)-by-3 output array; x fastest, then y, then z
function grid_points_3d_kernel(nx, ny, nz, px, py, pz, pts)
    nCols = nx + 1
    nRows = ny + 1
    for idx = 1:(nx + 1) * (ny + 1) * (nz + 1)
        rem1 = mod(idx - 1, nCols * nRows)
        k = div(idx - 1, nCols * nRows)
        j = div(rem1, nCols)
        i = mod(rem1, nCols)
        pts[idx, 1] = px[i + 1]
        pts[idx, 2] = py[j + 1]
        pts[idx, 3] = pz[k + 1]
    end
    return nothing
end

# kuhn_tets_kernel(nx, ny, nz, diagonalMode, permAxis, permParity, cells)
#
# Fill cells with 6 Kuhn tetrahedra per cube, over the point grid
# grid_points_3d_kernel fills for the same nx, ny, nz. Each tet is
# already positively oriented. The sign of the walk from the offset
# corner is the parity of that corner times the parity of the
# permutation. The closed form below corrects the sign in place, so no
# later pass is needed.
#
# nx, ny, nz: number of cubes along x, y, z
# diagonalMode: 0 main, 1 flip_x, 2 flip_y, 3 flip_z, 4 alternate_layers
# permAxis: 6-by-3 array. Row p lists the axis (1, 2, or 3) that flips
#   at each of the 3 steps of permutation p, walking corner to corner.
# permParity: length-6 array; +1 for an even permutation, -1 for odd
# cells: (6*nx*ny*nz)-by-4 output array, point rows 1-based
function kuhn_tets_kernel(nx, ny, nz, diagonalMode, permAxis, permParity, cells)
    nCols = nx + 1
    nRows = ny + 1
    for idx = 1:nx * ny * nz * 6
        cube = div(idx - 1, 6)
        p = mod(idx - 1, 6) + 1
        k = div(cube, nx * ny)
        rem1 = mod(cube, nx * ny)
        j = div(rem1, nx)
        i = mod(rem1, nx)

        oX = 0; oY = 0; oZ = 0
        if diagonalMode == 1
            oX = 1
        end
        if diagonalMode == 2
            oY = 1
        end
        if diagonalMode == 3
            oZ = 1
        end
        if diagonalMode == 4
            oX = mod(i, 2)
            oY = mod(j, 2)
            oZ = mod(k, 2)
        end

        vx0 = oX; vy0 = oY; vz0 = oZ
        axis1 = permAxis[p, 1]
        vx1 = vx0; vy1 = vy0; vz1 = vz0
        if axis1 == 1
            vx1 = 1 - vx0
        end
        if axis1 == 2
            vy1 = 1 - vy0
        end
        if axis1 == 3
            vz1 = 1 - vz0
        end

        axis2 = permAxis[p, 2]
        vx2 = vx1; vy2 = vy1; vz2 = vz1
        if axis2 == 1
            vx2 = 1 - vx1
        end
        if axis2 == 2
            vy2 = 1 - vy1
        end
        if axis2 == 3
            vz2 = 1 - vz1
        end

        axis3 = permAxis[p, 3]
        vx3 = vx2; vy3 = vy2; vz3 = vz2
        if axis3 == 1
            vx3 = 1 - vx2
        end
        if axis3 == 2
            vy3 = 1 - vy2
        end
        if axis3 == 3
            vz3 = 1 - vz2
        end

        row0 = (k + vz0) * nCols * nRows + (j + vy0) * nCols + (i + vx0) + 1
        row1 = (k + vz1) * nCols * nRows + (j + vy1) * nCols + (i + vx1) + 1
        row2 = (k + vz2) * nCols * nRows + (j + vy2) * nCols + (i + vx2) + 1
        row3 = (k + vz3) * nCols * nRows + (j + vy3) * nCols + (i + vx3) + 1

        parityO = 1 - 2 * mod(oX + oY + oZ, 2)
        pp = permParity[p]
        if parityO * pp < 0
            cells[idx, 1] = row0
            cells[idx, 2] = row1
            cells[idx, 3] = row3
            cells[idx, 4] = row2
        else
            cells[idx, 1] = row0
            cells[idx, 2] = row1
            cells[idx, 3] = row2
            cells[idx, 4] = row3
        end
    end
    return nothing
end

# point_markers_grid_2d_kernel(nx, ny, ptMk)
#
# Fill ptMk with a box-side marker per point of the tensor grid
# grid_points_2d_kernel fills for the same nx, ny: 1 bottom, 2 right, 3 top, 4 left, 0 interior. The
# lowest marker wins at a corner.
#
# nx: number of cells along x
# ny: number of cells along y
# ptMk: (nx+1)*(ny+1)-length output array
function point_markers_grid_2d_kernel(nx, ny, ptMk)
    nCols = nx + 1
    bigMk = 1000000
    for idx = 1:(ny + 1) * nCols
        j = div(idx - 1, nCols)
        i = mod(idx - 1, nCols)
        mkBottom = bigMk
        if j == 0
            mkBottom = 1
        end
        mkRight = bigMk
        if i == nx
            mkRight = 2
        end
        mkTop = bigMk
        if j == ny
            mkTop = 3
        end
        mkLeft = bigMk
        if i == 0
            mkLeft = 4
        end
        best = min(mkBottom, mkRight)
        best = min(best, mkTop)
        best = min(best, mkLeft)
        mkFinal = best
        if best == bigMk
            mkFinal = 0
        end
        ptMk[idx] = mkFinal
    end
    return nothing
end

# point_markers_grid_3d_kernel(nx, ny, nz, ptMk)
#
# Fill ptMk with a box-face marker per point of the tensor grid
# grid_points_3d_kernel fills for the same nx, ny, nz: 1 xmin, 2
# xmax, 3 ymin, 4 ymax, 5 zmin, 6 zmax, 0 interior. The lowest marker
# wins at a shared edge or corner.
#
# nx, ny, nz: number of cells along x, y, z
# ptMk: (nx+1)*(ny+1)*(nz+1)-length output array
function point_markers_grid_3d_kernel(nx, ny, nz, ptMk)
    nCols = nx + 1
    nRows = ny + 1
    bigMk = 1000000
    for idx = 1:(nx + 1) * (ny + 1) * (nz + 1)
        rem1 = mod(idx - 1, nCols * nRows)
        k = div(idx - 1, nCols * nRows)
        j = div(rem1, nCols)
        i = mod(rem1, nCols)
        mkXmin = bigMk
        if i == 0
            mkXmin = 1
        end
        mkXmax = bigMk
        if i == nx
            mkXmax = 2
        end
        mkYmin = bigMk
        if j == 0
            mkYmin = 3
        end
        mkYmax = bigMk
        if j == ny
            mkYmax = 4
        end
        mkZmin = bigMk
        if k == 0
            mkZmin = 5
        end
        mkZmax = bigMk
        if k == nz
            mkZmax = 6
        end
        best = min(mkXmin, mkXmax)
        best = min(best, mkYmin)
        best = min(best, mkYmax)
        best = min(best, mkZmin)
        best = min(best, mkZmax)
        mkFinal = best
        if best == bigMk
            mkFinal = 0
        end
        ptMk[idx] = mkFinal
    end
    return nothing
end

# ── Shared geometry kernels ───────────────────────────────────────────

# signed_measure_2d_kernel(cells, pts, n, meas)
#
# Fill meas with the signed area of each triangle, by the shoelace
# formula. A positive value marks the standard vertex order.
#
# cells: n-by-3 cell-to-vertex array
# pts: point array, x then y
# n: number of cells
# meas: length-n output array
function signed_measure_2d_kernel(cells, pts, n, meas)
    for c = 1:n
        i0 = cells[c, 1]; i1 = cells[c, 2]; i2 = cells[c, 3]
        x0 = pts[i0, 1]; y0 = pts[i0, 2]
        x1 = pts[i1, 1]; y1 = pts[i1, 2]
        x2 = pts[i2, 1]; y2 = pts[i2, 2]
        meas[c] = 0.5 * (x0 * (y1 - y2) + x1 * (y2 - y0) + x2 * (y0 - y1))
    end
    return nothing
end

# signed_measure_3d_kernel(cells, pts, n, meas)
#
# Fill meas with the signed volume of each tetrahedron, by the scalar
# triple product. A positive value marks the standard vertex order.
#
# cells: n-by-4 cell-to-vertex array
# pts: point array, x, y, z
# n: number of cells
# meas: length-n output array
function signed_measure_3d_kernel(cells, pts, n, meas)
    for c = 1:n
        i0 = cells[c, 1]; i1 = cells[c, 2]; i2 = cells[c, 3]; i3 = cells[c, 4]
        ax = pts[i1, 1] - pts[i0, 1]; ay = pts[i1, 2] - pts[i0, 2]; az = pts[i1, 3] - pts[i0, 3]
        bx = pts[i2, 1] - pts[i0, 1]; by = pts[i2, 2] - pts[i0, 2]; bz = pts[i2, 3] - pts[i0, 3]
        cx = pts[i3, 1] - pts[i0, 1]; cy = pts[i3, 2] - pts[i0, 2]; cz = pts[i3, 3] - pts[i0, 3]
        crossx = by * cz - bz * cy
        crossy = bz * cx - bx * cz
        crossz = bx * cy - by * cx
        meas[c] = (ax * crossx + ay * crossy + az * crossz) / 6.0
    end
    return nothing
end

# fix_orientation_kernel(meas, n, cellSize, cells)
#
# Force every cell onto the standard, positive vertex order.
# Swapping a cell's last two local vertices always flips its
# measure's sign, for both a triangle and a tetrahedron.
#
# meas: length-n signed measure, from signed_measure_2d_kernel or
#   signed_measure_3d_kernel on the same cells
# n: number of cells
# cellSize: local vertices per cell, 3 for a triangle, 4 for a
#   tetrahedron
# cells: n-by-cellSize cell-to-vertex array, corrected in place
function fix_orientation_kernel(meas, n, cellSize, cells)
    last = cellSize
    secondLast = cellSize - 1
    for c = 1:n
        if meas[c] < 0
            tmp = cells[c, secondLast]
            cells[c, secondLast] = cells[c, last]
            cells[c, last] = tmp
        end
    end
    return nothing
end

# cell_vertex_normal_2d_kernel(cells, pts, n, normal)
#
# Fill normal with one outward edge-normal vector per cell, per local
# vertex: Rot90 of the opposite edge. The sum of the 3 vectors for
# one cell equals zero. The third, z, component is left at zero.
#
# cells: n-by-3 cell-to-vertex array
# pts: point array, x then y
# n: number of cells
# normal: n-by-3-by-3 output array (cell, local vertex, xyz)
function cell_vertex_normal_2d_kernel(cells, pts, n, normal)
    for c = 1:n
        i0 = cells[c, 1]; i1 = cells[c, 2]; i2 = cells[c, 3]
        x0 = pts[i0, 1]; y0 = pts[i0, 2]
        x1 = pts[i1, 1]; y1 = pts[i1, 2]
        x2 = pts[i2, 1]; y2 = pts[i2, 2]
        e0x = x2 - x1; e0y = y2 - y1
        e1x = x0 - x2; e1y = y0 - y2
        e2x = x1 - x0; e2y = y1 - y0
        normal[c, 1, 1] = e0y; normal[c, 1, 2] = -e0x; normal[c, 1, 3] = 0.0
        normal[c, 2, 1] = e1y; normal[c, 2, 2] = -e1x; normal[c, 2, 3] = 0.0
        normal[c, 3, 1] = e2y; normal[c, 3, 2] = -e2x; normal[c, 3, 3] = 0.0
    end
    return nothing
end

# cell_vertex_normal_3d_kernel(cells, pts, n, normal)
#
# Fill normal with one dual-face-normal vector per cell, per local
# vertex, matching the positive-volume convention. The sum of the 4
# vectors for one cell equals zero.
#
# cells: n-by-4 cell-to-vertex array
# pts: point array, x, y, z
# n: number of cells
# normal: n-by-4-by-3 output array (cell, local vertex, xyz)
function cell_vertex_normal_3d_kernel(cells, pts, n, normal)
    for c = 1:n
        i0 = cells[c, 1]; i1 = cells[c, 2]; i2 = cells[c, 3]; i3 = cells[c, 4]
        ax = pts[i1, 1] - pts[i0, 1]; ay = pts[i1, 2] - pts[i0, 2]; az = pts[i1, 3] - pts[i0, 3]
        bx = pts[i2, 1] - pts[i0, 1]; by = pts[i2, 2] - pts[i0, 2]; bz = pts[i2, 3] - pts[i0, 3]
        cx = pts[i3, 1] - pts[i0, 1]; cy = pts[i3, 2] - pts[i0, 2]; cz = pts[i3, 3] - pts[i0, 3]
        s1x = 0.5 * (cy * bz - cz * by); s1y = 0.5 * (cz * bx - cx * bz); s1z = 0.5 * (cx * by - cy * bx)
        s2x = 0.5 * (ay * cz - az * cy); s2y = 0.5 * (az * cx - ax * cz); s2z = 0.5 * (ax * cy - ay * cx)
        s3x = 0.5 * (by * az - bz * ay); s3y = 0.5 * (bz * ax - bx * az); s3z = 0.5 * (bx * ay - by * ax)
        normal[c, 2, 1] = s1x; normal[c, 2, 2] = s1y; normal[c, 2, 3] = s1z
        normal[c, 3, 1] = s2x; normal[c, 3, 2] = s2y; normal[c, 3, 3] = s2z
        normal[c, 4, 1] = s3x; normal[c, 4, 2] = s3y; normal[c, 4, 3] = s3z
        normal[c, 1, 1] = -(s1x + s2x + s3x)
        normal[c, 1, 2] = -(s1y + s2y + s3y)
        normal[c, 1, 3] = -(s1z + s2z + s3z)
    end
    return nothing
end

# triangle_angle_kernel(cells, pts, n, minAng, maxAng)
#
# Fill minAng and maxAng with the smallest and the largest interior angle
# of each triangle, in degrees. The law of cosines gives each angle from
# the squared edge lengths. The cosine falls as the angle grows, so the
# largest cosine gives the smallest angle and the smallest cosine gives
# the largest angle. Only two acos calls are needed per cell.
#
# cells: n-by-3 cell-to-vertex array
# pts: point array, x then y
# n: number of cells
# minAng: length-n output array, smallest angle of each cell
# maxAng: length-n output array, largest angle of each cell
function triangle_angle_kernel(cells, pts, n, minAng, maxAng)
    for c = 1:n
        i0 = cells[c, 1]; i1 = cells[c, 2]; i2 = cells[c, 3]
        x0 = pts[i0, 1]; y0 = pts[i0, 2]
        x1 = pts[i1, 1]; y1 = pts[i1, 2]
        x2 = pts[i2, 1]; y2 = pts[i2, 2]
        # squared length of the edge opposite each vertex
        l0 = (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1)
        l1 = (x0 - x2) * (x0 - x2) + (y0 - y2) * (y0 - y2)
        l2 = (x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0)
        cos0 = (l1 + l2 - l0) / (2.0 * sqrt(l1 * l2))
        cos1 = (l0 + l2 - l1) / (2.0 * sqrt(l0 * l2))
        cos2 = (l0 + l1 - l2) / (2.0 * sqrt(l0 * l1))
        cosMax = min(1.0, max(cos0, max(cos1, cos2)))
        cosMin = max(-1.0, min(cos0, min(cos1, cos2)))
        minAng[c] = 57.29577951308232 * acos(cosMax)
        maxAng[c] = 57.29577951308232 * acos(cosMin)
    end
    return nothing
end

# vertex_volume_kernel(cells, cellVol, n, nPerCell, vertexVol)
#
# Split each cell's volume into an equal share for each of its own
# vertices. Then scatter-add the shares onto vertexVol. The target
# vertex differs by cell. This is a commutative accumulation, not a
# value-carrying loop. It still runs as one independent pass, one
# thread per (cell, local vertex) pair, per skill-stade rule 1.
#
# cells: n-by-nPerCell cell-to-vertex array
# cellVol: length-n array of cell volumes
# n: number of cells
# nPerCell: vertices per cell, 3 for a triangle, 4 for a tetrahedron
# vertexVol: output array, one entry per point; the caller must
#   pre-zero it before calling
function vertex_volume_kernel(cells, cellVol, n, nPerCell, vertexVol)
    for idx = 1:n * nPerCell
        c = div(idx - 1, nPerCell) + 1
        k = mod(idx - 1, nPerCell) + 1
        v = cells[c, k]
        share = cellVol[c] / nPerCell
        vertexVol[v] += share
    end
    return nothing
end

# count_negative_kernel(meas, n, out)
#
# Count the entries of meas that are below zero. The count does not
# depend on the order of the iterations.
#
# meas: length-n array of signed areas
# n: number of entries
# out: length-1 Int array. The value out[1] must be zero before the call.
function count_negative_kernel(meas, n, out)
    for c = 1:n
        if meas[c] < 0.0
            out[1] += 1
        end
    end
    return nothing
end

# ── 2-D Delaunay kernels ──────────────────────────────────────────────

# size_field_2d_kernel(px, py, n, obsKind, obsCx, obsCy, obsHx, obsHy,
#                      nObs, hMax, hObs, grading, h)
#
# Compute the target size at each query point. The target size is hObs
# on an obstacle boundary. It grows with the distance to the nearest
# obstacle boundary. It never exceeds hMax.
#
# px: x coordinate of each query point
# py: y coordinate of each query point
# n: number of query points
# obsKind: length-nObs array. 1 means circle. 2 means rectangle.
# obsCx: x coordinate of each obstacle center
# obsCy: y coordinate of each obstacle center
# obsHx: radius of a circle, or half the x extent of a rectangle
# obsHy: radius of a circle, or half the y extent of a rectangle
# nObs: number of obstacles
# hMax: target size far from obstacles
# hObs: target size on obstacle boundaries
# grading: growth of the target size per unit of distance
# h: length-n output array
function size_field_2d_kernel(px, py, n, obsKind, obsCx, obsCy, obsHx, obsHy, nObs, hMax, hObs, grading, h)
    for i = 1:n
        x = px[i]
        y = py[i]
        d = 1.0e30
        # The loop over obstacles carries a running minimum.
        for k = 1:nObs
            dx = x - obsCx[k]
            dy = y - obsCy[k]
            dk = abs(sqrt(dx * dx + dy * dy) - obsHx[k])
            if obsKind[k] == 2
                qx = abs(dx) - obsHx[k]
                qy = abs(dy) - obsHy[k]
                ox = max(qx, 0.0)
                oy = max(qy, 0.0)
                dk = abs(sqrt(ox * ox + oy * oy) + min(max(qx, qy), 0.0))
            end
            d = min(d, dk)
        end
        h[i] = min(hMax, hObs + grading * d)
    end
    return nothing
end

# box_side_vertices_kernel(nx, ny, px, py, vx, vy)
#
# Fill the vertices on the four box sides. The order is counterclockwise.
# It starts at the bottom-left corner. Each corner belongs to the side
# that starts at that corner. The kernel writes 2*nx + 2*ny vertices.
#
# nx: number of segments on the bottom side and on the top side
# ny: number of segments on the right side and on the left side
# px: length nx+1 breakpoint array along x, scaled to the box
# py: length ny+1 breakpoint array along y, scaled to the box
# vx: output array for the x coordinates
# vy: output array for the y coordinates
function box_side_vertices_kernel(nx, ny, px, py, vx, vy)
    # bottom side, from left to right
    for i = 1:nx
        vx[i] = px[i]
        vy[i] = py[1]
    end
    # right side, from bottom to top
    for j = 1:ny
        idx = nx + j
        vx[idx] = px[nx + 1]
        vy[idx] = py[j]
    end
    # top side, from right to left
    for i = 1:nx
        idx = nx + ny + i
        vx[idx] = px[nx - i + 2]
        vy[idx] = py[ny + 1]
    end
    # left side, from top to bottom
    for j = 1:ny
        idx = 2 * nx + ny + j
        vx[idx] = px[1]
        vy[idx] = py[ny - j + 2]
    end
    return nothing
end

# circle_vertices_kernel(cx, cy, r, nSeg, offset, vx, vy)
#
# Fill the vertices of one circle, in counterclockwise order. The first
# vertex lies on the positive x direction from the center.
#
# cx: x coordinate of the center
# cy: y coordinate of the center
# r: radius
# nSeg: number of vertices, equal to the number of segments
# offset: number of entries in vx and vy that come before this circle
# vx: output array for the x coordinates
# vy: output array for the y coordinates
function circle_vertices_kernel(cx, cy, r, nSeg, offset, vx, vy)
    twoPi = 6.283185307179586
    for k = 1:nSeg
        theta = twoPi * (k - 1) / nSeg
        vx[offset + k] = cx + r * cos(theta)
        vy[offset + k] = cy + r * sin(theta)
    end
    return nothing
end

# rectangle_vertices_kernel(cx, cy, hx, hy, mx, my, offset, vx, vy)
#
# Fill the vertices of one rectangle, in counterclockwise order. The
# first vertex is the bottom-left corner. Each corner belongs to the
# side that starts at that corner. The kernel writes 2*mx + 2*my vertices.
#
# cx: x coordinate of the center
# cy: y coordinate of the center
# hx: half the x extent
# hy: half the y extent
# mx: number of segments on the bottom side and on the top side
# my: number of segments on the right side and on the left side
# offset: number of entries in vx and vy that come before this rectangle
# vx: output array for the x coordinates
# vy: output array for the y coordinates
function rectangle_vertices_kernel(cx, cy, hx, hy, mx, my, offset, vx, vy)
    x0 = cx - hx
    x1 = cx + hx
    y0 = cy - hy
    y1 = cy + hy
    # bottom side, from left to right
    for i = 1:mx
        vx[offset + i] = x0 + (x1 - x0) * (i - 1) / mx
        vy[offset + i] = y0
    end
    # right side, from bottom to top
    for j = 1:my
        idx = offset + mx + j
        vx[idx] = x1
        vy[idx] = y0 + (y1 - y0) * (j - 1) / my
    end
    # top side, from right to left
    for i = 1:mx
        idx = offset + mx + my + i
        vx[idx] = x1 - (x1 - x0) * (i - 1) / mx
        vy[idx] = y1
    end
    # left side, from top to bottom
    for j = 1:my
        idx = offset + 2 * mx + my + j
        vx[idx] = x0
        vy[idx] = y1 - (y1 - y0) * (j - 1) / my
    end
    return nothing
end

# cell_circumcenter_2d_kernel(cellV, vx, vy, ids, nIds, ccx, ccy, ccR)
#
# Compute the circumcenter and the circumradius of each listed cell.
# Entry k of the outputs belongs to cell ids[k].
#
# cellV: cell-to-vertex array with 3 columns
# vx: x coordinate of each vertex
# vy: y coordinate of each vertex
# ids: list of cell numbers
# nIds: number of entries in ids
# ccx: output array for the x coordinate of the circumcenter
# ccy: output array for the y coordinate of the circumcenter
# ccR: output array for the circumradius
function cell_circumcenter_2d_kernel(cellV, vx, vy, ids, nIds, ccx, ccy, ccR)
    for k = 1:nIds
        c = ids[k]
        i0 = cellV[c, 1]; i1 = cellV[c, 2]; i2 = cellV[c, 3]
        x0 = vx[i0]
        y0 = vy[i0]
        # Move vertex 0 to the origin. Then solve two linear equations.
        bx = vx[i1] - x0
        by = vy[i1] - y0
        cx = vx[i2] - x0
        cy = vy[i2] - y0
        d = 2.0 * (bx * cy - by * cx)
        b2 = bx * bx + by * by
        c2 = cx * cx + cy * cy
        ux = (cy * b2 - by * c2) / d
        uy = (bx * c2 - cx * b2) / d
        ccx[k] = x0 + ux
        ccy[k] = y0 + uy
        ccR[k] = sqrt(ux * ux + uy * uy)
    end
    return nothing
end

# segment_encroach_2d_kernel(segV, vx, vy, ids, tx, ty, nIds, flag)
#
# Test whether a point encroaches on a segment. A point encroaches on a
# segment when it lies strictly inside the circle that has the segment
# as its diameter. This is the case when the two vectors from the point
# to the two segment ends have a negative dot product.
#
# segV: segment-to-vertex array with 2 columns
# vx: x coordinate of each vertex
# vy: y coordinate of each vertex
# ids: list of segment numbers
# tx: x coordinate of the test point for each entry
# ty: y coordinate of the test point for each entry
# nIds: number of entries in ids
# flag: output array. 1 means the point encroaches on the segment.
function segment_encroach_2d_kernel(segV, vx, vy, ids, tx, ty, nIds, flag)
    for k = 1:nIds
        s = ids[k]
        a = segV[s, 1]
        b = segV[s, 2]
        ax = vx[a] - tx[k]
        ay = vy[a] - ty[k]
        bx = vx[b] - tx[k]
        by = vy[b] - ty[k]
        flag[k] = 0
        if ax * bx + ay * by < 0.0
            flag[k] = 1
        end
    end
    return nothing
end

# segment_split_point_2d_kernel(segV, segCurve, vx, vy, obsCx, obsCy, obsHx, ids, nIds, sx, sy)
#
# Compute the split point of each listed segment. A straight segment
# splits at its midpoint. A segment of a circle splits at its midpoint
# moved onto the circle. The refined boundary then stays an inscribed
# polygon.
#
# segV: segment-to-vertex array with 2 columns
# segCurve: for each segment, 0 if the segment is straight, or the
#   number of the circle that carries the segment
# vx: x coordinate of each vertex
# vy: y coordinate of each vertex
# obsCx: x coordinate of each obstacle center
# obsCy: y coordinate of each obstacle center
# obsHx: radius of each circle
# ids: list of segment numbers
# nIds: number of entries in ids
# sx: output array for the x coordinate of the split point
# sy: output array for the y coordinate of the split point
function segment_split_point_2d_kernel(segV, segCurve, vx, vy, obsCx, obsCy, obsHx, ids, nIds, sx, sy)
    for k = 1:nIds
        s = ids[k]
        a = segV[s, 1]
        b = segV[s, 2]
        mx = 0.5 * (vx[a] + vx[b])
        my = 0.5 * (vy[a] + vy[b])
        sx[k] = mx
        sy[k] = my
        curve = segCurve[s]
        if curve > 0
            cx = obsCx[curve]
            cy = obsCy[curve]
            rad = obsHx[curve]
            dx = mx - cx
            dy = my - cy
            dist = sqrt(dx * dx + dy * dy)
            scale = rad / dist
            sx[k] = cx + dx * scale
            sy[k] = cy + dy * scale
        end
    end
    return nothing
end

# point_markers_2d_domain_kernel(pts, n, Lx, Ly, tol, obsKind, obsCx,
#                                obsCy, obsHx, obsHy, nObs, ptMk)
#
# Set the marker of each vertex from its coordinates. The markers are
# 1 bottom, 2 right, 3 top, 4 left, and 10 + i for obstacle number i,
# counted from 0. An interior vertex gets 0. Each write replaces the
# earlier one, and the writes go from the highest marker to the lowest.
# The lowest marker wins at a shared vertex, and a box side wins over
# an obstacle.
#
# pts: vertex array with x in column 1 and y in column 2
# n: number of vertices
# Lx: box length along x
# Ly: box length along y
# tol: distance below which a vertex lies on a boundary
# obsKind: 1 means circle. 2 means rectangle.
# obsCx: x coordinate of each obstacle center
# obsCy: y coordinate of each obstacle center
# obsHx: radius of a circle, or half the x extent of a rectangle
# obsHy: radius of a circle, or half the y extent of a rectangle
# nObs: number of obstacles
# ptMk: length-n output array
function point_markers_2d_domain_kernel(pts, n, Lx, Ly, tol, obsKind, obsCx, obsCy, obsHx, obsHy, nObs, ptMk)
    for v = 1:n
        x = pts[v, 1]
        y = pts[v, 2]
        mk = 0
        # Obstacles first. The highest obstacle number goes first. This
        # loop carries a value, because the last write wins.
        for k = 1:nObs
            i = nObs - k + 1
            dx = x - obsCx[i]
            dy = y - obsCy[i]
            val = abs(sqrt(dx * dx + dy * dy) - obsHx[i])
            if obsKind[i] == 2
                val = abs(max(abs(dx) - obsHx[i], abs(dy) - obsHy[i]))
            end
            if val < tol
                mk = 9 + i
            end
        end
        # Box sides last. Marker 1 is the last write, so it wins.
        if abs(x) < tol
            mk = 4
        end
        if abs(y - Ly) < tol
            mk = 3
        end
        if abs(x - Lx) < tol
            mk = 2
        end
        if abs(y) < tol
            mk = 1
        end
        ptMk[v] = mk
    end
    return nothing
end

# edge_markers_2d_domain_kernel(pts, edges, isBnd, nEdges, Lx, Ly, tol,
#                               obsKind, obsCx, obsCy, obsHx, obsHy, nObs, edgeMk)
#
# Set the marker of each boundary edge. An edge gets the marker of an
# entity when both of its end vertices lie on that entity. The order of
# the writes is the same as in point_markers_2d_domain_kernel. An edge
# that is not a boundary edge gets 0.
#
# pts: vertex array with x in column 1 and y in column 2
# edges: nEdges-by-2 edge-to-vertex array, 1-based
# isBnd: 1 for a boundary edge, 0 for an interior edge
# nEdges: number of edges
# Lx: box length along x
# Ly: box length along y
# tol: distance below which a vertex lies on a boundary
# obsKind: 1 means circle. 2 means rectangle.
# obsCx: x coordinate of each obstacle center
# obsCy: y coordinate of each obstacle center
# obsHx: radius of a circle, or half the x extent of a rectangle
# obsHy: radius of a circle, or half the y extent of a rectangle
# nObs: number of obstacles
# edgeMk: length-nEdges output array
function edge_markers_2d_domain_kernel(pts, edges, isBnd, nEdges, Lx, Ly, tol, obsKind, obsCx, obsCy, obsHx, obsHy, nObs, edgeMk)
    for e = 1:nEdges
        u = edges[e, 1]
        w = edges[e, 2]
        xu = pts[u, 1]
        yu = pts[u, 2]
        xw = pts[w, 1]
        yw = pts[w, 2]
        mk = 0
        if isBnd[e] == 1
            # This loop carries a value, because the last write wins.
            for k = 1:nObs
                i = nObs - k + 1
                dxu = xu - obsCx[i]
                dyu = yu - obsCy[i]
                dxw = xw - obsCx[i]
                dyw = yw - obsCy[i]
                valU = abs(sqrt(dxu * dxu + dyu * dyu) - obsHx[i])
                valW = abs(sqrt(dxw * dxw + dyw * dyw) - obsHx[i])
                if obsKind[i] == 2
                    valU = abs(max(abs(dxu) - obsHx[i], abs(dyu) - obsHy[i]))
                    valW = abs(max(abs(dxw) - obsHx[i], abs(dyw) - obsHy[i]))
                end
                # Both end vertices must lie on the obstacle.
                if max(valU, valW) < tol
                    mk = 9 + i
                end
            end
            if max(abs(xu), abs(xw)) < tol
                mk = 4
            end
            if max(abs(yu - Ly), abs(yw - Ly)) < tol
                mk = 3
            end
            if max(abs(xu - Lx), abs(xw - Lx)) < tol
                mk = 2
            end
            if max(abs(yu), abs(yw)) < tol
                mk = 1
            end
        end
        edgeMk[e] = mk
    end
    return nothing
end

# ── 3-D Delaunay kernels ──────────────────────────────────────────────

# size_field_3d_kernel(px, py, pz, n, obsKind, obsCx, obsCy, obsCz,
#                      obsHx, obsHy, obsHz, nObs, hMax, hObs, grading, h)
#
# Compute the target size at each query point. The target size is hObs
# on an obstacle surface. It grows with the distance to the nearest
# obstacle surface. It never exceeds hMax.
#
# px: x coordinate of each query point
# py: y coordinate of each query point
# pz: z coordinate of each query point
# n: number of query points
# obsKind: length-nObs array. 1 means sphere. 2 means box.
# obsCx: x coordinate of each obstacle center
# obsCy: y coordinate of each obstacle center
# obsCz: z coordinate of each obstacle center
# obsHx: radius of a sphere, or half the x extent of a box
# obsHy: radius of a sphere, or half the y extent of a box
# obsHz: radius of a sphere, or half the z extent of a box
# nObs: number of obstacles
# hMax: target size far from obstacles
# hObs: target size on obstacle surfaces
# grading: growth of the target size per unit of distance
# h: length-n output array
function size_field_3d_kernel(px, py, pz, n, obsKind, obsCx, obsCy, obsCz, obsHx, obsHy, obsHz, nObs, hMax, hObs, grading, h)
    for i = 1:n
        x = px[i]
        y = py[i]
        z = pz[i]
        d = 1.0e30
        # The loop over obstacles carries a running minimum.
        for k = 1:nObs
            dx = x - obsCx[k]
            dy = y - obsCy[k]
            dz = z - obsCz[k]
            dk = abs(sqrt(dx * dx + dy * dy + dz * dz) - obsHx[k])
            if obsKind[k] == 2
                qx = abs(dx) - obsHx[k]
                qy = abs(dy) - obsHy[k]
                qz = abs(dz) - obsHz[k]
                ox = max(qx, 0.0)
                oy = max(qy, 0.0)
                oz = max(qz, 0.0)
                dk = abs(sqrt(ox * ox + oy * oy + oz * oz) + min(max(max(qx, qy), qz), 0.0))
            end
            d = min(d, dk)
        end
        h[i] = min(hMax, hObs + grading * d)
    end
    return nothing
end

# line_points_kernel(x0, y0, z0, x1, y1, z1, s, i0, i1, offset, vx, vy, vz)
#
# Fill the points along a straight line. Point k of the line has the
# parameter s[k] in [0,1]. The kernel writes the points for k = i0 up to
# k = i1. The caller keeps the two end points out of the range, so they
# stay exact.
#
# x0, y0, z0: first end of the line
# x1, y1, z1: second end of the line
# s: breakpoint array in [0,1]
# i0: index of the first breakpoint to use
# i1: index of the last breakpoint to use
# offset: number of entries in vx, vy, vz that come before these points
# vx, vy, vz: output arrays for the coordinates
function line_points_kernel(x0, y0, z0, x1, y1, z1, s, i0, i1, offset, vx, vy, vz)
    for k = i0:i1
        t = s[k]
        idx = offset + k - i0 + 1
        vx[idx] = x0 + t * (x1 - x0)
        vy[idx] = y0 + t * (y1 - y0)
        vz[idx] = z0 + t * (z1 - z0)
    end
    return nothing
end

# rect_grid_points_kernel(ox, oy, oz, ux, uy, uz, wx, wy, wz, s1, s2,
#                         n1, n2, offset, vx, vy, vz)
#
# Fill the interior grid points of one rectangle. The rectangle has the
# corner o and the two side vectors u and w. The grid has n1 cells along
# u and n2 cells along w. The kernel writes (n1-1)*(n2-1) points. The
# point with the grid indices i and j, counted from 2, goes to the entry
# offset + (j-2)*(n1-1) + (i-1).
#
# ox, oy, oz: corner of the rectangle
# ux, uy, uz: side vector along the first axis of the grid
# wx, wy, wz: side vector along the second axis of the grid
# s1: length n1+1 breakpoint array along u
# s2: length n2+1 breakpoint array along w
# n1: number of cells along u
# n2: number of cells along w
# offset: number of entries in vx, vy, vz that come before these points
# vx, vy, vz: output arrays for the coordinates
function rect_grid_points_kernel(ox, oy, oz, ux, uy, uz, wx, wy, wz, s1, s2, n1, n2, offset, vx, vy, vz)
    for idx = 1:(n1 - 1) * (n2 - 1)
        i = mod(idx - 1, n1 - 1) + 2
        j = div(idx - 1, n1 - 1) + 2
        a = s1[i]
        b = s2[j]
        vx[offset + idx] = ox + a * ux + b * wx
        vy[offset + idx] = oy + a * uy + b * wy
        vz[offset + idx] = oz + a * uz + b * wz
    end
    return nothing
end

# sphere_project_kernel(cx, cy, cz, r, n, rx, ry, rz, offset, vx, vy, vz)
#
# Move n direction vectors onto a sphere. Vertex k gets the position
# center + r * (direction / length of direction).
#
# cx, cy, cz: center of the sphere
# r: radius of the sphere
# n: number of vectors
# rx, ry, rz: direction vectors
# offset: number of entries in vx, vy, vz that come before these vertices
# vx, vy, vz: output arrays for the coordinates
function sphere_project_kernel(cx, cy, cz, r, n, rx, ry, rz, offset, vx, vy, vz)
    for k = 1:n
        len = sqrt(rx[k] * rx[k] + ry[k] * ry[k] + rz[k] * rz[k])
        scale = r / len
        vx[offset + k] = cx + rx[k] * scale
        vy[offset + k] = cy + ry[k] * scale
        vz[offset + k] = cz + rz[k] * scale
    end
    return nothing
end

# cell_circumcenter_3d_kernel(cellV, vx, vy, vz, ids, nIds, ccx, ccy, ccz, ccR)
#
# Compute the circumcenter and the circumradius of each listed cell.
# Entry k of the outputs belongs to cell ids[k].
#
# cellV: cell-to-vertex array with 4 columns
# vx, vy, vz: coordinates of each vertex
# ids: list of cell numbers
# nIds: number of entries in ids
# ccx, ccy, ccz: output arrays for the circumcenter
# ccR: output array for the circumradius
function cell_circumcenter_3d_kernel(cellV, vx, vy, vz, ids, nIds, ccx, ccy, ccz, ccR)
    for k = 1:nIds
        c = ids[k]
        i0 = cellV[c, 1]; i1 = cellV[c, 2]; i2 = cellV[c, 3]; i3 = cellV[c, 4]
        x0 = vx[i0]
        y0 = vy[i0]
        z0 = vz[i0]
        # Move vertex 0 to the origin. The vectors b, c, d go to the other three.
        bx = vx[i1] - x0; by = vy[i1] - y0; bz = vz[i1] - z0
        cx = vx[i2] - x0; cy = vy[i2] - y0; cz = vz[i2] - z0
        dx = vx[i3] - x0; dy = vy[i3] - y0; dz = vz[i3] - z0
        b2 = bx * bx + by * by + bz * bz
        c2 = cx * cx + cy * cy + cz * cz
        d2 = dx * dx + dy * dy + dz * dz
        cdx = cy * dz - cz * dy; cdy = cz * dx - cx * dz; cdz = cx * dy - cy * dx
        dbx = dy * bz - dz * by; dby = dz * bx - dx * bz; dbz = dx * by - dy * bx
        bcx = by * cz - bz * cy; bcy = bz * cx - bx * cz; bcz = bx * cy - by * cx
        den = 2.0 * (bx * cdx + by * cdy + bz * cdz)
        ux = (b2 * cdx + c2 * dbx + d2 * bcx) / den
        uy = (b2 * cdy + c2 * dby + d2 * bcy) / den
        uz = (b2 * cdz + c2 * dbz + d2 * bcz) / den
        ccx[k] = x0 + ux
        ccy[k] = y0 + uy
        ccz[k] = z0 + uz
        ccR[k] = sqrt(ux * ux + uy * uy + uz * uz)
    end
    return nothing
end

# subfacet_circumcenter_kernel(sfV, vx, vy, vz, ids, nIds, ccx, ccy, ccz, ccR)
#
# Compute the circumcenter and the circumradius of each listed subfacet.
# The circumcenter lies in the plane of the subfacet. It is the center
# of the diametral ball.
#
# sfV: subfacet-to-vertex array with 3 columns
# vx, vy, vz: coordinates of each vertex
# ids: list of subfacet numbers
# nIds: number of entries in ids
# ccx, ccy, ccz: output arrays for the circumcenter
# ccR: output array for the circumradius
function subfacet_circumcenter_kernel(sfV, vx, vy, vz, ids, nIds, ccx, ccy, ccz, ccR)
    for k = 1:nIds
        s = ids[k]
        a = sfV[s, 1]; b = sfV[s, 2]; c = sfV[s, 3]
        # u and v are the two sides from vertex a. w is their cross product.
        ux = vx[b] - vx[a]; uy = vy[b] - vy[a]; uz = vz[b] - vz[a]
        vvx = vx[c] - vx[a]; vvy = vy[c] - vy[a]; vvz = vz[c] - vz[a]
        wx = uy * vvz - uz * vvy
        wy = uz * vvx - ux * vvz
        wz = ux * vvy - uy * vvx
        u2 = ux * ux + uy * uy + uz * uz
        v2 = vvx * vvx + vvy * vvy + vvz * vvz
        w2 = wx * wx + wy * wy + wz * wz
        # The center is a + (v2 * (w x u) + u2 * (v x w)) / (2 * w2).
        wux = wy * uz - wz * uy; wuy = wz * ux - wx * uz; wuz = wx * uy - wy * ux
        vwx = vvy * wz - vvz * wy; vwy = vvz * wx - vvx * wz; vwz = vvx * wy - vvy * wx
        px = (v2 * wux + u2 * vwx) / (2.0 * w2)
        py = (v2 * wuy + u2 * vwy) / (2.0 * w2)
        pz = (v2 * wuz + u2 * vwz) / (2.0 * w2)
        ccx[k] = vx[a] + px
        ccy[k] = vy[a] + py
        ccz[k] = vz[a] + pz
        ccR[k] = sqrt(px * px + py * py + pz * pz)
    end
    return nothing
end

# segment_encroach_3d_kernel(segV, vx, vy, vz, ids, tx, ty, tz, nIds, flag)
#
# Test whether a point encroaches on a segment. A point encroaches on a
# segment when it lies strictly inside the ball that has the segment as
# its diameter. This is the case when the two vectors from the point to
# the two segment ends have a negative dot product. The test uses a
# relative margin of 1e-9. A point on the ball, as in a regular lattice,
# does not count.
#
# segV: segment-to-vertex array with 2 columns
# vx, vy, vz: coordinates of each vertex
# ids: list of segment numbers
# tx, ty, tz: coordinates of the test point for each entry
# nIds: number of entries in ids
# flag: output array. 1 means the point encroaches on the segment.
function segment_encroach_3d_kernel(segV, vx, vy, vz, ids, tx, ty, tz, nIds, flag)
    for k = 1:nIds
        s = ids[k]
        a = segV[s, 1]
        b = segV[s, 2]
        ax = vx[a] - tx[k]; ay = vy[a] - ty[k]; az = vz[a] - tz[k]
        bx = vx[b] - tx[k]; by = vy[b] - ty[k]; bz = vz[b] - tz[k]
        lx = vx[a] - vx[b]; ly = vy[a] - vy[b]; lz = vz[a] - vz[b]
        flag[k] = 0
        if ax * bx + ay * by + az * bz < -1.0e-9 * (lx * lx + ly * ly + lz * lz)
            flag[k] = 1
        end
    end
    return nothing
end

# subfacet_encroach_kernel(sfV, vx, vy, vz, ids, tx, ty, tz, nIds, flag)
#
# Test whether a point encroaches on a subfacet. A point encroaches on a
# subfacet when it lies strictly inside the diametral ball of the
# subfacet. The test uses a relative margin of 1e-9. A point on the ball
# does not count.
#
# sfV: subfacet-to-vertex array with 3 columns
# vx, vy, vz: coordinates of each vertex
# ids: list of subfacet numbers
# tx, ty, tz: coordinates of the test point for each entry
# nIds: number of entries in ids
# flag: output array. 1 means the point encroaches on the subfacet.
function subfacet_encroach_kernel(sfV, vx, vy, vz, ids, tx, ty, tz, nIds, flag)
    for k = 1:nIds
        s = ids[k]
        a = sfV[s, 1]; b = sfV[s, 2]; c = sfV[s, 3]
        ux = vx[b] - vx[a]; uy = vy[b] - vy[a]; uz = vz[b] - vz[a]
        vvx = vx[c] - vx[a]; vvy = vy[c] - vy[a]; vvz = vz[c] - vz[a]
        wx = uy * vvz - uz * vvy
        wy = uz * vvx - ux * vvz
        wz = ux * vvy - uy * vvx
        u2 = ux * ux + uy * uy + uz * uz
        v2 = vvx * vvx + vvy * vvy + vvz * vvz
        w2 = wx * wx + wy * wy + wz * wz
        wux = wy * uz - wz * uy; wuy = wz * ux - wx * uz; wuz = wx * uy - wy * ux
        vwx = vvy * wz - vvz * wy; vwy = vvz * wx - vvx * wz; vwz = vvx * wy - vvy * wx
        px = (v2 * wux + u2 * vwx) / (2.0 * w2)
        py = (v2 * wuy + u2 * vwy) / (2.0 * w2)
        pz = (v2 * wuz + u2 * vwz) / (2.0 * w2)
        r2 = px * px + py * py + pz * pz
        dx = tx[k] - vx[a] - px
        dy = ty[k] - vy[a] - py
        dz = tz[k] - vz[a] - pz
        flag[k] = 0
        if dx * dx + dy * dy + dz * dz < r2 * (1.0 - 1.0e-9)
            flag[k] = 1
        end
    end
    return nothing
end

# segment_midpoint_kernel(segV, vx, vy, vz, ids, nIds, sx, sy, sz)
#
# Compute the midpoint of each listed segment.
#
# segV: segment-to-vertex array with 2 columns
# vx, vy, vz: coordinates of each vertex
# ids: list of segment numbers
# nIds: number of entries in ids
# sx, sy, sz: output arrays for the midpoint
function segment_midpoint_kernel(segV, vx, vy, vz, ids, nIds, sx, sy, sz)
    for k = 1:nIds
        s = ids[k]
        a = segV[s, 1]
        b = segV[s, 2]
        sx[k] = 0.5 * (vx[a] + vx[b])
        sy[k] = 0.5 * (vy[a] + vy[b])
        sz[k] = 0.5 * (vz[a] + vz[b])
    end
    return nothing
end

# ── Smoothing and quality kernels ─────────────────────────────────────

# cell_quality_2d_kernel(cellV, vx, vy, ids, nIds, ccR, hAt, ratio, sizeRatio)
#
# Compute two quality numbers for each listed cell. The first is the
# circumradius divided by the shortest edge. A small value means large
# angles. The second is the circumradius of the cell divided by the
# circumradius of an equilateral cell that has the target size.
#
# cellV: cell-to-vertex array with 3 columns
# vx: x coordinate of each vertex
# vy: y coordinate of each vertex
# ids: list of cell numbers
# nIds: number of entries in ids
# ccR: circumradius of each listed cell, from cell_circumcenter_2d_kernel
# hAt: target size at the circumcenter of each listed cell
# ratio: output array for the radius-edge ratio
# sizeRatio: output array for the size ratio. A value above 1 means
#   the cell is too large.
function cell_quality_2d_kernel(cellV, vx, vy, ids, nIds, ccR, hAt, ratio, sizeRatio)
    sqrt3 = 1.7320508075688772
    for k = 1:nIds
        c = ids[k]
        i0 = cellV[c, 1]; i1 = cellV[c, 2]; i2 = cellV[c, 3]
        ex = vx[i1] - vx[i0]; ey = vy[i1] - vy[i0]
        l01 = ex * ex + ey * ey
        ex = vx[i2] - vx[i1]; ey = vy[i2] - vy[i1]
        l12 = ex * ex + ey * ey
        ex = vx[i0] - vx[i2]; ey = vy[i0] - vy[i2]
        l20 = ex * ex + ey * ey
        shortest = sqrt(min(min(l01, l12), l20))
        ratio[k] = ccR[k] / shortest
        sizeRatio[k] = ccR[k] * sqrt3 / hAt[k]
    end
    return nothing
end

# smooth_jacobi_2d_kernel(pts, nbrOff, nbrIdx, n, avgPts)
#
# Compute the average position of the neighbors of each vertex. The
# kernel reads pts and writes avgPts. This is a Jacobi step, so no
# iteration reads a value that another iteration wrote. A vertex
# without neighbors keeps its position.
#
# pts: vertex array with x in column 1 and y in column 2
# nbrOff: length n+1 array. The neighbors of vertex v are
#   nbrIdx[nbrOff[v]] up to nbrIdx[nbrOff[v+1]-1].
# nbrIdx: list of neighbor vertex numbers
# n: number of vertices
# avgPts: n-by-2 output array
function smooth_jacobi_2d_kernel(pts, nbrOff, nbrIdx, n, avgPts)
    for v = 1:n
        k0 = nbrOff[v]
        k1 = nbrOff[v + 1]
        cnt = k1 - k0
        sx = 0.0
        sy = 0.0
        # The loop over the neighbors of one vertex carries a running sum.
        for k = k0:k1 - 1
            w = nbrIdx[k]
            sx += pts[w, 1]
            sy += pts[w, 2]
        end
        avgPts[v, 1] = pts[v, 1]
        avgPts[v, 2] = pts[v, 2]
        if cnt > 0
            avgPts[v, 1] = sx / cnt
            avgPts[v, 2] = sy / cnt
        end
    end
    return nothing
end

# blend_positions_2d_kernel(oldPts, avgPts, ptMk, n, omega, newPts)
#
# Move each interior vertex by the fraction omega of the way from its
# old position to its neighbor average. A vertex with a marker above 0
# is on the boundary and does not move.
#
# oldPts: n-by-2 array of the old positions
# avgPts: n-by-2 array of the neighbor averages
# ptMk: length-n array of vertex markers
# n: number of vertices
# omega: fraction of the move to apply
# newPts: n-by-2 output array
function blend_positions_2d_kernel(oldPts, avgPts, ptMk, n, omega, newPts)
    for v = 1:n
        x = oldPts[v, 1]
        y = oldPts[v, 2]
        if ptMk[v] == 0
            x += omega * (avgPts[v, 1] - x)
            y += omega * (avgPts[v, 2] - y)
        end
        newPts[v, 1] = x
        newPts[v, 2] = y
    end
    return nothing
end

# cell_quality_3d_kernel(cellV, vx, vy, vz, ids, nIds, ccR, hAt, ratio, sizeRatio)
#
# Compute two quality numbers for each listed cell. The first is the
# circumradius divided by the shortest edge. The second is the
# circumradius of the cell divided by the circumradius of a regular cell
# that has the target size as its edge length.
#
# cellV: cell-to-vertex array with 4 columns
# vx, vy, vz: coordinates of each vertex
# ids: list of cell numbers
# nIds: number of entries in ids
# ccR: circumradius of each listed cell, from cell_circumcenter_3d_kernel
# hAt: target size at the circumcenter of each listed cell
# ratio: output array for the radius-edge ratio
# sizeRatio: output array for the size ratio. A value above 1 means
#   the cell is too large.
function cell_quality_3d_kernel(cellV, vx, vy, vz, ids, nIds, ccR, hAt, ratio, sizeRatio)
    regular = 0.6123724356957945
    for k = 1:nIds
        c = ids[k]
        i0 = cellV[c, 1]; i1 = cellV[c, 2]; i2 = cellV[c, 3]; i3 = cellV[c, 4]
        ex = vx[i1] - vx[i0]; ey = vy[i1] - vy[i0]; ez = vz[i1] - vz[i0]
        m = ex * ex + ey * ey + ez * ez
        ex = vx[i2] - vx[i0]; ey = vy[i2] - vy[i0]; ez = vz[i2] - vz[i0]
        m = min(m, ex * ex + ey * ey + ez * ez)
        ex = vx[i3] - vx[i0]; ey = vy[i3] - vy[i0]; ez = vz[i3] - vz[i0]
        m = min(m, ex * ex + ey * ey + ez * ez)
        ex = vx[i2] - vx[i1]; ey = vy[i2] - vy[i1]; ez = vz[i2] - vz[i1]
        m = min(m, ex * ex + ey * ey + ez * ez)
        ex = vx[i3] - vx[i1]; ey = vy[i3] - vy[i1]; ez = vz[i3] - vz[i1]
        m = min(m, ex * ex + ey * ey + ez * ez)
        ex = vx[i3] - vx[i2]; ey = vy[i3] - vy[i2]; ez = vz[i3] - vz[i2]
        m = min(m, ex * ex + ey * ey + ez * ez)
        ratio[k] = ccR[k] / sqrt(m)
        sizeRatio[k] = ccR[k] / (regular * hAt[k])
    end
    return nothing
end

# tet_dihedral_kernel(cells, pts, n, minD, maxD)
#
# Compute the smallest and the largest dihedral angle of each cell, in
# radians. The cells must have the standard vertex order. The outward
# normals of the four faces give the angles.
#
# cells: n-by-4 cell-to-vertex array
# pts: vertex array with x, y, z in columns 1 to 3
# n: number of cells
# minD: length-n output array for the smallest angle
# maxD: length-n output array for the largest angle
function tet_dihedral_kernel(cells, pts, n, minD, maxD)
    for c = 1:n
        i0 = cells[c, 1]; i1 = cells[c, 2]; i2 = cells[c, 3]; i3 = cells[c, 4]
        ax = pts[i1, 1] - pts[i0, 1]; ay = pts[i1, 2] - pts[i0, 2]; az = pts[i1, 3] - pts[i0, 3]
        bx = pts[i2, 1] - pts[i0, 1]; by = pts[i2, 2] - pts[i0, 2]; bz = pts[i2, 3] - pts[i0, 3]
        cx = pts[i3, 1] - pts[i0, 1]; cy = pts[i3, 2] - pts[i0, 2]; cz = pts[i3, 3] - pts[i0, 3]
        # Outward normals, scaled by twice the face area.
        # Face 1 is opposite vertex 0 and so on.
        n1x = (by - ay) * (cz - az) - (bz - az) * (cy - ay)
        n1y = (bz - az) * (cx - ax) - (bx - ax) * (cz - az)
        n1z = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
        n2x = cy * bz - cz * by
        n2y = cz * bx - cx * bz
        n2z = cx * by - cy * bx
        n3x = ay * cz - az * cy
        n3y = az * cx - ax * cz
        n3z = ax * cy - ay * cx
        n4x = by * az - bz * ay
        n4y = bz * ax - bx * az
        n4z = bx * ay - by * ax
        l1 = sqrt(n1x * n1x + n1y * n1y + n1z * n1z)
        l2 = sqrt(n2x * n2x + n2y * n2y + n2z * n2z)
        l3 = sqrt(n3x * n3x + n3y * n3y + n3z * n3z)
        l4 = sqrt(n4x * n4x + n4y * n4y + n4z * n4z)
        # cos of the dihedral angle between two faces is minus the dot
        # product of their unit outward normals.
        d12 = -(n1x * n2x + n1y * n2y + n1z * n2z) / (l1 * l2)
        d13 = -(n1x * n3x + n1y * n3y + n1z * n3z) / (l1 * l3)
        d14 = -(n1x * n4x + n1y * n4y + n1z * n4z) / (l1 * l4)
        d23 = -(n2x * n3x + n2y * n3y + n2z * n3z) / (l2 * l3)
        d24 = -(n2x * n4x + n2y * n4y + n2z * n4z) / (l2 * l4)
        d34 = -(n3x * n4x + n3y * n4y + n3z * n4z) / (l3 * l4)
        cmax = max(max(max(d12, d13), max(d14, d23)), max(d24, d34))
        cmin = min(min(min(d12, d13), min(d14, d23)), min(d24, d34))
        minD[c] = acos(min(max(cmax, -1.0), 1.0))
        maxD[c] = acos(min(max(cmin, -1.0), 1.0))
    end
    return nothing
end

# group_extrema_kernel(vals, idx, start, count, nG, outMin, outMax, outSum)
#
# Reduce each group of numbers to its smallest value, its largest value
# and its sum. Group g holds the numbers vals[idx[k]], with k from
# start[g] up to start[g] + count[g] - 1. Each group has at least one
# member. The groups are independent of each other.
#
# vals: array of numbers
# idx: list of positions in vals, group after group
# start: length-nG array, position in idx of the first member of each group
# count: length-nG array, number of members of each group
# nG: number of groups
# outMin: length-nG output array for the smallest value of each group
# outMax: length-nG output array for the largest value of each group
# outSum: length-nG output array for the sum of each group
function group_extrema_kernel(vals, idx, start, count, nG, outMin, outMax, outSum)
    for g = 1:nG
        k0 = start[g]
        j0 = idx[k0]
        lo = vals[j0]
        hi = vals[j0]
        total = 0.0
        # running reduction over the members of one group, carried
        # between the inner iterations
        for k = k0:k0 + count[g] - 1
            j = idx[k]
            v = vals[j]
            lo = min(lo, v)
            hi = max(hi, v)
            total += v
        end
        outMin[g] = lo
        outMax[g] = hi
        outSum[g] = total
    end
    return nothing
end

# edge_midpoint_kernel(pts, edges, nV, nE, nd, out)
#
# Write the midpoint of every edge into out. The midpoint of edge e goes
# to row nV + e. The rows 1 to nV of out stay as the caller set them.
#
# pts: vertex array with nd columns
# edges: nE-by-2 edge-to-vertex array, 1-based
# nV: number of vertices in pts
# nE: number of edges
# nd: number of columns of pts, 2 or 3
# out: (nV + nE)-by-nd output array
function edge_midpoint_kernel(pts, edges, nV, nE, nd, out)
    # One fused loop over the edges and the columns.
    for idx = 1:nE * nd
        e = div(idx - 1, nd) + 1
        d = mod(idx - 1, nd) + 1
        a = edges[e, 1]
        b = edges[e, 2]
        out[nV + e, d] = 0.5 * (pts[a, d] + pts[b, d])
    end
    return nothing
end

# refine_triangles_2d_kernel(cells, c2e, nV, nCells, out)
#
# Split each triangle into 4 triangles. Three corner triangles and one
# center triangle make the children. The children of cell k are the rows
# 4k-3 to 4k of out. A midpoint has the vertex number nV + 1 + c2e[k, j],
# where c2e holds 0-based edge numbers. The local edges 1, 2 and 3 join
# the local vertices (1, 2), (1, 3) and (2, 3).
#
# cells: nCells-by-3 cell-to-vertex array, 1-based
# c2e: nCells-by-3 cell-to-edge array, 0-based
# nV: number of vertices of the coarse mesh
# nCells: number of coarse cells
# out: 4 nCells-by-3 output array
function refine_triangles_2d_kernel(cells, c2e, nV, nCells, out)
    for k = 1:nCells
        a = cells[k, 1]
        b = cells[k, 2]
        c = cells[k, 3]
        mab = nV + 1 + c2e[k, 1]
        mac = nV + 1 + c2e[k, 2]
        mbc = nV + 1 + c2e[k, 3]
        r = 4 * (k - 1)
        # corner at a
        out[r + 1, 1] = a
        out[r + 1, 2] = mab
        out[r + 1, 3] = mac
        # corner at b
        out[r + 2, 1] = mab
        out[r + 2, 2] = b
        out[r + 2, 3] = mbc
        # corner at c
        out[r + 3, 1] = mac
        out[r + 3, 2] = mbc
        out[r + 3, 3] = c
        # center
        out[r + 4, 1] = mab
        out[r + 4, 2] = mbc
        out[r + 4, 3] = mac
    end
    return nothing
end

# tet_shape_kernel(cells, pts, n, out)
#
# Write the cube of the mean-ratio shape measure of each tetrahedron.
# The measure is 12 (3 V)^(2/3) divided by the sum of the squared edge
# lengths. It is 1 for a regular tetrahedron and goes to 0 for a sliver.
# The cube is 432 D^2 / S^3. D is the determinant of the 3 edge vectors
# from vertex 1. S is the sum of the 6 squared edge lengths. The cube
# keeps the order of the measure and needs no fractional power. The sign
# of D does not matter.
#
# cells: n-by-4 cell-to-vertex array, 1-based
# pts: vertex array with x, y, z in columns 1 to 3
# n: number of cells
# out: length-n output array
function tet_shape_kernel(cells, pts, n, out)
    for c = 1:n
        a = cells[c, 1]
        b = cells[c, 2]
        d = cells[c, 3]
        e = cells[c, 4]
        # edge vectors from vertex 1
        ux = pts[b, 1] - pts[a, 1]
        uy = pts[b, 2] - pts[a, 2]
        uz = pts[b, 3] - pts[a, 3]
        vx = pts[d, 1] - pts[a, 1]
        vy = pts[d, 2] - pts[a, 2]
        vz = pts[d, 3] - pts[a, 3]
        wx = pts[e, 1] - pts[a, 1]
        wy = pts[e, 2] - pts[a, 2]
        wz = pts[e, 3] - pts[a, 3]
        det = ux * (vy * wz - vz * wy) - uy * (vx * wz - vz * wx) + uz * (vx * wy - vy * wx)
        # sum of the 6 squared edge lengths
        s = ux * ux + uy * uy + uz * uz
        s += vx * vx + vy * vy + vz * vz
        s += wx * wx + wy * wy + wz * wz
        s += (vx - ux) * (vx - ux) + (vy - uy) * (vy - uy) + (vz - uz) * (vz - uz)
        s += (wx - ux) * (wx - ux) + (wy - uy) * (wy - uy) + (wz - uz) * (wz - uz)
        s += (wx - vx) * (wx - vx) + (wy - vy) * (wy - vy) + (wz - vz) * (wz - vz)
        out[c] = 432.0 * det * det / (s * s * s)
    end
    return nothing
end

# refine_tets_3d_kernel(cells, c2e, bey, dgSel, nV, nCells, out)
#
# Split each tetrahedron into 8 tetrahedra (the method of Bey). Four
# corner children keep one vertex of the parent each. An octahedron stays
# between them. The diagonal dgSel[k] cuts it into 4 children. The
# children of cell k are the rows 8k-7 to 8k of out. The table bey holds
# a code for each vertex of each child. The codes 1 to 4 are the vertices
# of the parent. The codes 5 to 10 are the midpoints of the local edges 1
# to 6. A midpoint has the vertex number nV + 1 + c2e[k, code - 4].
# The orientation of a child is not fixed here.
#
# cells: nCells-by-4 cell-to-vertex array, 1-based
# c2e: nCells-by-6 cell-to-edge array, 0-based
# bey: 3-by-8-by-4 table of codes (diagonal, child, vertex)
# dgSel: diagonal (1, 2 or 3) of each cell
# nV: number of vertices of the coarse mesh
# nCells: number of coarse cells
# out: 8 nCells-by-4 output array
function refine_tets_3d_kernel(cells, c2e, bey, dgSel, nV, nCells, out)
    # One fused loop over the cells, the 8 children and the 4 vertices.
    for idx = 1:nCells * 32
        k = div(idx - 1, 32) + 1
        r = mod(idx - 1, 32)
        j = div(r, 4) + 1
        i = mod(r, 4) + 1
        dg = dgSel[k]
        code = bey[dg, j, i]
        # 0 for a vertex code, 1 for a midpoint code
        isMid = min(div(code - 1, 4), 1)
        cv = min(code, 4)
        cm = max(code - 4, 1)
        va = cells[k, cv]
        me = c2e[k, cm]
        out[8 * (k - 1) + j, i] = (1 - isMid) * va + isMid * (nV + 1 + me)
    end
    return nothing
end

# diagonal_choice_kernel(eta1, eta2, eta3, nCells, dgSel)
#
# Choose the diagonal of each octahedron by the max-min rule. For each
# diagonal, take the smallest shape value of the 4 interior children.
# Keep the diagonal with the largest of these values. A tie goes to the
# lowest diagonal number.
#
# eta1: shape values of the 8 children of each cell, for diagonal 1
# eta2: the same for diagonal 2
# eta3: the same for diagonal 3
# nCells: number of coarse cells
# dgSel: length-nCells output array with the chosen diagonal
function diagonal_choice_kernel(eta1, eta2, eta3, nCells, dgSel)
    for k = 1:nCells
        r = 8 * (k - 1)
        # the children 5 to 8 are the interior children
        w1 = min(min(eta1[r + 5], eta1[r + 6]), min(eta1[r + 7], eta1[r + 8]))
        w2 = min(min(eta2[r + 5], eta2[r + 6]), min(eta2[r + 7], eta2[r + 8]))
        w3 = min(min(eta3[r + 5], eta3[r + 6]), min(eta3[r + 7], eta3[r + 8]))
        best = 1
        wBest = w1
        if w2 > wBest
            best = 2
            wBest = w2
        end
        if w3 > wBest
            best = 3
        end
        dgSel[k] = best
    end
    return nothing
end

# smooth_jacobi_3d_kernel(pts, nbrOff, nbrIdx, n, avgPts)
#
# Compute the average position of the neighbors of each vertex. The
# kernel reads pts and writes avgPts. This is a Jacobi step, so no
# iteration reads a value that another iteration wrote. A vertex
# without neighbors keeps its position.
#
# pts: vertex array with x, y, z in columns 1 to 3
# nbrOff: length n+1 array. The neighbors of vertex v are
#   nbrIdx[nbrOff[v]] up to nbrIdx[nbrOff[v+1]-1].
# nbrIdx: list of neighbor vertex numbers
# n: number of vertices
# avgPts: n-by-3 output array
function smooth_jacobi_3d_kernel(pts, nbrOff, nbrIdx, n, avgPts)
    for v = 1:n
        k0 = nbrOff[v]
        k1 = nbrOff[v + 1]
        cnt = k1 - k0
        sx = 0.0
        sy = 0.0
        sz = 0.0
        # The loop over the neighbors of one vertex carries a running sum.
        for k = k0:k1 - 1
            w = nbrIdx[k]
            sx += pts[w, 1]
            sy += pts[w, 2]
            sz += pts[w, 3]
        end
        avgPts[v, 1] = pts[v, 1]
        avgPts[v, 2] = pts[v, 2]
        avgPts[v, 3] = pts[v, 3]
        if cnt > 0
            avgPts[v, 1] = sx / cnt
            avgPts[v, 2] = sy / cnt
            avgPts[v, 3] = sz / cnt
        end
    end
    return nothing
end

# smooth_accept_3d_kernel(oldPts, avgPts, ptMk, v2cOff, v2cIdx, cells, n, omega, newPts)
#
# Move each interior vertex toward its neighbor average, but only when
# the move does not make its cells worse. The kernel moves the vertex by
# the fraction omega of the way. It then measures the shape of each cell
# around the vertex before and after the move. The shape measure is
# 432 * W * abs(W) / S^3. W is six times the volume. S is the sum of
# the six squared edge lengths. The measure is 1 for a regular cell and
# near 0 for a flat cell. It is negative for an inverted cell. The vertex
# keeps the move when the worst measure of its cells does not fall. The
# worst measure must also stay positive. A vertex with a marker above 0 is on the boundary and does
# not move. The neighbors use their old positions.
#
# oldPts: n-by-3 array of the old positions
# avgPts: n-by-3 array of the neighbor averages
# ptMk: length-n array of vertex markers
# v2cOff: length n+1 array. The cells of vertex v are
#   v2cIdx[v2cOff[v] + 1] up to v2cIdx[v2cOff[v + 1]]. Each stored
#   cell number is 0-based.
# v2cIdx: list of cell numbers
# cells: cell-to-vertex array with 4 columns
# n: number of vertices
# omega: fraction of the move to try
# newPts: n-by-3 output array
function smooth_accept_3d_kernel(oldPts, avgPts, ptMk, v2cOff, v2cIdx, cells, n, omega, newPts)
    for v = 1:n
        x = oldPts[v, 1]
        y = oldPts[v, 2]
        z = oldPts[v, 3]
        qx = x
        qy = y
        qz = z
        if ptMk[v] == 0
            qx = x + omega * (avgPts[v, 1] - x)
            qy = y + omega * (avgPts[v, 2] - y)
            qz = z + omega * (avgPts[v, 3] - z)
            qOld = 1.0e30
            qNew = 1.0e30
            # The loop over the cells of one vertex carries two running minima.
            for k = v2cOff[v] + 1:v2cOff[v + 1]
                c = v2cIdx[k] + 1
                a0 = cells[c, 1]; a1 = cells[c, 2]; a2 = cells[c, 3]; a3 = cells[c, 4]
                for pass = 1:2
                    sx = x
                    sy = y
                    sz = z
                    if pass == 2
                        sx = qx
                        sy = qy
                        sz = qz
                    end
                    x0 = oldPts[a0, 1]; y0 = oldPts[a0, 2]; z0 = oldPts[a0, 3]
                    x1 = oldPts[a1, 1]; y1 = oldPts[a1, 2]; z1 = oldPts[a1, 3]
                    x2 = oldPts[a2, 1]; y2 = oldPts[a2, 2]; z2 = oldPts[a2, 3]
                    x3 = oldPts[a3, 1]; y3 = oldPts[a3, 2]; z3 = oldPts[a3, 3]
                    if a0 == v
                        x0 = sx; y0 = sy; z0 = sz
                    end
                    if a1 == v
                        x1 = sx; y1 = sy; z1 = sz
                    end
                    if a2 == v
                        x2 = sx; y2 = sy; z2 = sz
                    end
                    if a3 == v
                        x3 = sx; y3 = sy; z3 = sz
                    end
                    ax = x1 - x0; ay = y1 - y0; az = z1 - z0
                    bx = x2 - x0; by = y2 - y0; bz = z2 - z0
                    cx = x3 - x0; cy = y3 - y0; cz = z3 - z0
                    vol6 = ax * (by * cz - bz * cy) + ay * (bz * cx - bx * cz) + az * (bx * cy - by * cx)
                    s = ax * ax + ay * ay + az * az + bx * bx + by * by + bz * bz + cx * cx + cy * cy + cz * cz
                    s = s + (bx - ax) * (bx - ax) + (by - ay) * (by - ay) + (bz - az) * (bz - az)
                    s = s + (cx - ax) * (cx - ax) + (cy - ay) * (cy - ay) + (cz - az) * (cz - az)
                    s = s + (cx - bx) * (cx - bx) + (cy - by) * (cy - by) + (cz - bz) * (cz - bz)
                    quality = 432.0 * vol6 * abs(vol6) / (s * s * s)
                    if pass == 1
                        qOld = min(qOld, quality)
                    else
                        qNew = min(qNew, quality)
                    end
                end
            end
            if qNew < qOld
                qx = x
                qy = y
                qz = z
            end
            if qNew <= 0.0
                qx = x
                qy = y
                qz = z
            end
        end
        newPts[v, 1] = qx
        newPts[v, 2] = qy
        newPts[v, 3] = qz
    end
    return nothing
end

# =====================================================================
# end of computing kernels
# =====================================================================

# =====================================================================
# shared driver code -- no package needed
# =====================================================================

# ── Exact predicates ─────────────────────────────────────────────────────
#
# The 2-D mesher needs two predicates: orient2d and incircle. Each one
# returns the exact sign as an Int. A plain Float64 test gives wrong
# signs on the collinear and cocircular vertices of a structured
# boundary. Each predicate has two stages. The filter stage computes the
# value in Float64 and compares it with an error bound. The exact stage
# computes the value with floating-point expansions.
#
# An expansion is a vector of Float64 values. The values do not overlap.
# They increase in magnitude. Their exact sum is the number that the
# expansion stands for. The sign of the number is the sign of the last
# value.

# The unit round-off of Float64 with round-to-nearest, 2^-53. The error
# bounds of the filtered predicates are small functions of it. The
# compiler folds each call to a number.
@inline _eps53() = 1.1102230246251565e-16
@inline _o2d_errbound() = (3.0 + 16.0 * _eps53()) * _eps53()
@inline _icc_errbound() = (10.0 + 96.0 * _eps53()) * _eps53()

# Return x and y such that x + y equals a + b exactly. x is the rounded sum.
@inline function _two_sum(a::Float64, b::Float64)
    x = a + b
    bv = x - a
    av = x - bv
    return x, (a - av) + (b - bv)
end

# Return x and y such that x + y equals a * b exactly. x is the rounded product.
@inline function _two_prod(a::Float64, b::Float64)
    x = a * b
    return x, fma(a, b, -x)
end

# Add one Float64 to an expansion. The result is an expansion.
function _exp_grow(e::Vector{Float64}, b::Float64)
    h = Float64[]
    q = b
    for ei in e
        q, hh = _two_sum(q, ei)
        hh != 0.0 && push!(h, hh)
    end
    (q != 0.0 || isempty(h)) && push!(h, q)
    return h
end

# Add two expansions.
function _exp_sum(e::Vector{Float64}, f::Vector{Float64})
    h = e
    for fi in f
        h = _exp_grow(h, fi)
    end
    return h
end

# Multiply an expansion by one Float64.
function _exp_scale(e::Vector{Float64}, b::Float64)
    h = Float64[]
    isempty(e) && return [0.0]
    q, hh = _two_prod(e[1], b)
    hh != 0.0 && push!(h, hh)
    for i in 2:length(e)
        t1, t0 = _two_prod(e[i], b)
        s, err = _two_sum(q, t0)
        err != 0.0 && push!(h, err)
        q, err2 = _two_sum(t1, s)
        err2 != 0.0 && push!(h, err2)
    end
    (q != 0.0 || isempty(h)) && push!(h, q)
    return h
end

# Multiply two expansions.
function _exp_mul(e::Vector{Float64}, f::Vector{Float64})
    h = Float64[0.0]
    for fj in f
        h = _exp_sum(h, _exp_scale(e, fj))
    end
    return h
end

# Return the expansion of a - b, exact.
function _exp_diff(a::Float64, b::Float64)
    x, y = _two_sum(a, -b)
    return y != 0.0 ? [y, x] : [x]
end

# Return the sign of an expansion as -1, 0, or 1.
function _exp_sign(e::Vector{Float64})
    s = e[end]
    return s > 0.0 ? 1 : (s < 0.0 ? -1 : 0)
end

# Exact stage of orient2d.
function _orient2d_exact(ax, ay, bx, by, cx, cy)
    acx = _exp_diff(ax, cx)
    bcy = _exp_diff(by, cy)
    acy = _exp_diff(ay, cy)
    bcx = _exp_diff(bx, cx)
    left = _exp_mul(acx, bcy)
    right = _exp_mul(acy, bcx)
    return _exp_sign(_exp_sum(left, -right))
end

# Return 1 if a, b, c turn counterclockwise, -1 if clockwise, and 0 if
# the three points are collinear.
function _orient2d(ax::Float64, ay::Float64, bx::Float64, by::Float64, cx::Float64, cy::Float64)
    detleft = (ax - cx) * (by - cy)
    detright = (ay - cy) * (bx - cx)
    # Both products are zero when one factor of each is exactly zero.
    # Then the determinant is exactly zero.
    if detleft == 0.0 && detright == 0.0
        return 0
    end
    det = detleft - detright
    errbound = _o2d_errbound() * (abs(detleft) + abs(detright))
    det > errbound && return 1
    -det > errbound && return -1
    return _orient2d_exact(ax, ay, bx, by, cx, cy)
end

# Exact stage of incircle.
function _incircle_exact(ax, ay, bx, by, cx, cy, dx, dy)
    adx = _exp_diff(ax, dx); ady = _exp_diff(ay, dy)
    bdx = _exp_diff(bx, dx); bdy = _exp_diff(by, dy)
    cdx = _exp_diff(cx, dx); cdy = _exp_diff(cy, dy)
    alift = _exp_sum(_exp_mul(adx, adx), _exp_mul(ady, ady))
    blift = _exp_sum(_exp_mul(bdx, bdx), _exp_mul(bdy, bdy))
    clift = _exp_sum(_exp_mul(cdx, cdx), _exp_mul(cdy, cdy))
    t1 = _exp_sum(_exp_mul(bdx, cdy), -_exp_mul(cdx, bdy))
    t2 = _exp_sum(_exp_mul(cdx, ady), -_exp_mul(adx, cdy))
    t3 = _exp_sum(_exp_mul(adx, bdy), -_exp_mul(bdx, ady))
    det = _exp_sum(_exp_sum(_exp_mul(alift, t1), _exp_mul(blift, t2)), _exp_mul(clift, t3))
    return _exp_sign(det)
end

# Return 1 if d lies strictly inside the circle through a, b, c. Return
# -1 if d lies strictly outside. Return 0 if d lies on the circle. The
# points a, b, c must turn counterclockwise.
function _incircle(ax::Float64, ay::Float64, bx::Float64, by::Float64, cx::Float64, cy::Float64, dx::Float64, dy::Float64)
    adx = ax - dx; ady = ay - dy
    bdx = bx - dx; bdy = by - dy
    cdx = cx - dx; cdy = cy - dy
    bdxcdy = bdx * cdy; cdxbdy = cdx * bdy
    cdxady = cdx * ady; adxcdy = adx * cdy
    adxbdy = adx * bdy; bdxady = bdx * ady
    alift = adx * adx + ady * ady
    blift = bdx * bdx + bdy * bdy
    clift = cdx * cdx + cdy * cdy
    det = alift * (bdxcdy - cdxbdy) + blift * (cdxady - adxcdy) + clift * (adxbdy - bdxady)
    permanent = (abs(bdxcdy) + abs(cdxbdy)) * alift +
                (abs(cdxady) + abs(adxcdy)) * blift +
                (abs(adxbdy) + abs(bdxady)) * clift
    errbound = _icc_errbound() * permanent
    det > errbound && return 1
    -det > errbound && return -1
    return _incircle_exact(ax, ay, bx, by, cx, cy, dx, dy)
end

# ── 3-D predicates ───────────────────────────────────────────────────────
#
# The 3-D mesher needs two more predicates: orient3d and insphere. Each
# one has the same two stages as the 2-D predicates. The sign convention
# matches the standard vertex order of a cell. The function _orient3d
# returns 1 when the determinant of (b - a, c - a, d - a) is positive.
# This is the sign of the signed volume of the cell (a, b, c, d). The
# function _insphere returns 1 when e lies strictly inside the sphere
# through a, b, c, d. The four points a, b, c, d must have a positive
# _orient3d.

@inline _o3d_errbound() = (7.0 + 56.0 * _eps53()) * _eps53()
@inline _isp_errbound() = (16.0 + 224.0 * _eps53()) * _eps53()

# Exact stage of orient3d. The function returns the sign of the
# determinant of (a - d, b - d, c - d). This is the opposite sign of the
# standard orientation.
function _orient3d_exact(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz)
    adx = _exp_diff(ax, dx); ady = _exp_diff(ay, dy); adz = _exp_diff(az, dz)
    bdx = _exp_diff(bx, dx); bdy = _exp_diff(by, dy); bdz = _exp_diff(bz, dz)
    cdx = _exp_diff(cx, dx); cdy = _exp_diff(cy, dy); cdz = _exp_diff(cz, dz)
    t1 = _exp_sum(_exp_mul(bdx, cdy), -_exp_mul(cdx, bdy))
    t2 = _exp_sum(_exp_mul(cdx, ady), -_exp_mul(adx, cdy))
    t3 = _exp_sum(_exp_mul(adx, bdy), -_exp_mul(bdx, ady))
    det = _exp_sum(_exp_sum(_exp_mul(adz, t1), _exp_mul(bdz, t2)), _exp_mul(cdz, t3))
    return _exp_sign(det)
end

# Return 1 when the cell (a, b, c, d) has a positive volume. Return -1
# when it has a negative volume. Return 0 when the four points lie in
# one plane.
function _orient3d(ax::Float64, ay::Float64, az::Float64, bx::Float64, by::Float64, bz::Float64,
                   cx::Float64, cy::Float64, cz::Float64, dx::Float64, dy::Float64, dz::Float64)
    adx = ax - dx; ady = ay - dy; adz = az - dz
    bdx = bx - dx; bdy = by - dy; bdz = bz - dz
    cdx = cx - dx; cdy = cy - dy; cdz = cz - dz
    bdxcdy = bdx * cdy; cdxbdy = cdx * bdy
    cdxady = cdx * ady; adxcdy = adx * cdy
    adxbdy = adx * bdy; bdxady = bdx * ady
    det = adz * (bdxcdy - cdxbdy) + bdz * (cdxady - adxcdy) + cdz * (adxbdy - bdxady)
    permanent = (abs(bdxcdy) + abs(cdxbdy)) * abs(adz) +
                (abs(cdxady) + abs(adxcdy)) * abs(bdz) +
                (abs(adxbdy) + abs(bdxady)) * abs(cdz)
    errbound = _o3d_errbound() * permanent
    det > errbound && return -1
    -det > errbound && return 1
    return -_orient3d_exact(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz)
end

# Exact stage of insphere. The function returns the sign of the
# determinant of Shewchuk. This sign is positive when e lies inside the
# sphere and (a, b, c, d) has a negative standard orientation.
function _insphere_exact(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, ex, ey, ez)
    aex = _exp_diff(ax, ex); aey = _exp_diff(ay, ey); aez = _exp_diff(az, ez)
    bex = _exp_diff(bx, ex); bey = _exp_diff(by, ey); bez = _exp_diff(bz, ez)
    cex = _exp_diff(cx, ex); cey = _exp_diff(cy, ey); cez = _exp_diff(cz, ez)
    dex = _exp_diff(dx, ex); dey = _exp_diff(dy, ey); dez = _exp_diff(dz, ez)
    ab = _exp_sum(_exp_mul(aex, bey), -_exp_mul(bex, aey))
    bc = _exp_sum(_exp_mul(bex, cey), -_exp_mul(cex, bey))
    cd = _exp_sum(_exp_mul(cex, dey), -_exp_mul(dex, cey))
    da = _exp_sum(_exp_mul(dex, aey), -_exp_mul(aex, dey))
    ac = _exp_sum(_exp_mul(aex, cey), -_exp_mul(cex, aey))
    bd = _exp_sum(_exp_mul(bex, dey), -_exp_mul(dex, bey))
    abc = _exp_sum(_exp_sum(_exp_mul(aez, bc), -_exp_mul(bez, ac)), _exp_mul(cez, ab))
    bcd = _exp_sum(_exp_sum(_exp_mul(bez, cd), -_exp_mul(cez, bd)), _exp_mul(dez, bc))
    cda = _exp_sum(_exp_sum(_exp_mul(cez, da), _exp_mul(dez, ac)), _exp_mul(aez, cd))
    dab = _exp_sum(_exp_sum(_exp_mul(dez, ab), _exp_mul(aez, bd)), _exp_mul(bez, da))
    alift = _exp_sum(_exp_sum(_exp_mul(aex, aex), _exp_mul(aey, aey)), _exp_mul(aez, aez))
    blift = _exp_sum(_exp_sum(_exp_mul(bex, bex), _exp_mul(bey, bey)), _exp_mul(bez, bez))
    clift = _exp_sum(_exp_sum(_exp_mul(cex, cex), _exp_mul(cey, cey)), _exp_mul(cez, cez))
    dlift = _exp_sum(_exp_sum(_exp_mul(dex, dex), _exp_mul(dey, dey)), _exp_mul(dez, dez))
    left = _exp_sum(_exp_mul(dlift, abc), -_exp_mul(clift, dab))
    right = _exp_sum(_exp_mul(blift, cda), -_exp_mul(alift, bcd))
    return _exp_sign(_exp_sum(left, right))
end

# Floating-point stage of the insphere test. Return 1 when e lies
# strictly inside the sphere through a, b, c, d. Return -1 when e lies
# strictly outside. Return 2 when the rounding error bound cannot decide.
function _insphere_fast(ax::Float64, ay::Float64, az::Float64, bx::Float64, by::Float64, bz::Float64,
                   cx::Float64, cy::Float64, cz::Float64, dx::Float64, dy::Float64, dz::Float64,
                   ex::Float64, ey::Float64, ez::Float64)
    aex = ax - ex; bex = bx - ex; cex = cx - ex; dex = dx - ex
    aey = ay - ey; bey = by - ey; cey = cy - ey; dey = dy - ey
    aez = az - ez; bez = bz - ez; cez = cz - ez; dez = dz - ez
    aexbey = aex * bey; bexaey = bex * aey; ab = aexbey - bexaey
    bexcey = bex * cey; cexbey = cex * bey; bc = bexcey - cexbey
    cexdey = cex * dey; dexcey = dex * cey; cd = cexdey - dexcey
    dexaey = dex * aey; aexdey = aex * dey; da = dexaey - aexdey
    aexcey = aex * cey; cexaey = cex * aey; ac = aexcey - cexaey
    bexdey = bex * dey; dexbey = dex * bey; bd = bexdey - dexbey
    abc = aez * bc - bez * ac + cez * ab
    bcd = bez * cd - cez * bd + dez * bc
    cda = cez * da + dez * ac + aez * cd
    dab = dez * ab + aez * bd + bez * da
    alift = aex * aex + aey * aey + aez * aez
    blift = bex * bex + bey * bey + bez * bez
    clift = cex * cex + cey * cey + cez * cez
    dlift = dex * dex + dey * dey + dez * dez
    det = (dlift * abc - clift * dab) + (blift * cda - alift * bcd)
    permanent = ((abs(cexdey) + abs(dexcey)) * abs(bez) +
                 (abs(dexbey) + abs(bexdey)) * abs(cez) +
                 (abs(bexcey) + abs(cexbey)) * abs(dez)) * alift +
                ((abs(dexaey) + abs(aexdey)) * abs(cez) +
                 (abs(aexcey) + abs(cexaey)) * abs(dez) +
                 (abs(cexdey) + abs(dexcey)) * abs(aez)) * blift +
                ((abs(aexbey) + abs(bexaey)) * abs(dez) +
                 (abs(bexdey) + abs(dexbey)) * abs(aez) +
                 (abs(dexaey) + abs(aexdey)) * abs(bez)) * clift +
                ((abs(bexcey) + abs(cexbey)) * abs(aez) +
                 (abs(cexaey) + abs(aexcey)) * abs(bez) +
                 (abs(aexbey) + abs(bexaey)) * abs(cez)) * dlift
    errbound = _isp_errbound() * permanent
    det > errbound && return -1
    -det > errbound && return 1
    return 2
end

# Return 1 when e lies strictly inside the sphere through a, b, c, d.
# Return -1 when e lies strictly outside. Return 0 when e lies on the
# sphere. The cell (a, b, c, d) must have a positive volume.
function _insphere(ax::Float64, ay::Float64, az::Float64, bx::Float64, by::Float64, bz::Float64,
                   cx::Float64, cy::Float64, cz::Float64, dx::Float64, dy::Float64, dz::Float64,
                   ex::Float64, ey::Float64, ez::Float64)
    s = _insphere_fast(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, ex, ey, ez)
    s != 2 && return s
    return -_insphere_exact(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, ex, ey, ez)
end

# ── Queue ────────────────────────────────────────────────────────────────

# A first-in first-out queue of Int values. The queue is a NamedTuple.
# The field items holds the values. The field head is a Ref that holds
# the index of the front value.
function _new_queue()
    return (items = Int[], head = Ref(1))
end
_qpush!(q::NamedTuple, x::Int) = push!(q.items, x)
_qempty(q::NamedTuple) = q.head[] > length(q.items)
function _qpop!(q::NamedTuple)
    x = q.items[q.head[]]
    q.head[] += 1
    # Remove the consumed front part when it becomes large.
    if q.head[] > 100000 && 2 * q.head[] > length(q.items)
        deleteat!(q.items, 1:q.head[] - 1)
        q.head[] = 1
    end
    return x
end

# Return a copy of A with newRows rows. The first rows hold the old values.
function _grow_matrix(A::Matrix{Int}, newRows::Int)
    B = zeros(Int, newRows, size(A, 2))
    B[1:size(A, 1), :] = A
    return B
end

# ── Edges, cell-to-edge, vertex-to-cell ──────────────────────────────────
#
# These functions are ordinary driver code, not kernels. They sort and
# count integer keys. This is not the fixed-shape loop style of
# skill-stade. The edge numbers that they store are 0-based. The cell
# numbers and vertex numbers are 1-based. The order of the edges is the
# order in which they first appear when the loop goes over the cells,
# then over the local edges. One loop serves a triangle (3 vertices per
# cell) and a tetrahedron (4 vertices per cell).

# Build the edge-to-vertex array and the cell-to-edge array of a mesh of
# triangles (dim 2) or of tetrahedra (dim 3). A sort of the edge keys
# replaces a Dict. The sort is stable, so the first member of each group
# is the first appearance. The vertex numbers must be less than 2^32.
function _build_edge_registry(cells::Matrix{Int}, dim::Int)
    localEdges = dim == 2 ? _tri_local_edges() : _tet_local_edges()
    nLocal = length(localEdges)
    nC = size(cells, 1)
    m = nLocal * nC
    edgeKey = Vector{Int}(undef, m)
    for c in 1:nC
        for (li, (a, b)) in enumerate(localEdges)
            va = cells[c, a]
            vb = cells[c, b]
            edgeKey[nLocal * (c - 1) + li] = va < vb ? (va << 32) | vb : (vb << 32) | va
        end
    end
    perm = sortperm(edgeKey; alg = MergeSort)
    groupOf = Vector{Int}(undef, m)
    firstOcc = Int[]
    repKey = Int[]
    g = 0
    for j in 1:m
        if j == 1 || edgeKey[perm[j]] != edgeKey[perm[j - 1]]
            g += 1
            push!(firstOcc, perm[j])
            push!(repKey, edgeKey[perm[j]])
        end
        groupOf[perm[j]] = g
    end
    order = sortperm(firstOcc)
    rank = Vector{Int}(undef, g)
    for (r, gi) in enumerate(order)
        rank[gi] = r - 1
    end
    edges = Matrix{Int}(undef, g, 2)
    for gi in 1:g
        r = rank[gi] + 1
        edges[r, 1] = repKey[gi] >> 32
        edges[r, 2] = repKey[gi] & 0xffffffff
    end
    c2e = Matrix{Int}(undef, nC, nLocal)
    for c in 1:nC, li in 1:nLocal
        c2e[c, li] = rank[groupOf[nLocal * (c - 1) + li]]
    end
    return edges, c2e
end

# Build the vertex-to-cell lists in compressed form. The function counts
# the cells of each vertex, then fills one flat array. The stored cell
# numbers are 0-based, to match the HDF5 file.
function _vertex_to_cell(nPts::Int, cells::Matrix{Int})
    counts = zeros(Int, nPts)
    for c in 1:size(cells, 1), k in 1:size(cells, 2)
        counts[cells[c, k]] += 1
    end
    offsets = zeros(Int, nPts + 1)
    for i in 1:nPts
        offsets[i + 1] = offsets[i] + counts[i]
    end
    fill = copy(offsets[1:nPts])
    indices = Vector{Int}(undef, offsets[nPts + 1])
    for c in 1:size(cells, 1), k in 1:size(cells, 2)
        v = cells[c, k]
        fill[v] += 1
        indices[fill[v]] = c - 1
    end
    return offsets, indices
end

# =====================================================================
# 2-D mesher core -- driver code, no package needed
# =====================================================================

# ── The mesher state ───────────────────────────────────────────────
#
# Local vertex i of a cell is opposite to local edge i. Local edge i
# joins local vertices i+1 and i+2, counting around 1, 2, 3. The cell
# lists its vertices in counterclockwise order.
#   cellV[c, i]  vertex number of local vertex i
#   cellN[c, i]  neighbor cell across local edge i. 0 means none.
#   cellS[c, i]  segment number on local edge i. 0 means none.
#   cellIn[c]    1 if the cell is in the domain, 0 if not
# The vertices 1, 2, 3 are the corners of the bounding triangle.

# Return a new, empty mesher with room for the given counts. The mesher is a
# NamedTuple. A field that the code rebinds later is a Ref. Read it as
# M.name[] and set it as M.name[] = value.
function _new_mesher2(capV::Int, capC::Int, capS::Int, maxVertices::Int)
    return (
        nV = Ref{Int}(0),
        nC = Ref{Int}(0),
        nS = Ref{Int}(0),
        capV = Ref{Int}(capV),
        capC = Ref{Int}(capC),
        capS = Ref{Int}(capS),
        maxVertices = maxVertices,
        vx = zeros(Float64, capV),
        vy = zeros(Float64, capV),
        vCell = zeros(Int, capV),
        tmpStart = zeros(Int, capV),
        lastVertex = Ref{Int}(0),
        cellV = Ref{Matrix{Int}}(zeros(Int, capC, 3)),
        cellN = Ref{Matrix{Int}}(zeros(Int, capC, 3)),
        cellS = Ref{Matrix{Int}}(zeros(Int, capC, 3)),
        cellAlive = zeros(Int, capC),
        cellIn = zeros(Int, capC),
        cellMark = zeros(Int, capC),
        freeCells = Int[],
        segV = Ref{Matrix{Int}}(zeros(Int, capS, 2)),
        segMate = zeros(Int, capS),
        segCurve = zeros(Int, capS),
        segAxis = zeros(Int, capS),
        segAlive = zeros(Int, capS),
        obsKind = Ref{Vector{Int}}(Int[]),
        obsCx = Ref{Vector{Float64}}(Float64[]),
        obsCy = Ref{Vector{Float64}}(Float64[]),
        obsHx = Ref{Vector{Float64}}(Float64[]),
        obsHy = Ref{Vector{Float64}}(Float64[]),
        nObs = Ref{Int}(0),
        hMax = Ref{Float64}(1.0),
        hObs = Ref{Float64}(1.0),
        grading = Ref{Float64}(0.0),
        ratioBound = Ref{Float64}(1.0),
        doBad = Ref{Bool}(false),
        stamp = Ref{Int}(0),
        cav = Int[],
        stack = Int[],
        viol = Int[],
        newCells = Int[],
        qA = _new_queue(),
        qB = _new_queue(),
        wSeg = Int[],
        wCell = Int[],
        wFlag = Int[],
        wTx = Float64[],
        wTy = Float64[],
        ccx = Float64[],
        ccy = Float64[],
        ccR = Float64[],
        hAt = Float64[],
        ratio = Float64[],
        sizeRatio = Float64[],
        iterLimit = 100 * maxVertices + 10_000_000,
    )
end


# Make sure the work arrays of the kernels hold at least n entries.
function _ws2!(M::NamedTuple, n::Int)
    if length(M.wSeg) < n
        m = max(n, 2 * length(M.wSeg), 64)
        resize!(M.wSeg, m); resize!(M.wCell, m); resize!(M.wFlag, m)
        resize!(M.wTx, m); resize!(M.wTy, m)
        resize!(M.ccx, m); resize!(M.ccy, m); resize!(M.ccR, m)
        resize!(M.hAt, m); resize!(M.ratio, m); resize!(M.sizeRatio, m)
    end
    return nothing
end

# Add one vertex. Return its number.
function _add_vertex2!(M::NamedTuple, x::Float64, y::Float64)
    if M.nV[] - 3 >= M.maxVertices
        error("the mesh needs more than max_vertices = $(M.maxVertices) vertices. " *
              "Increase mesh_size or increase max_vertices.")
    end
    M.nV[] += 1
    if M.nV[] > M.capV[]
        M.capV[] = 2 * M.capV[]
        resize!(M.vx, M.capV[]); resize!(M.vy, M.capV[])
        resize!(M.vCell, M.capV[]); resize!(M.tmpStart, M.capV[])
    end
    M.vx[M.nV[]] = x
    M.vy[M.nV[]] = y
    M.vCell[M.nV[]] = 0
    return M.nV[]
end

# Take one cell slot. Reuse a free slot if one exists.
function _alloc_cell2!(M::NamedTuple)
    if !isempty(M.freeCells)
        c = pop!(M.freeCells)
    else
        M.nC[] += 1
        if M.nC[] > M.capC[]
            newCap = 2 * M.capC[]
            M.cellV[] = _grow_matrix(M.cellV[], newCap)
            M.cellN[] = _grow_matrix(M.cellN[], newCap)
            M.cellS[] = _grow_matrix(M.cellS[], newCap)
            resize!(M.cellAlive, newCap)
            resize!(M.cellIn, newCap)
            resize!(M.cellMark, newCap)
            M.capC[] = newCap
        end
        c = M.nC[]
    end
    M.cellAlive[c] = 1
    M.cellIn[c] = 1
    M.cellMark[c] = 0
    for i in 1:3
        M.cellN[][c, i] = 0
        M.cellS[][c, i] = 0
    end
    return c
end

# Add one segment from vertex a to vertex b. Return its number.
function _add_segment2!(M::NamedTuple, a::Int, b::Int, curve::Int, axis::Int)
    M.nS[] += 1
    if M.nS[] > M.capS[]
        newCap = 2 * M.capS[]
        M.segV[] = _grow_matrix(M.segV[], newCap)
        resize!(M.segMate, newCap); resize!(M.segCurve, newCap)
        resize!(M.segAxis, newCap); resize!(M.segAlive, newCap)
        M.capS[] = newCap
    end
    s = M.nS[]
    M.segV[][s, 1] = a
    M.segV[][s, 2] = b
    M.segMate[s] = 0
    M.segCurve[s] = curve
    M.segAxis[s] = axis
    M.segAlive[s] = 1
    return s
end

# Return the local index of the neighbor slot of cell n that points to cell c.
function _slot_of2(M::NamedTuple, n::Int, c::Int)
    for j in 1:3
        M.cellN[][n, j] == c && return j
    end
    error("internal error: the cell adjacency is broken")
end

# ── Point location ───────────────────────────────────────────────────────

# Walk from cell start to the cell that holds the point (x, y). The walk
# crosses an edge when the point lies on the other side of that edge.
# Return a code and a number.
#   1, c   the point lies in cell c or on its boundary
#   2, s   the walk reached an edge of segment s and stopAtSeg is true
#   3, v   the point equals vertex v
#   4, 0   the point lies outside the bounding triangle
function _locate2(M::NamedTuple, x::Float64, y::Float64, start::Int, stopAtSeg::Bool)
    c = start
    steps = 0
    limit = 10 * M.nC[] + 1000
    while true
        steps += 1
        steps > limit && error("internal error: point location does not end")
        moved = false
        for k in 0:2
            i = mod(steps + k, 3) + 1
            u = M.cellV[][c, i % 3 + 1]
            w = M.cellV[][c, (i + 1) % 3 + 1]
            if _orient2d(M.vx[u], M.vy[u], M.vx[w], M.vy[w], x, y) < 0
                n = M.cellN[][c, i]
                n == 0 && return 4, 0
                if stopAtSeg && M.cellS[][c, i] != 0
                    return 2, M.cellS[][c, i]
                end
                c = n
                moved = true
                break
            end
        end
        if !moved
            for i in 1:3
                v = M.cellV[][c, i]
                if M.vx[v] == x && M.vy[v] == y
                    return 3, v
                end
            end
            return 1, c
        end
    end
end

# Find the edge from vertex a to vertex b. Return the cell and the local
# edge index. Return 0, 0 when the edge does not exist.
function _find_edge2(M::NamedTuple, a::Int, b::Int)
    c0 = M.vCell[a]
    c0 == 0 && return 0, 0
    c = c0
    while true
        k = M.cellV[][c, 1] == a ? 1 : (M.cellV[][c, 2] == a ? 2 : 3)
        # Edge (a, u) is opposite the third vertex of the cell.
        if M.cellV[][c, k % 3 + 1] == b
            return c, (k + 1) % 3 + 1
        end
        # Edge (w, a) is opposite the second vertex of the cell.
        if M.cellV[][c, (k + 1) % 3 + 1] == b
            return c, k % 3 + 1
        end
        # Turn counterclockwise around a, across the edge (w, a).
        c = M.cellN[][c, k % 3 + 1]
        (c == 0 || c == c0) && return 0, 0
    end
end

# Set the segment number on both sides of the edge of segment s. Return
# false when the edge does not exist.
function _mark_segment2!(M::NamedTuple, s::Int)
    c, i = _find_edge2(M, M.segV[][s, 1], M.segV[][s, 2])
    c == 0 && return false
    M.cellS[][c, i] = s
    n = M.cellN[][c, i]
    if n != 0
        M.cellS[][n, _slot_of2(M, n, c)] = s
    end
    return true
end

# ── Insertion of one vertex (Bowyer-Watson) ──────────────────────────────

# Return the sign of the incircle test of cell c and the point (x, y).
function _incircle_cell2(M::NamedTuple, c::Int, x::Float64, y::Float64)
    a = M.cellV[][c, 1]; b = M.cellV[][c, 2]; d = M.cellV[][c, 3]
    return _incircle(M.vx[a], M.vy[a], M.vx[b], M.vy[b], M.vx[d], M.vy[d], x, y)
end

# Find the cavity of the point (x, y). The cavity is the set of cells
# whose circumcircle holds the point strictly. The function fills M.cav.
# mode 0: plain insertion. The function never rejects.
# mode 1: insertion of a circumcenter. The function rejects in two
#   cases. A segment lies inside the cavity. Or the point encroaches on
#   a segment on the cavity boundary.
# mode 2: insertion of the split point of segment splitSeg. The function
#   rejects when a segment other than splitSeg lies inside the cavity.
# A rejection returns 1. M.viol then holds the segments to split first.
# A success returns 0.
function _cavity2!(M::NamedTuple, x::Float64, y::Float64, c0::Int, mode::Int, splitSeg::Int)
    empty!(M.viol)
    empty!(M.cav)
    empty!(M.stack)
    M.stamp[] += 1
    st = M.stamp[]
    push!(M.cav, c0)
    push!(M.stack, c0)
    M.cellMark[c0] = st
    while !isempty(M.stack)
        c = pop!(M.stack)
        for i in 1:3
            n = M.cellN[][c, i]
            n == 0 && continue
            mk = M.cellMark[n]
            if mk != st && mk != -st
                if _incircle_cell2(M, n, x, y) > 0
                    M.cellMark[n] = st
                    push!(M.cav, n)
                    push!(M.stack, n)
                    mk = st
                else
                    M.cellMark[n] = -st
                    mk = -st
                end
            end
            # A segment between two cells of the cavity would disappear.
            if mk == st
                s = M.cellS[][c, i]
                if s != 0 && (mode == 1 || (mode == 2 && s != splitSeg))
                    push!(M.viol, s)
                end
            end
        end
    end
    if mode == 1
        # Test the segments on the cavity boundary for encroachment.
        nb = 0
        for c in M.cav
            for i in 1:3
                s = M.cellS[][c, i]
                s == 0 && continue
                n = M.cellN[][c, i]
                if n == 0 || M.cellMark[n] != st
                    nb += 1
                    _ws2!(M, nb)
                    M.wSeg[nb] = s
                    M.wTx[nb] = x
                    M.wTy[nb] = y
                end
            end
        end
        if nb > 0
            segment_encroach_2d_kernel(M.segV[], M.vx, M.vy, M.wSeg, M.wTx, M.wTy, nb, M.wFlag)
            for k in 1:nb
                M.wFlag[k] == 1 && push!(M.viol, M.wSeg[k])
            end
        end
    end
    return isempty(M.viol) ? 0 : 1
end

# Replace the cavity in M.cav by a fan of new cells around a new vertex
# at (x, y). The new cells are in M.newCells. A new cell keeps the
# in-domain flag of the cavity cell that it replaces. One case differs.
# The edge of the segment splitSeg can lie on the cavity boundary. This
# happens when the split point lies slightly off the old edge. The new
# cell next to that edge then takes the flag of the cell on the other
# side of the old edge. The old edge stops being a segment.
function _commit2!(M::NamedTuple, x::Float64, y::Float64, splitSeg::Int)
    st = M.stamp[]
    v = _add_vertex2!(M, x, y)
    empty!(M.newCells)
    for c in M.cav
        for i in 1:3
            n = M.cellN[][c, i]
            if n != 0 && M.cellMark[n] == st
                continue
            end
            a = M.cellV[][c, i % 3 + 1]
            b = M.cellV[][c, (i + 1) % 3 + 1]
            s = M.cellS[][c, i]
            flag = M.cellIn[c]
            slot = n == 0 ? 0 : _slot_of2(M, n, c)
            if s != 0 && s == splitSeg
                s = 0
                if n != 0
                    flag = M.cellIn[n]
                    M.cellS[][n, slot] = 0
                end
            end
            nc = _alloc_cell2!(M)
            M.cellV[][nc, 1] = v
            M.cellV[][nc, 2] = a
            M.cellV[][nc, 3] = b
            M.cellN[][nc, 1] = n
            M.cellS[][nc, 1] = s
            M.cellIn[nc] = flag
            if n != 0
                M.cellN[][n, slot] = nc
            end
            M.tmpStart[a] = nc
            push!(M.newCells, nc)
        end
    end
    # Link the new cells to each other around the new vertex.
    for nc in M.newCells
        b = M.cellV[][nc, 3]
        t = M.tmpStart[b]
        M.cellN[][nc, 2] = t
        M.cellN[][t, 3] = nc
    end
    for nc in M.newCells
        M.vCell[M.cellV[][nc, 1]] = nc
        M.vCell[M.cellV[][nc, 2]] = nc
        M.vCell[M.cellV[][nc, 3]] = nc
    end
    for c in M.cav
        M.cellAlive[c] = 0
        push!(M.freeCells, c)
    end
    M.lastVertex[] = v
    return v
end

# Locate the point and find its cavity. Do not change the mesh. Return 0
# when the insertion can happen and 1 when it cannot. In the second case
# M.viol holds the segments to split first. The cell number of the
# located cell is in M.stack[1] only during the call.
function _dry_run2!(M::NamedTuple, x::Float64, y::Float64, start::Int, mode::Int, splitSeg::Int)
    code, id = _locate2(M, x, y, start, mode == 1)
    if code == 2
        empty!(M.viol)
        push!(M.viol, id)
        return 1
    elseif code == 3
        error("two vertices coincide at ($x, $y). Increase mesh_size or check the obstacles.")
    elseif code == 4
        error("internal error: a vertex lies outside the bounding triangle")
    end
    return _cavity2!(M, x, y, id, mode, splitSeg)
end

# Insert the vertex (x, y). Return 0 on success and 1 on rejection. The
# new vertex number is M.lastVertex after a success.
function _insert2!(M::NamedTuple, x::Float64, y::Float64, start::Int, mode::Int, splitSeg::Int)
    _dry_run2!(M, x, y, start, mode, splitSeg) != 0 && return 1
    _commit2!(M, x, y, splitSeg)
    return 0
end

# ── Scan new cells for work ──────────────────────────────────────────────

# Examine the cells ids[1:n]. Queue each segment that an apex vertex
# encroaches on. When M.doBad is true, queue each cell in the domain
# that is too large or has an angle that is too small.
function _scan_cells2!(M::NamedTuple, ids::Vector{Int}, n::Int)
    nPair = 0
    for k in 1:n
        c = ids[k]
        for i in 1:3
            s = M.cellS[][c, i]
            if s != 0
                nPair += 1
                _ws2!(M, nPair)
                apex = M.cellV[][c, i]
                M.wSeg[nPair] = s
                M.wTx[nPair] = M.vx[apex]
                M.wTy[nPair] = M.vy[apex]
            end
        end
    end
    if nPair > 0
        segment_encroach_2d_kernel(M.segV[], M.vx, M.vy, M.wSeg, M.wTx, M.wTy, nPair, M.wFlag)
        for k in 1:nPair
            M.wFlag[k] == 1 && _qpush!(M.qA, 2 * M.wSeg[k])
        end
    end
    if M.doBad[]
        nIn = 0
        for k in 1:n
            c = ids[k]
            if M.cellAlive[c] == 1 && M.cellIn[c] == 1
                nIn += 1
                _ws2!(M, nIn)
                M.wCell[nIn] = c
            end
        end
        if nIn > 0
            _bad_cells2!(M, nIn, true)
        end
    end
    return nothing
end

# Compute the quality of the cells M.wCell[1:n]. Queue the bad cells
# when push is true. Return the number of bad cells.
function _bad_cells2!(M::NamedTuple, n::Int, push::Bool)
    cell_circumcenter_2d_kernel(M.cellV[], M.vx, M.vy, M.wCell, n, M.ccx, M.ccy, M.ccR)
    size_field_2d_kernel(M.ccx, M.ccy, n, M.obsKind[], M.obsCx[], M.obsCy[], M.obsHx[], M.obsHy[],
                         M.nObs[], M.hMax[], M.hObs[], M.grading[], M.hAt)
    cell_quality_2d_kernel(M.cellV[], M.vx, M.vy, M.wCell, n, M.ccR, M.hAt, M.ratio, M.sizeRatio)
    nBad = 0
    for k in 1:n
        if M.ratio[k] > M.ratioBound[] || M.sizeRatio[k] > 1.0
            nBad += 1
            push && _qpush!(M.qB, M.wCell[k])
        end
    end
    return nBad
end

# ── Segments ─────────────────────────────────────────────────────────────

# Return true when a vertex encroaches on segment s. Test the apex
# vertex of each of the two cells next to the segment.
function _segment_encroached2(M::NamedTuple, s::Int)
    c, i = _find_edge2(M, M.segV[][s, 1], M.segV[][s, 2])
    c == 0 && return false
    _ws2!(M, 2)
    nPair = 1
    M.wSeg[1] = s
    M.wTx[1] = M.vx[M.cellV[][c, i]]
    M.wTy[1] = M.vy[M.cellV[][c, i]]
    n = M.cellN[][c, i]
    if n != 0
        nPair = 2
        apex = M.cellV[][n, _slot_of2(M, n, c)]
        M.wSeg[2] = s
        M.wTx[2] = M.vx[apex]
        M.wTy[2] = M.vy[apex]
    end
    segment_encroach_2d_kernel(M.segV[], M.vx, M.vy, M.wSeg, M.wTx, M.wTy, nPair, M.wFlag)
    return M.wFlag[1] == 1 || (nPair == 2 && M.wFlag[2] == 1)
end

# Return a cell to start a walk toward a point on segment s.
function _start_cell2(M::NamedTuple, s::Int)
    a = M.segV[][s, 1]
    c, _ = _find_edge2(M, a, M.segV[][s, 2])
    return c != 0 ? c : M.vCell[a]
end

# Return the coordinate of vertex v along the axis of a periodic pair.
# Axis 1 means x. Axis 2 means y.
_along(M::NamedTuple, v::Int, axis::Int) = axis == 1 ? M.vx[v] : M.vy[v]

# Split segment s at the point (px, py). The first half keeps the number
# s. The second half gets a new number. Return the number of the second
# half. The function inserts the vertex, then sets the segment flags on
# the edges, then scans the new cells.
function _split_commit2!(M::NamedTuple, s::Int, px::Float64, py::Float64)
    a = M.segV[][s, 1]
    b = M.segV[][s, 2]
    start = _start_cell2(M, s)
    _insert2!(M, px, py, start, 2, s) != 0 &&
        error("internal error: a split insertion failed after a successful dry run")
    v = M.lastVertex[]
    sB = _add_segment2!(M, v, b, M.segCurve[s], M.segAxis[s])
    M.segV[][s, 2] = v
    foundA, foundB = _mark_halves2!(M, v, a, b, s, sB)
    # A half that is not an edge yet returns to the queue. The entry is forced.
    foundA || _qpush!(M.qA, 2 * s + 1)
    foundB || _qpush!(M.qA, 2 * sB + 1)
    _scan_cells2!(M, M.newCells, length(M.newCells))
    return sB
end

# Set the segment numbers sA and sB on the new edges that join the new
# vertex v to the old ends a and b. Return whether each edge exists.
function _mark_halves2!(M::NamedTuple, v::Int, a::Int, b::Int, sA::Int, sB::Int)
    foundA = false
    foundB = false
    for nc in M.newCells
        x = M.cellV[][nc, 2]
        y = M.cellV[][nc, 3]
        # Local edge 3 joins v and x. Local edge 2 joins y and v.
        if x == a
            M.cellS[][nc, 3] = sA; foundA = true
        elseif x == b
            M.cellS[][nc, 3] = sB; foundB = true
        end
        if y == a
            M.cellS[][nc, 2] = sA; foundA = true
        elseif y == b
            M.cellS[][nc, 2] = sB; foundB = true
        end
    end
    return foundA, foundB
end

# Process one entry of queue A. The entry is 2*s + force. A forced
# entry always splits. An entry that is not forced splits only when a
# vertex still encroaches on the segment, or when the segment is not
# an edge.
function _process_segment2!(M::NamedTuple, item::Int)
    s = div(item, 2)
    force = mod(item, 2)
    M.segAlive[s] == 1 || return nothing
    present = _mark_segment2!(M, s)
    if present && force == 0 && !_segment_encroached2(M, s)
        return nothing
    end
    # Compute the split point of segment s, and of its mate when it has one.
    _ws2!(M, 2)
    M.wSeg[1] = s
    segment_split_point_2d_kernel(M.segV[], M.segCurve, M.vx, M.vy, M.obsCx[], M.obsCy[], M.obsHx[], M.wSeg, 1, M.wTx, M.wTy)
    px = M.wTx[1]
    py = M.wTy[1]
    m = M.segMate[s]
    hasMate = m != 0 && M.segAlive[m] == 1
    qx = 0.0
    qy = 0.0
    if hasMate
        mv = M.segV[][m, 1]
        if M.segAxis[s] == 1
            qx = px
            qy = M.vy[mv]
        else
            qx = M.vx[mv]
            qy = py
        end
    end
    # Test both insertions before the first one changes the mesh.
    viol = Int[]
    code = _dry_run2!(M, px, py, _start_cell2(M, s), 2, s)
    append!(viol, M.viol)
    if hasMate
        code2 = _dry_run2!(M, qx, qy, _start_cell2(M, m), 2, m)
        append!(viol, M.viol)
        code = max(code, code2)
    end
    if code != 0
        # Split the segments that block the insertion first. Try again later.
        for t in unique(viol)
            _qpush!(M.qA, 2 * t + 1)
        end
        _qpush!(M.qA, item)
        return nothing
    end
    sB = _split_commit2!(M, s, px, py)
    if hasMate
        mB = _split_commit2!(M, m, qx, qy)
        _pair_halves2!(M, s, sB, m, mB)
    end
    return nothing
end

# Set the mate links between the halves of two segments that split in
# the same way. The halves that cover the same range along the axis are
# mates. The numbers sA and mA are the first halves. The numbers sB and
# mB are the second halves.
function _pair_halves2!(M::NamedTuple, sA::Int, sB::Int, mA::Int, mB::Int)
    axis = M.segAxis[sA]
    outerS = _along(M, M.segV[][sA, 1], axis)
    outerM = _along(M, M.segV[][mA, 1], axis)
    if abs(outerS - outerM) <= 1.0e-12 * (1.0 + abs(outerS))
        M.segMate[sA] = mA; M.segMate[mA] = sA
        M.segMate[sB] = mB; M.segMate[mB] = sB
    else
        M.segMate[sA] = mB; M.segMate[mB] = sA
        M.segMate[sB] = mA; M.segMate[mA] = sB
    end
    return nothing
end

# ── Bad cells ────────────────────────────────────────────────────────────

# Process one cell of queue B. Insert its circumcenter when the cell is
# still bad. When the insertion is not allowed, queue the segments that
# block it and put the cell back at the end of queue B.
function _process_cell2!(M::NamedTuple, c::Int)
    (M.cellAlive[c] == 1 && M.cellIn[c] == 1) || return nothing
    _ws2!(M, 1)
    M.wCell[1] = c
    _bad_cells2!(M, 1, false) == 0 && return nothing
    px = M.ccx[1]
    py = M.ccy[1]
    if !(isfinite(px) && isfinite(py))
        return nothing
    end
    if _insert2!(M, px, py, c, 1, 0) == 0
        _scan_cells2!(M, M.newCells, length(M.newCells))
    else
        for t in unique(M.viol)
            _qpush!(M.qA, 2 * t + 1)
        end
        _qpush!(M.qB, c)
    end
    return nothing
end

# Process queue A until it is empty.
function _run_segment_queue2!(M::NamedTuple)
    iter = 0
    while !_qempty(M.qA)
        _process_segment2!(M, _qpop!(M.qA))
        iter += 1
        iter > M.iterLimit && error("segment recovery does not end. Check the obstacles and mesh_size.")
    end
    return nothing
end

# Process both queues until both are empty. Queue A has priority.
function _run_refinement2!(M::NamedTuple)
    iter = 0
    while true
        if !_qempty(M.qA)
            _process_segment2!(M, _qpop!(M.qA))
        elseif !_qempty(M.qB)
            _process_cell2!(M, _qpop!(M.qB))
        else
            break
        end
        iter += 1
        iter > M.iterLimit && error("refinement does not end. Check the obstacles, mesh_size, and min_angle.")
    end
    return nothing
end

# ── Domain flags ─────────────────────────────────────────────────────────

# Spread the flag 0 from the cells on the stack to their neighbors. The
# spread does not cross an edge of a segment.
function _flood2!(M::NamedTuple, stack::Vector{Int})
    while !isempty(stack)
        c = pop!(stack)
        for i in 1:3
            n = M.cellN[][c, i]
            (n == 0 || M.cellS[][c, i] != 0) && continue
            if M.cellIn[n] == 1
                M.cellIn[n] = 0
                push!(stack, n)
            end
        end
    end
    return nothing
end

# Set cellIn to 1 for each cell in the domain and to 0 for each cell
# outside it. A cell is outside when the flood from the bounding
# triangle or from an obstacle center reaches it.
function _classify_cells2!(M::NamedTuple)
    for c in 1:M.nC[]
        M.cellIn[c] = M.cellAlive[c]
    end
    stack = Int[]
    for c in 1:M.nC[]
        M.cellAlive[c] == 1 || continue
        if M.cellV[][c, 1] <= 3 || M.cellV[][c, 2] <= 3 || M.cellV[][c, 3] <= 3
            M.cellIn[c] = 0
            push!(stack, c)
        end
    end
    _flood2!(M, stack)
    for k in 1:M.nObs[]
        code, id = _locate2(M, M.obsCx[][k], M.obsCy[][k], M.vCell[4], false)
        c = code == 3 ? M.vCell[id] : id
        if M.cellIn[c] == 1
            M.cellIn[c] = 0
            push!(stack, c)
            _flood2!(M, stack)
        end
    end
    return nothing
end

# ── Output of the mesher ─────────────────────────────────────────────────

# Collect the cells in the domain. Renumber the vertices that these cells
# use. Return the vertex array and the cell array.
function _compact2(M::NamedTuple)
    used = zeros(Int, M.nV[])
    nCells = 0
    for c in 1:M.nC[]
        (M.cellAlive[c] == 1 && M.cellIn[c] == 1) || continue
        nCells += 1
        for k in 1:3
            used[M.cellV[][c, k]] = 1
        end
    end
    (used[1] + used[2] + used[3]) == 0 || error("internal error: a cell in the domain uses a corner of the bounding triangle")
    newIdx = zeros(Int, M.nV[])
    nPts = 0
    for v in 1:M.nV[]
        if used[v] == 1
            nPts += 1
            newIdx[v] = nPts
        end
    end
    pts = zeros(Float64, nPts, 2)
    for v in 1:M.nV[]
        if used[v] == 1
            pts[newIdx[v], 1] = M.vx[v]
            pts[newIdx[v], 2] = M.vy[v]
        end
    end
    cells = zeros(Int, nCells, 3)
    k = 0
    for c in 1:M.nC[]
        (M.cellAlive[c] == 1 && M.cellIn[c] == 1) || continue
        k += 1
        for j in 1:3
            cells[k, j] = newIdx[M.cellV[][c, j]]
        end
    end
    return pts, cells
end

# ── Options and obstacles ────────────────────────────────────────────────

# Check the mesh_options record. Set the defaults. Return a NamedTuple.
function _resolve_mesh_options2(raw)
    raw isa AbstractDict || error("'mesh_options' must be a mapping")
    known = ("mesh_size", "mesh_size_obstacle", "grading", "min_angle", "smoothing_sweeps", "max_vertices")
    for k in keys(raw)
        k in known || error("mesh_options: unknown key '$k'. Valid keys: $(join(known, ", "))")
    end
    haskey(raw, "mesh_size") || error("mesh_options: 'mesh_size' is required")
    num(key, default) = begin
        v = get(raw, key, default)
        v isa Real || error("mesh_options['$key'] must be a number, got $(repr(v))")
        Float64(v)
    end
    meshSize = num("mesh_size", nothing)
    meshSize > 0 || error("mesh_options['mesh_size'] must be greater than 0")
    meshSizeObs = num("mesh_size_obstacle", meshSize)
    meshSizeObs > 0 || error("mesh_options['mesh_size_obstacle'] must be greater than 0")
    grading = num("grading", 0.3)
    grading >= 0 || error("mesh_options['grading'] must be 0 or more")
    minAngle = num("min_angle", 20.0)
    (minAngle > 0 && minAngle <= 20.7) || error("mesh_options['min_angle'] must be in (0, 20.7]")
    sweeps = num("smoothing_sweeps", 3)
    (sweeps >= 0 && sweeps == floor(sweeps)) || error("mesh_options['smoothing_sweeps'] must be a whole number, 0 or more")
    maxVertices = num("max_vertices", 2000000)
    (maxVertices >= 100 && maxVertices == floor(maxVertices)) || error("mesh_options['max_vertices'] must be a whole number, 100 or more")
    return (meshSize = meshSize, meshSizeObstacle = meshSizeObs, grading = grading,
            minAngle = minAngle, smoothingSweeps = Int(sweeps), maxVertices = Int(maxVertices))
end

# Check each obstacle record. Raise an error for a bad field, for an
# obstacle that crosses or touches the box, and for two obstacles whose
# bounding boxes overlap. Warn when a gap is smaller than hObs.
function _validate_obstacles2(logBuf::IOBuffer, obstacles, box::Tuple, hObs::Float64)
    isempty(obstacles) && return
    lo = Vector{Vector{Float64}}()
    hi = Vector{Vector{Float64}}()
    for (i, obs) in enumerate(obstacles)
        idx = i - 1
        shape = get(obs, "shape", nothing)
        if !(shape in _obstacle_shapes_2d())
            error("obstacle $idx: shape must be one of $(_obstacle_shapes_2d()) for a 2-D mesh, got '$shape'")
        end
        c = get(obs, "center", nothing)
        (c === nothing || length(c) != 2) && error("obstacle $idx: 'center' must have 2 entries")
        c = Float64.(c)
        if shape == "circle"
            r = get(obs, "radius", nothing)
            (r === nothing || r <= 0) && error("obstacle $idx: 'radius' must be a positive number")
            half = fill(Float64(r), 2)
        else
            he = get(obs, "half_extents", nothing)
            (he === nothing || length(he) != 2 || any(<=(0), he)) && error("obstacle $idx: 'half_extents' must have 2 positive entries")
            half = Float64.(he)
        end
        for d in 1:2
            if c[d] - half[d] <= 0 || c[d] + half[d] >= box[d]
                error("obstacle $idx: crosses or touches the box boundary along axis $(_dirs()[d])")
            end
            gap = min(c[d] - half[d], box[d] - c[d] - half[d])
            gap < hObs && _log(logBuf, "  [warn] obstacle $idx: the gap to the box along axis $(_dirs()[d]) is smaller than mesh_size_obstacle")
        end
        push!(lo, c .- half)
        push!(hi, c .+ half)
    end
    n = length(lo)
    for i in 1:n, j in (i + 1):n
        separation = maximum(max(lo[j][d] - hi[i][d], lo[i][d] - hi[j][d]) for d in 1:2)
        if separation <= 0
            error("obstacles $(i-1) and $(j-1): the bounding boxes overlap")
        elseif separation < hObs
            _log(logBuf, "  [warn] obstacles $(i-1) and $(j-1): the gap is smaller than mesh_size_obstacle")
        end
    end
end

# Build the record of one obstacle for the HDF5 file.
function _obstacle_record2(obs, marker::Int)
    center = Tuple(Float64.(obs["center"]))
    if obs["shape"] == "circle"
        return (shape=obs["shape"], center=center, radius=Float64(obs["radius"]), marker=marker)
    else
        return (shape=obs["shape"], center=center, half_extents=Tuple(Float64.(obs["half_extents"])), marker=marker)
    end
end

# ── The 2-D mesh builder ─────────────────────────────────────────────────

# Return the largest radius-edge ratio of all cells. The caller passes the
# vertex array and the cell array.
function _max_ratio2(pts::Matrix{Float64}, cells::Matrix{Int})
    m = size(cells, 1)
    vx = pts[:, 1]
    vy = pts[:, 2]
    ids = collect(1:m)
    ccx = zeros(Float64, m); ccy = zeros(Float64, m); ccR = zeros(Float64, m)
    hAt = ones(Float64, m); ratio = zeros(Float64, m); sizeRatio = zeros(Float64, m)
    cell_circumcenter_2d_kernel(cells, vx, vy, ids, m, ccx, ccy, ccR)
    cell_quality_2d_kernel(cells, vx, vy, ids, m, ccR, hAt, ratio, sizeRatio)
    return maximum(ratio)
end

# Build the triangle mesh of a box with obstacles. Return the vertex
# array, the cell array, the edge array, the cell-to-edge array, the
# vertex markers, the edge markers, and a NamedTuple of statistics. The function needs no package.
#
# Lx, Ly: box lengths
# perX, perY: true when the box is periodic along x or y
# obsKind, obsCx, obsCy, obsHx, obsHy: flat obstacle arrays. Kind 1 is a
#   circle. Kind 2 is a rectangle.
# hMax, hObs, grading, minAngle, sweeps, maxVertices: the mesh options
function _triangulate_domain2(logBuf::IOBuffer, Lx::Float64, Ly::Float64, perX::Bool, perY::Bool,
                             obsKind::Vector{Int}, obsCx::Vector{Float64}, obsCy::Vector{Float64},
                             obsHx::Vector{Float64}, obsHy::Vector{Float64},
                             hMax::Float64, hObs::Float64, grading::Float64,
                             minAngle::Float64, sweeps::Int, maxVertices::Int)
    nObs = length(obsKind)
    hObsUsed = nObs == 0 ? hMax : hObs

    # ---- Count the boundary vertices.
    nx = max(1, ceil(Int, Lx / hMax))
    ny = max(1, ceil(Int, Ly / hMax))
    nBox = 2 * nx + 2 * ny
    nPer = zeros(Int, nObs)
    mxs = zeros(Int, nObs)
    mys = zeros(Int, nObs)
    for k in 1:nObs
        if obsKind[k] == 1
            nPer[k] = max(8, ceil(Int, 2.0 * pi * obsHx[k] / hObsUsed))
        else
            mxs[k] = max(1, ceil(Int, 2.0 * obsHx[k] / hObsUsed))
            mys[k] = max(1, ceil(Int, 2.0 * obsHy[k] / hObsUsed))
            nPer[k] = 2 * mxs[k] + 2 * mys[k]
        end
    end
    nBnd = nBox + sum(nPer)

    # ---- Fill the boundary vertices with the kernels.
    px = zeros(Float64, nx + 1)
    py = zeros(Float64, ny + 1)
    stretch_positions_kernel(nx, 1.0, 0, px)
    stretch_positions_kernel(ny, 1.0, 0, py)
    for i in 1:nx + 1
        px[i] *= Lx
    end
    for j in 1:ny + 1
        py[j] *= Ly
    end
    bx = zeros(Float64, nBnd)
    by = zeros(Float64, nBnd)
    box_side_vertices_kernel(nx, ny, px, py, bx, by)
    offset = nBox
    for k in 1:nObs
        if obsKind[k] == 1
            circle_vertices_kernel(obsCx[k], obsCy[k], obsHx[k], nPer[k], offset, bx, by)
        else
            rectangle_vertices_kernel(obsCx[k], obsCy[k], obsHx[k], obsHy[k], mxs[k], mys[k], offset, bx, by)
        end
        offset += nPer[k]
    end

    # ---- Start the mesher with the bounding triangle.
    M = _new_mesher2(max(1024, 8 * nBnd), max(2048, 16 * nBnd), max(256, 4 * nBnd), maxVertices)
    M.obsKind[] = obsKind; M.obsCx[] = obsCx; M.obsCy[] = obsCy; M.obsHx[] = obsHx; M.obsHy[] = obsHy
    M.nObs[] = nObs
    M.hMax[] = hMax; M.hObs[] = hObsUsed; M.grading[] = grading
    M.ratioBound[] = 1.0 / (2.0 * sind(minAngle))
    ccx = 0.5 * Lx
    ccy = 0.5 * Ly
    bigR = 1000.0 * hypot(Lx, Ly)
    _add_vertex2!(M, ccx - 4.0 * bigR, ccy - 3.0 * bigR)
    _add_vertex2!(M, ccx + 4.0 * bigR, ccy - 3.0 * bigR)
    _add_vertex2!(M, ccx, ccy + 5.0 * bigR)
    c1 = _alloc_cell2!(M)
    M.cellV[][c1, 1] = 1; M.cellV[][c1, 2] = 2; M.cellV[][c1, 3] = 3
    M.vCell[1] = c1; M.vCell[2] = c1; M.vCell[3] = c1

    # ---- Insert the boundary vertices. Vertex 3 + i is boundary vertex i.
    start = c1
    for i in 1:nBnd
        _insert2!(M, bx[i], by[i], start, 0, 0) == 0 || error("internal error: plain insertion failed")
        start = M.vCell[M.lastVertex[]]
    end

    # ---- Build the segments. Each loop closes on its first vertex.
    for k in 1:nBox
        a = 3 + k
        b = 3 + (k == nBox ? 1 : k + 1)
        axis = (k <= nx || (k > nx + ny && k <= 2 * nx + ny)) ? 1 : 2
        _add_segment2!(M, a, b, 0, axis)
    end
    if perY
        for k in 1:nx
            s = k
            m = nx + ny + (nx - k + 1)
            M.segMate[s] = m; M.segMate[m] = s
        end
    end
    if perX
        for j in 1:ny
            s = nx + j
            m = 2 * nx + ny + (ny - j + 1)
            M.segMate[s] = m; M.segMate[m] = s
        end
    end
    base = nBox
    for k in 1:nObs
        for i in 1:nPer[k]
            a = 3 + base + i
            b = 3 + base + (i == nPer[k] ? 1 : i + 1)
            _add_segment2!(M, a, b, obsKind[k] == 1 ? k : 0, 0)
        end
        base += nPer[k]
    end

    # ---- Make every segment an edge of the mesh. Remove encroachment.
    for s in 1:M.nS[]
        _mark_segment2!(M, s) || _qpush!(M.qA, 2 * s + 1)
    end
    ids = Int[c for c in 1:M.nC[] if M.cellAlive[c] == 1]
    _scan_cells2!(M, ids, length(ids))
    while true
        _run_segment_queue2!(M)
        nMissing = 0
        for s in 1:M.nS[]
            M.segAlive[s] == 1 || continue
            if !_mark_segment2!(M, s)
                _qpush!(M.qA, 2 * s + 1)
                nMissing += 1
            end
        end
        nMissing == 0 && break
    end

    # ---- Set the domain flags. Then refine the cells in the domain.
    _classify_cells2!(M)
    M.doBad[] = true
    ids = Int[c for c in 1:M.nC[] if M.cellAlive[c] == 1 && M.cellIn[c] == 1]
    isempty(ids) && error("the domain has no cell. Check the box and the obstacles.")
    _scan_cells2!(M, ids, length(ids))
    _run_refinement2!(M)
    for s in 1:M.nS[]
        M.segAlive[s] == 1 || continue
        _mark_segment2!(M, s) || error("internal error: a segment is not an edge after refinement")
    end

    # ---- Remove the cells outside the domain.
    pts, cells = _compact2(M)
    nPts = size(pts, 1)
    nCells = size(cells, 1)
    meas = zeros(Float64, nCells)
    signed_measure_2d_kernel(cells, pts, nCells, meas)
    nNeg = zeros(Int, 1)
    count_negative_kernel(meas, nCells, nNeg)
    if nNeg[1] > 0
        _log(logBuf, "  [warn] triangle: $(nNeg[1]) cells had a negative area. The writer reordered them.")
        fix_orientation_kernel(meas, nCells, 3, cells)
    end

    # ---- Edges, boundary edges, and vertex markers.
    edges, c2e = _build_edge_registry(cells, 2)
    nEdges = size(edges, 1)
    edgeCount = zeros(Int, nEdges)
    for ci in 1:nCells, li in 1:3
        edgeCount[c2e[ci, li] + 1] += 1
    end
    isBnd = zeros(Int, nEdges)
    for e in 1:nEdges
        edgeCount[e] == 1 && (isBnd[e] = 1)
    end
    tol = 1.0e-9 * max(Lx, Ly)
    ptMk = zeros(Int32, nPts)
    point_markers_2d_domain_kernel(pts, nPts, Lx, Ly, tol, obsKind, obsCx, obsCy, obsHx, obsHy, nObs, ptMk)

    # ---- Smooth the interior vertices with Jacobi sweeps. A sweep stays
    # only when no cell flips and the worst radius-edge ratio does not grow.
    deg = zeros(Int, nPts)
    for e in 1:nEdges
        deg[edges[e, 1]] += 1
        deg[edges[e, 2]] += 1
    end
    nbrOff = ones(Int, nPts + 1)
    for v in 1:nPts
        nbrOff[v + 1] = nbrOff[v] + deg[v]
    end
    nbrIdx = zeros(Int, nbrOff[nPts + 1] - 1)
    cursor = copy(nbrOff[1:nPts])
    for e in 1:nEdges
        u = edges[e, 1]; w = edges[e, 2]
        nbrIdx[cursor[u]] = w; cursor[u] += 1
        nbrIdx[cursor[w]] = u; cursor[w] += 1
    end
    avgPts = zeros(Float64, nPts, 2)
    newPts = zeros(Float64, nPts, 2)
    worst = _max_ratio2(pts, cells)
    sweepsDone = 0
    for _ in 1:sweeps
        smooth_jacobi_2d_kernel(pts, nbrOff, nbrIdx, nPts, avgPts)
        blend_positions_2d_kernel(pts, avgPts, ptMk, nPts, 0.5, newPts)
        signed_measure_2d_kernel(cells, newPts, nCells, meas)
        nNeg[1] = 0
        count_negative_kernel(meas, nCells, nNeg)
        nNeg[1] > 0 && break
        worstNew = _max_ratio2(newPts, cells)
        worstNew > worst && break
        pts, newPts = newPts, pts
        worst = worstNew
        sweepsDone += 1
    end

    # ---- Edge markers.
    edgeMk = zeros(Int32, nEdges)
    edge_markers_2d_domain_kernel(pts, edges, isBnd, nEdges, Lx, Ly, tol, obsKind, obsCx, obsCy, obsHx, obsHy, nObs, edgeMk)

    info = (nVertices = nPts, nCells = nCells, nEdges = nEdges, maxRatio = worst,
            minAngleDeg = asind(1.0 / (2.0 * worst)), sweepsDone = sweepsDone)
    return pts, cells, edges, c2e, ptMk, edgeMk, info
end

# =====================================================================
# 3-D mesher core -- driver code, no package needed
# =====================================================================

# ── The mesher state ───────────────────────────────────────────────
#
# Local vertex i of a cell is opposite to local face i. The cell lists
# its vertices so that _orient3d is positive. The function _face3 lists the
# vertices of each face. Local face i has the vertices _face3(i) in an
# order that gives its normal a direction away from the cell.
#   cellV[c, i]  vertex number of local vertex i
#   cellN[c, i]  neighbor cell across local face i. 0 means none.
#   cellF[c, i]  number of the subfacet on local face i. 0 means none.
#   cellIn[c]    1 if the cell is in the domain, 0 if not
# The vertices 1 to 4 are the corners of the bounding tetrahedron.

# Return the local vertices of local face i (1 to 4) of a tetrahedron.
@inline function _face3(i::Int)
    i == 1 && return (2, 3, 4)
    i == 2 && return (1, 4, 3)
    i == 3 && return (1, 2, 4)
    return (1, 3, 2)
end

# Return a new, empty mesher with room for the given counts. The mesher is a
# NamedTuple. A field that the code rebinds later is a Ref. Read it as
# M.name[] and set it as M.name[] = value.
function _new_mesher3(capV::Int, capC::Int, capS::Int, capSF::Int, maxVertices::Int)
    return (
        insCount = zeros(Int, 2),
        nV = Ref{Int}(0),
        nC = Ref{Int}(0),
        nS = Ref{Int}(0),
        nSF = Ref{Int}(0),
        capV = Ref{Int}(capV),
        capC = Ref{Int}(capC),
        capS = Ref{Int}(capS),
        capSF = Ref{Int}(capSF),
        maxVertices = maxVertices,
        vx = zeros(Float64, capV),
        vy = zeros(Float64, capV),
        vz = zeros(Float64, capV),
        vCell = zeros(Int, capV),
        vLabel = zeros(Int, capV),
        vOnSeg = zeros(Int, capV),
        tmpSF = zeros(Int, capV),
        lastVertex = Ref{Int}(0),
        cellV = Ref{Matrix{Int}}(zeros(Int, capC, 4)),
        cellN = Ref{Matrix{Int}}(zeros(Int, capC, 4)),
        cellF = Ref{Matrix{Int}}(zeros(Int, capC, 4)),
        cellAlive = zeros(Int, capC),
        cellIn = zeros(Int, capC),
        cellMark = zeros(Int, capC),
        freeCells = Int[],
        segV = Ref{Matrix{Int}}(zeros(Int, capS, 2)),
        segFac = Ref{Matrix{Int}}(zeros(Int, capS, 2)),
        segAxis = zeros(Int, capS),
        segGroup = zeros(Int, capS),
        segAlive = zeros(Int, capS),
        segOf = Dict{Tuple{Int,Int},Int}(),
        groups = Vector{Int}[],
        facLabel = Int[],
        facDrop = Int[],
        facMate = Int[],
        facMateAxis = Int[],
        facMateCoord = Float64[],
        facSF = Int[],
        sfV = Ref{Matrix{Int}}(zeros(Int, capSF, 3)),
        sfN = Ref{Matrix{Int}}(zeros(Int, capSF, 3)),
        sfFac = zeros(Int, capSF),
        sfAlive = zeros(Int, capSF),
        sfMark = zeros(Int, capSF),
        obsKind = Ref{Vector{Int}}(Int[]),
        obsCx = Ref{Vector{Float64}}(Float64[]),
        obsCy = Ref{Vector{Float64}}(Float64[]),
        obsCz = Ref{Vector{Float64}}(Float64[]),
        obsHx = Ref{Vector{Float64}}(Float64[]),
        obsHy = Ref{Vector{Float64}}(Float64[]),
        obsHz = Ref{Vector{Float64}}(Float64[]),
        nObs = Ref{Int}(0),
        hMax = Ref{Float64}(1.0),
        hObs = Ref{Float64}(1.0),
        grading = Ref{Float64}(0.0),
        ratioBound = Ref{Float64}(2.0),
        doBad = Ref{Bool}(false),
        stamp = Ref{Int}(0),
        cav = Int[],
        stack = Int[],
        star = Int[],
        violSeg = Int[],
        violSF = Int[],
        lostSeg = Int[],
        lostSF = Int[],
        newCells = Int[],
        newSF = Int[],
        cavSF = Int[],
        hBits = Ref{Int}(10),
        hK = Ref{Vector{Int}}(zeros(Int, 1 << 10)),
        hT = Ref{Vector{Int}}(zeros(Int, 1 << 10)),
        hL = Ref{Vector{Int}}(zeros(Int, 1 << 10)),
        hS = Ref{Vector{Int}}(zeros(Int, 1 << 10)),
        hStamp = Ref{Int}(0),
        qA = _new_queue(),
        qB = _new_queue(),
        qC = _new_queue(),
        wSeg = Int[],
        wCell = Int[],
        wFlag = Int[],
        wTx = Float64[],
        wTy = Float64[],
        wTz = Float64[],
        ccx = Float64[],
        ccy = Float64[],
        ccz = Float64[],
        ccR = Float64[],
        hAt = Float64[],
        ratio = Float64[],
        sizeRatio = Float64[],
        iterLimit = 200 * maxVertices + 10_000_000,
        qFine = Ref{Float64}(0.0),
    )
end

# Round a coordinate to the fine grid. The grid spacing is a power of two.
# The fine grid is much finer than the lattice of the input, but much
# coarser than the rounding error of a computed circumcenter. A value
# that is exactly on the grid in exact arithmetic stays exactly on it.
# Mathematically degenerate configurations, such as four coplanar or
# cospherical lattice points, then stay exactly degenerate in floating
# point. A flat cell cannot appear from a rounding error.
@inline function _quantize3(M::NamedTuple, x::Float64)
    return round(x / M.qFine[]) * M.qFine[]
end

# Round the three coordinates of a point to the fine grid.
@inline function _quantize_point3(M::NamedTuple, x::Float64, y::Float64, z::Float64)
    return _quantize3(M, x), _quantize3(M, y), _quantize3(M, z)
end


# Make the hash table for the linking of new cells hold 2^bits slots.
function _hash_resize3!(M::NamedTuple, bits::Int)
    M.hBits[] = bits
    n = 1 << bits
    M.hK[] = zeros(Int, n)
    M.hT[] = zeros(Int, n)
    M.hL[] = zeros(Int, n)
    M.hS[] = zeros(Int, n)
    return nothing
end

# Make sure the work arrays of the kernels hold at least n entries.
function _ws3!(M::NamedTuple, n::Int)
    if length(M.wSeg) < n
        m = max(n, 2 * length(M.wSeg), 64)
        resize!(M.wSeg, m); resize!(M.wCell, m); resize!(M.wFlag, m)
        resize!(M.wTx, m); resize!(M.wTy, m); resize!(M.wTz, m)
        resize!(M.ccx, m); resize!(M.ccy, m); resize!(M.ccz, m); resize!(M.ccR, m)
        resize!(M.hAt, m); resize!(M.ratio, m); resize!(M.sizeRatio, m)
    end
    return nothing
end

# Add one vertex. Return its number.
function _add_vertex3!(M::NamedTuple, x::Float64, y::Float64, z::Float64, label::Int, onSeg::Int)
    if M.nV[] - 4 >= M.maxVertices
        error("the mesh needs more than max_vertices = $(M.maxVertices) vertices. " *
              "Increase mesh_size or increase max_vertices.")
    end
    M.nV[] += 1
    if M.nV[] > M.capV[]
        M.capV[] = 2 * M.capV[]
        resize!(M.vx, M.capV[]); resize!(M.vy, M.capV[]); resize!(M.vz, M.capV[])
        resize!(M.vCell, M.capV[]); resize!(M.vLabel, M.capV[])
        resize!(M.vOnSeg, M.capV[]); resize!(M.tmpSF, M.capV[])
    end
    v = M.nV[]
    M.vx[v] = x
    M.vy[v] = y
    M.vz[v] = z
    M.vCell[v] = 0
    M.vLabel[v] = label
    M.vOnSeg[v] = onSeg
    M.tmpSF[v] = 0
    return v
end

# Take one cell slot. Reuse a free slot if one exists.
function _alloc_cell3!(M::NamedTuple)
    if !isempty(M.freeCells)
        c = pop!(M.freeCells)
    else
        M.nC[] += 1
        if M.nC[] > M.capC[]
            newCap = 2 * M.capC[]
            M.cellV[] = _grow_matrix(M.cellV[], newCap)
            M.cellN[] = _grow_matrix(M.cellN[], newCap)
            M.cellF[] = _grow_matrix(M.cellF[], newCap)
            resize!(M.cellAlive, newCap)
            resize!(M.cellIn, newCap)
            resize!(M.cellMark, newCap)
            M.capC[] = newCap
        end
        c = M.nC[]
    end
    M.cellAlive[c] = 1
    M.cellIn[c] = 1
    M.cellMark[c] = 0
    for i in 1:4
        M.cellN[][c, i] = 0
        M.cellF[][c, i] = 0
    end
    return c
end

# Return the local index of the neighbor slot of cell n that points to cell c.
function _slot_of3(M::NamedTuple, n::Int, c::Int)
    for j in 1:4
        M.cellN[][n, j] == c && return j
    end
    error("internal error: the cell adjacency is broken")
end

# Return the number of the live subfacet on local face i of cell c.
# Return 0 when the face has none. A number that points to a dead
# subfacet is stale. The function clears a stale number.
function _live_sf3(M::NamedTuple, c::Int, i::Int)
    s = M.cellF[][c, i]
    s == 0 && return 0
    if M.sfAlive[s] == 1
        return s
    end
    M.cellF[][c, i] = 0
    return 0
end

# Return the sign of orient3d for the face (a, b, d) of vertices and the point (x, y, z).
@inline function _o3v3(M::NamedTuple, a::Int, b::Int, d::Int, x::Float64, y::Float64, z::Float64)
    return _orient3d(M.vx[a], M.vy[a], M.vz[a], M.vx[b], M.vy[b], M.vz[b],
                     M.vx[d], M.vy[d], M.vz[d], x, y, z)
end

# Return the sign of the insphere test of cell c and the point (x, y, z).
@inline function _insphere_cell3(M::NamedTuple, c::Int, x::Float64, y::Float64, z::Float64)
    a = M.cellV[][c, 1]; b = M.cellV[][c, 2]; d = M.cellV[][c, 3]; e = M.cellV[][c, 4]
    # insCount[1] counts the calls. insCount[2] counts the calls that
    # needed the exact stage.
    M.insCount[1] += 1
    s = _insphere_fast(M.vx[a], M.vy[a], M.vz[a], M.vx[b], M.vy[b], M.vz[b],
                       M.vx[d], M.vy[d], M.vz[d], M.vx[e], M.vy[e], M.vz[e], x, y, z)
    s != 2 && return s
    M.insCount[2] += 1
    return -_insphere_exact(M.vx[a], M.vy[a], M.vz[a], M.vx[b], M.vy[b], M.vz[b],
                            M.vx[d], M.vy[d], M.vz[d], M.vx[e], M.vy[e], M.vz[e], x, y, z)
end

# ── Point location ───────────────────────────────────────────────────────

# Walk from cell start to the cell that holds the point (x, y, z). The
# walk crosses a face when the point lies on the outer side of that face.
# Return a code and a number.
#   1, c   the point lies in cell c or on its boundary
#   2, s   the walk reached a face of subfacet s and stopAtFace is true
#   3, v   the point equals vertex v
#   4, 0   the point lies outside the bounding tetrahedron
function _locate3(M::NamedTuple, x::Float64, y::Float64, z::Float64, start::Int, stopAtFace::Bool)
    c = start
    steps = 0
    limit = 10 * M.nC[] + 1000
    while true
        steps += 1
        steps > limit && error("internal error: point location does not end")
        moved = false
        for k in 0:3
            i = mod(steps + k, 4) + 1
            f = _face3(i)
            if _o3v3(M, M.cellV[][c, f[1]], M.cellV[][c, f[2]], M.cellV[][c, f[3]], x, y, z) > 0
                n = M.cellN[][c, i]
                n == 0 && return 4, 0
                if stopAtFace
                    s = _live_sf3(M, c, i)
                    s != 0 && return 2, s
                end
                c = n
                moved = true
                break
            end
        end
        if !moved
            for i in 1:4
                v = M.cellV[][c, i]
                if M.vx[v] == x && M.vy[v] == y && M.vz[v] == z
                    return 3, v
                end
            end
            return 1, c
        end
    end
end

# Collect the cells that hold vertex a. The result is in M.star.
function _star3!(M::NamedTuple, a::Int)
    empty!(M.star)
    c0 = M.vCell[a]
    c0 == 0 && return nothing
    M.stamp[] += 1
    st = M.stamp[]
    push!(M.star, c0)
    M.cellMark[c0] = st
    k = 1
    while k <= length(M.star)
        c = M.star[k]
        k += 1
        for i in 1:4
            M.cellV[][c, i] == a && continue
            n = M.cellN[][c, i]
            if n != 0 && M.cellMark[n] != st
                M.cellMark[n] = st
                push!(M.star, n)
            end
        end
    end
    return nothing
end

# Find a cell that holds both vertices a and b. Return 0 if none exists.
function _find_edge3(M::NamedTuple, a::Int, b::Int)
    _star3!(M, a)
    for c in M.star
        for i in 1:4
            M.cellV[][c, i] == b && return c
        end
    end
    return 0
end

# Find the cell face with the vertices a, b, c. Return the cell and the
# local index of the face. Return 0, 0 if no such face exists.
function _find_face3(M::NamedTuple, a::Int, b::Int, d::Int)
    _star3!(M, a)
    for c in M.star
        hasB = 0
        hasD = 0
        for i in 1:4
            v = M.cellV[][c, i]
            v == b && (hasB = i)
            v == d && (hasD = i)
        end
        if hasB != 0 && hasD != 0
            # The face is opposite the one vertex that is not a, b, or d.
            for i in 1:4
                v = M.cellV[][c, i]
                (v != a && v != b && v != d) && return c, i
            end
        end
    end
    return 0, 0
end

# Return the segment number of the edge (u, w). Return 0 if the edge is
# not a segment. The flags filter out nearly all edges before the lookup.
@inline function _seg_of_edge3(M::NamedTuple, u::Int, w::Int)
    (M.vOnSeg[u] == 1 && M.vOnSeg[w] == 1) || return 0
    key = u < w ? (u, w) : (w, u)
    return get(M.segOf, key, 0)
end

# Return true when every cell around the edge (u, w) belongs to the
# cavity of the current stamp st. The cell c holds the edge.
function _ring_interior3(M::NamedTuple, c::Int, u::Int, w::Int, st::Int)
    cur = c
    x = 0
    y = 0
    for i in 1:4
        v = M.cellV[][cur, i]
        if v != u && v != w
            x == 0 ? (x = v) : (y = v)
        end
    end
    count = 0
    while true
        count += 1
        count > 1000 && error("internal error: the ring around an edge does not close")
        lx = 0
        for i in 1:4
            M.cellV[][cur, i] == x && (lx = i)
        end
        n = M.cellN[][cur, lx]
        n == 0 && return false
        M.cellMark[n] != st && return false
        n == c && return true
        z = 0
        for i in 1:4
            v = M.cellV[][n, i]
            (v != u && v != w && v != y) && (z = v)
        end
        x = y
        y = z
        cur = n
    end
end

# ── Insertion of one vertex (Bowyer-Watson) ──────────────────────────────

# Make a place for the keys of n new faces in the hash table.
function _hash_prepare3!(M::NamedTuple, n::Int)
    bits = M.hBits[]
    while (1 << bits) < 4 * n + 16
        bits += 1
    end
    bits != M.hBits[] && _hash_resize3!(M, bits)
    M.hStamp[] += 1
    return nothing
end

# Link the cell nc through the face li to the other new cell that holds
# the same edge (a, b). The first cell that holds an edge stores it. The
# second cell finds it and links both.
function _hash_link3!(M::NamedTuple, nc::Int, li::Int, a::Int, b::Int)
    key = a < b ? (a << 32) | b : (b << 32) | a
    mask = (1 << M.hBits[]) - 1
    h = Int((UInt64(key) * 0x9E3779B97F4A7C15) >> (64 - M.hBits[])) & mask
    while true
        if M.hS[][h + 1] != M.hStamp[]
            M.hS[][h + 1] = M.hStamp[]
            M.hK[][h + 1] = key
            M.hT[][h + 1] = nc
            M.hL[][h + 1] = li
            return nothing
        end
        if M.hK[][h + 1] == key
            other = M.hT[][h + 1]
            ol = M.hL[][h + 1]
            M.cellN[][nc, li] = other
            M.cellN[][other, ol] = nc
            return nothing
        end
        h = (h + 1) & mask
    end
end

# Find the cavity of the point (x, y, z). The cavity is the set of cells
# whose circumsphere holds the point strictly. The function fills M.cav.
# mode 0: plain insertion. The function never rejects.
# mode 1: circumcenter of a cell. A deleted segment or subfacet is
#   forbidden. The function also rejects when the point encroaches on a
#   segment or subfacet on the cavity boundary.
# mode 2: midpoint of the segment splitSeg. The function never rejects.
#   The insertion has priority. A segment or a subfacet that the cavity
#   deletes goes to the lists M.lostSeg and M.lostSF. The caller queues
#   them for recovery. A deleted subfacet is a subfacet that the point
#   encroaches on.
# mode 3: circumcenter of a subfacet of the facet fa1. A deleted segment
#   is forbidden. The function also rejects when the point encroaches on
#   a segment on the cavity boundary. A deleted subfacet goes to the list
#   M.lostSF.
# A rejection returns 1. The lists M.violSeg and M.violSF then hold the
# segments and subfacets to split first. A success returns 0.
function _cavity3!(M::NamedTuple, x::Float64, y::Float64, z::Float64, c0::Int,
                   mode::Int, splitSeg::Int, fa1::Int, fa2::Int)
    empty!(M.violSeg)
    empty!(M.violSF)
    empty!(M.lostSeg)
    empty!(M.lostSF)
    empty!(M.cav)
    empty!(M.stack)
    M.stamp[] += 1
    st = M.stamp[]
    push!(M.cav, c0)
    push!(M.stack, c0)
    M.cellMark[c0] = st
    while !isempty(M.stack)
        c = pop!(M.stack)
        for i in 1:4
            n = M.cellN[][c, i]
            n == 0 && continue
            mk = M.cellMark[n]
            if mk != st && mk != -st
                if _insphere_cell3(M, n, x, y, z) > 0
                    M.cellMark[n] = st
                    push!(M.cav, n)
                    push!(M.stack, n)
                    mk = st
                else
                    M.cellMark[n] = -st
                    mk = -st
                end
            end
        end
    end
    mode == 0 && return 0
    # Deleted subfacets: a face with a subfacet between two cavity cells.
    for c in M.cav
        for i in 1:4
            n = M.cellN[][c, i]
            (n != 0 && M.cellMark[n] == st) || continue
            s = _live_sf3(M, c, i)
            s == 0 && continue
            mode == 1 ? push!(M.violSF, s) : push!(M.lostSF, s)
        end
    end
    nb = 0
    for c in M.cav
        for (la, lb) in _tet_local_edges()
            u = M.cellV[][c, la]
            w = M.cellV[][c, lb]
            s = _seg_of_edge3(M, u, w)
            s == 0 && continue
            if _ring_interior3(M, c, u, w, st)
                if mode == 2
                    s == splitSeg || push!(M.lostSeg, s)
                else
                    push!(M.violSeg, s)
                end
            elseif mode == 1 || mode == 3
                nb += 1
                _ws3!(M, nb)
                M.wSeg[nb] = s
                M.wTx[nb] = x; M.wTy[nb] = y; M.wTz[nb] = z
            end
        end
    end
    if nb > 0
        segment_encroach_3d_kernel(M.segV[], M.vx, M.vy, M.vz, M.wSeg, M.wTx, M.wTy, M.wTz, nb, M.wFlag)
        for k in 1:nb
            M.wFlag[k] == 1 && push!(M.violSeg, M.wSeg[k])
        end
    end
    if mode == 1
        # Subfacets on the cavity boundary that the point encroaches on.
        nb = 0
        for c in M.cav
            for i in 1:4
                n = M.cellN[][c, i]
                (n == 0 || M.cellMark[n] != st) || continue
                s = _live_sf3(M, c, i)
                s == 0 && continue
                nb += 1
                _ws3!(M, nb)
                M.wSeg[nb] = s
                M.wTx[nb] = x; M.wTy[nb] = y; M.wTz[nb] = z
            end
        end
        if nb > 0
            subfacet_encroach_kernel(M.sfV[], M.vx, M.vy, M.vz, M.wSeg, M.wTx, M.wTy, M.wTz, nb, M.wFlag)
            for k in 1:nb
                M.wFlag[k] == 1 && push!(M.violSF, M.wSeg[k])
            end
        end
    end
    return (isempty(M.violSeg) && isempty(M.violSF)) ? 0 : 1
end

# Replace the cavity in M.cav by a fan of new cells around the vertex v.
# The new cells are in M.newCells. A new cell keeps the in-domain flag
# and the subfacet number of the cavity cell that it replaces.
function _commit3!(M::NamedTuple, v::Int)
    st = M.stamp[]
    x = M.vx[v]; y = M.vy[v]; z = M.vz[v]
    empty!(M.newCells)
    _hash_prepare3!(M, 4 * length(M.cav))
    for c in M.cav
        for i in 1:4
            n = M.cellN[][c, i]
            if n != 0 && M.cellMark[n] == st
                continue
            end
            f = _face3(i)
            a = M.cellV[][c, f[1]]
            b = M.cellV[][c, f[2]]
            d = M.cellV[][c, f[3]]
            s = _live_sf3(M, c, i)
            flag = M.cellIn[c]
            slot = n == 0 ? 0 : _slot_of3(M, n, c)
            nc = _alloc_cell3!(M)
            M.cellV[][nc, 1] = v
            M.cellV[][nc, 2] = a
            M.cellV[][nc, 3] = b
            M.cellV[][nc, 4] = d
            M.cellN[][nc, 1] = n
            M.cellF[][nc, 1] = s
            M.cellIn[nc] = flag
            if n != 0
                M.cellN[][n, slot] = nc
            end
            # The faces next to v. Face 2 holds the edge (b, d). Face 3
            # holds the edge (a, d). Face 4 holds the edge (a, b).
            _hash_link3!(M, nc, 2, b, d)
            _hash_link3!(M, nc, 3, a, d)
            _hash_link3!(M, nc, 4, a, b)
            push!(M.newCells, nc)
        end
    end
    for nc in M.newCells
        M.vCell[M.cellV[][nc, 1]] = nc
        M.vCell[M.cellV[][nc, 2]] = nc
        M.vCell[M.cellV[][nc, 3]] = nc
        M.vCell[M.cellV[][nc, 4]] = nc
    end
    for c in M.cav
        M.cellAlive[c] = 0
        push!(M.freeCells, c)
    end
    M.lastVertex[] = v
    return nothing
end

# Locate the point and find its cavity. Do not change the mesh. Return 0
# when the insertion can happen and 1 when it cannot. In the second case
# the lists M.violSeg and M.violSF hold the items to split first.
function _dry3!(M::NamedTuple, x::Float64, y::Float64, z::Float64, start::Int,
                mode::Int, splitSeg::Int, fa1::Int, fa2::Int)
    code, id = _locate3(M, x, y, z, start, mode == 1)
    if code == 2
        empty!(M.violSeg)
        empty!(M.violSF)
        push!(M.violSF, id)
        return 1
    elseif code == 3
        error("two vertices coincide at ($x, $y, $z). Increase mesh_size or check the obstacles.")
    elseif code == 4
        error("internal error: a vertex lies outside the bounding tetrahedron")
    end
    return _cavity3!(M, x, y, z, id, mode, splitSeg, fa1, fa2)
end

# Insert the existing vertex v. Return 0 on success and 1 on rejection.
function _insert3!(M::NamedTuple, v::Int, start::Int, mode::Int, splitSeg::Int, fa1::Int, fa2::Int)
    _dry3!(M, M.vx[v], M.vy[v], M.vz[v], start, mode, splitSeg, fa1, fa2) != 0 && return 1
    _commit3!(M, v)
    return 0
end

# Start the mesher with a bounding tetrahedron around the box. The box
# has the corner lo and the corner hi.
function _start_bounding_tetra3!(M::NamedTuple, lox, loy, loz, hix, hiy, hiz)
    cx = 0.5 * (lox + hix); cy = 0.5 * (loy + hiy); cz = 0.5 * (loz + hiz)
    diag = sqrt((hix - lox)^2 + (hiy - loy)^2 + (hiz - loz)^2)
    s = 600.0 * diag
    v1 = _add_vertex3!(M, cx + s, cy + s, cz + s, 0, 0)
    v2 = _add_vertex3!(M, cx + s, cy - s, cz - s, 0, 0)
    v3 = _add_vertex3!(M, cx - s, cy + s, cz - s, 0, 0)
    v4 = _add_vertex3!(M, cx - s, cy - s, cz + s, 0, 0)
    c1 = _alloc_cell3!(M)
    M.cellV[][c1, 1] = v1; M.cellV[][c1, 2] = v2; M.cellV[][c1, 3] = v3; M.cellV[][c1, 4] = v4
    if _o3v3(M, v1, v2, v3, M.vx[v4], M.vy[v4], M.vz[v4]) < 0
        M.cellV[][c1, 3] = v4
        M.cellV[][c1, 4] = v3
    end
    for v in 1:4
        M.vCell[v] = c1
    end
    return c1
end

# ── The facet layer ──────────────────────────────────────────────────────
#
# Each facet has its own triangulation. A triangle of that triangulation
# is a subfacet. The layer projects a facet to 2-D by dropping one axis.
# The axis facDrop[f] is the axis with the largest normal component.
# The drop values 1, 2, and 3 keep (y, z), (x, z), and (x, y).
# The subfacets list their vertices in counterclockwise order in the
# projection. Local edge i of a subfacet joins local vertices i+1 and
# i+2. The neighbor sfN[k, i] lies across local edge i. The value 0
# means a facet boundary. The boundary edges are segments.
# The layer never reuses a subfacet number. A dead subfacet keeps its
# number, so a stale number on a cell face stays harmless.

# Return the first projected coordinate of the point (x, y, z) in facet f.
@inline function _proj_u3(M::NamedTuple, f::Int, x::Float64, y::Float64, z::Float64)
    d = M.facDrop[f]
    return d == 1 ? y : x
end

# Return the second projected coordinate of the point (x, y, z) in facet f.
@inline function _proj_w3(M::NamedTuple, f::Int, x::Float64, y::Float64, z::Float64)
    d = M.facDrop[f]
    return d == 3 ? y : z
end

# Return the sign of the 2-D orientation of the projected vertices a, b and
# the projected point (pu, pw) in facet f.
@inline function _o2f3(M::NamedTuple, f::Int, a::Int, b::Int, pu::Float64, pw::Float64)
    return _orient2d(_proj_u3(M, f, M.vx[a], M.vy[a], M.vz[a]), _proj_w3(M, f, M.vx[a], M.vy[a], M.vz[a]),
                     _proj_u3(M, f, M.vx[b], M.vy[b], M.vz[b]), _proj_w3(M, f, M.vx[b], M.vy[b], M.vz[b]),
                     pu, pw)
end

# Return the sign of the 2-D incircle test of subfacet k and the projected
# point (pu, pw).
@inline function _incircle_sf3(M::NamedTuple, k::Int, pu::Float64, pw::Float64)
    f = M.sfFac[k]
    a = M.sfV[][k, 1]; b = M.sfV[][k, 2]; c = M.sfV[][k, 3]
    return _incircle(_proj_u3(M, f, M.vx[a], M.vy[a], M.vz[a]), _proj_w3(M, f, M.vx[a], M.vy[a], M.vz[a]),
                     _proj_u3(M, f, M.vx[b], M.vy[b], M.vz[b]), _proj_w3(M, f, M.vx[b], M.vy[b], M.vz[b]),
                     _proj_u3(M, f, M.vx[c], M.vy[c], M.vz[c]), _proj_w3(M, f, M.vx[c], M.vy[c], M.vz[c]),
                     pu, pw)
end

# Return true when the projected point (pu, pw) is very near the line
# through the vertices a and b of facet f. The distance limit is 1e-10
# of the edge length.
# Two halves of a split segment are collinear in exact arithmetic. A
# rounding error can still put a new point on the wrong side of one
# half. This test treats such a point as lying on the boundary.
function _near_line3(M::NamedTuple, f::Int, a::Int, b::Int, pu::Float64, pw::Float64)
    au = _proj_u3(M, f, M.vx[a], M.vy[a], M.vz[a]); aw = _proj_w3(M, f, M.vx[a], M.vy[a], M.vz[a])
    bu = _proj_u3(M, f, M.vx[b], M.vy[b], M.vz[b]); bw = _proj_w3(M, f, M.vx[b], M.vy[b], M.vz[b])
    len = hypot(bu - au, bw - aw)
    cross = (bu - au) * (pw - aw) - (bw - aw) * (pu - au)
    return abs(cross) <= 1.0e-10 * len * len
end

# Add one facet. Return its number.
function _add_facet3!(M::NamedTuple, label::Int, drop::Int)
    push!(M.facLabel, label)
    push!(M.facDrop, drop)
    push!(M.facMate, 0)
    push!(M.facMateAxis, 0)
    push!(M.facMateCoord, 0.0)
    push!(M.facSF, 0)
    return length(M.facLabel)
end

# Add one subfacet with the vertices a, b, c to facet f. The function
# reorders the vertices so that they run counterclockwise in the
# projection. Return the subfacet number.
function _new_sf3!(M::NamedTuple, a::Int, b::Int, c::Int, f::Int)
    M.nSF[] += 1
    if M.nSF[] > M.capSF[]
        newCap = 2 * M.capSF[]
        M.sfV[] = _grow_matrix(M.sfV[], newCap)
        M.sfN[] = _grow_matrix(M.sfN[], newCap)
        resize!(M.sfFac, newCap); resize!(M.sfAlive, newCap); resize!(M.sfMark, newCap)
        M.capSF[] = newCap
    end
    k = M.nSF[]
    pu = _proj_u3(M, f, M.vx[c], M.vy[c], M.vz[c])
    pw = _proj_w3(M, f, M.vx[c], M.vy[c], M.vz[c])
    if _o2f3(M, f, a, b, pu, pw) < 0
        b, c = c, b
    end
    M.sfV[][k, 1] = a; M.sfV[][k, 2] = b; M.sfV[][k, 3] = c
    M.sfN[][k, 1] = 0; M.sfN[][k, 2] = 0; M.sfN[][k, 3] = 0
    M.sfFac[k] = f
    M.sfAlive[k] = 1
    M.sfMark[k] = 0
    return k
end

# Add one segment from vertex a to vertex b. The facets fa and fb share it.
# For an axis-aligned segment the function stores the vertex with the
# smaller coordinate first. Return the segment number.
function _add_segment3!(M::NamedTuple, a::Int, b::Int, fa::Int, fb::Int)
    axis = 0
    dx = M.vx[a] != M.vx[b]
    dy = M.vy[a] != M.vy[b]
    dz = M.vz[a] != M.vz[b]
    if dx && !dy && !dz
        axis = 1
        M.vx[a] > M.vx[b] && ((a, b) = (b, a))
    elseif dy && !dx && !dz
        axis = 2
        M.vy[a] > M.vy[b] && ((a, b) = (b, a))
    elseif dz && !dx && !dy
        axis = 3
        M.vz[a] > M.vz[b] && ((a, b) = (b, a))
    end
    M.nS[] += 1
    if M.nS[] > M.capS[]
        newCap = 2 * M.capS[]
        M.segV[] = _grow_matrix(M.segV[], newCap)
        M.segFac[] = _grow_matrix(M.segFac[], newCap)
        resize!(M.segAxis, newCap); resize!(M.segGroup, newCap); resize!(M.segAlive, newCap)
        M.capS[] = newCap
    end
    s = M.nS[]
    M.segV[][s, 1] = a; M.segV[][s, 2] = b
    M.segFac[][s, 1] = fa; M.segFac[][s, 2] = fb
    M.segAxis[s] = axis
    M.segGroup[s] = 0
    M.segAlive[s] = 1
    M.segOf[a < b ? (a, b) : (b, a)] = s
    M.vOnSeg[a] = 1
    M.vOnSeg[b] = 1
    return s
end

# Build the subfacets of facet f from a list of vertex triples. Link the
# subfacets across shared edges. Make a segment for each boundary edge,
# or add the facet to the segment that exists.
function _facet_build3!(M::NamedTuple, f::Int, tris::Vector{NTuple{3,Int}})
    first = 0
    edgeOwner = Dict{Tuple{Int,Int},Tuple{Int,Int}}()
    ids = Int[]
    for (a, b, c) in tris
        k = _new_sf3!(M, a, b, c, f)
        first == 0 && (first = k)
        push!(ids, k)
    end
    for k in ids
        for i in 1:3
            u = M.sfV[][k, i % 3 + 1]
            w = M.sfV[][k, (i + 1) % 3 + 1]
            key = u < w ? (u, w) : (w, u)
            if haskey(edgeOwner, key)
                k2, i2 = edgeOwner[key]
                M.sfN[][k, i] = k2
                M.sfN[][k2, i2] = k
            else
                edgeOwner[key] = (k, i)
            end
        end
    end
    M.facSF[f] = first
    for k in ids
        for i in 1:3
            M.sfN[][k, i] == 0 || continue
            u = M.sfV[][k, i % 3 + 1]
            w = M.sfV[][k, (i + 1) % 3 + 1]
            key = u < w ? (u, w) : (w, u)
            s = get(M.segOf, key, 0)
            if s == 0
                _add_segment3!(M, u, w, f, 0)
            else
                M.segFac[][s, 2] == 0 || error("internal error: a segment has more than two facets")
                M.segFac[][s, 2] = f
            end
        end
    end
    return nothing
end

# Walk inside facet f from subfacet start to the subfacet that holds the
# projected point (pu, pw). The walk crosses an edge when the point lies
# beyond it. A boundary edge with the end vertices hintA and hintB does
# not stop the walk as an error. The point lies on that edge.
# Return a code and a number.
#   1, k   the point lies in subfacet k or on its boundary
#   2, s   the walk reached a boundary edge of segment s
#   3, v   the point equals a vertex v
function _facet_locate3(M::NamedTuple, f::Int, pu::Float64, pw::Float64, start::Int, hintA::Int, hintB::Int)
    k = start
    steps = 0
    limit = 10 * M.nSF[] + 1000
    while true
        steps += 1
        steps > limit && error("internal error: facet location does not end")
        moved = false
        for j in 0:2
            i = mod(steps + j, 3) + 1
            u = M.sfV[][k, i % 3 + 1]
            w = M.sfV[][k, (i + 1) % 3 + 1]
            if _o2f3(M, f, u, w, pu, pw) < 0
                n = M.sfN[][k, i]
                if n == 0
                    if (u == hintA && w == hintB) || (u == hintB && w == hintA) || _near_line3(M, f, u, w, pu, pw)
                        continue
                    end
                    key = u < w ? (u, w) : (w, u)
                    return 2, get(M.segOf, key, 0)
                end
                k = n
                moved = true
                break
            end
        end
        if !moved
            for i in 1:3
                v = M.sfV[][k, i]
                if _proj_u3(M, f, M.vx[v], M.vy[v], M.vz[v]) == pu && _proj_w3(M, f, M.vx[v], M.vy[v], M.vz[v]) == pw
                    return 3, v
                end
            end
            return 1, k
        end
    end
end

# Insert the vertex v into facet f. The subfacet k0 holds the point or is
# next to it. A boundary edge with the end vertices hintA and hintB holds
# the point exactly. The function never makes a flat subfacet on a
# boundary edge that holds the point. The new subfacets are in M.newSF.
function _facet_insert3!(M::NamedTuple, f::Int, v::Int, k0::Int, hintA::Int, hintB::Int)
    pu = _proj_u3(M, f, M.vx[v], M.vy[v], M.vz[v])
    pw = _proj_w3(M, f, M.vx[v], M.vy[v], M.vz[v])
    empty!(M.cavSF)
    empty!(M.stack)
    M.stamp[] += 1
    st = M.stamp[]
    push!(M.cavSF, k0)
    push!(M.stack, k0)
    M.sfMark[k0] = st
    while !isempty(M.stack)
        k = pop!(M.stack)
        for i in 1:3
            n = M.sfN[][k, i]
            n == 0 && continue
            mk = M.sfMark[n]
            (mk == st || mk == -st) && continue
            if _incircle_sf3(M, n, pu, pw) > 0
                M.sfMark[n] = st
                push!(M.cavSF, n)
                push!(M.stack, n)
            else
                M.sfMark[n] = -st
            end
        end
    end
    empty!(M.newSF)
    for k in M.cavSF
        for i in 1:3
            n = M.sfN[][k, i]
            if n != 0 && M.sfMark[n] == st
                continue
            end
            a = M.sfV[][k, i % 3 + 1]
            b = M.sfV[][k, (i + 1) % 3 + 1]
            if n == 0
                onHint = (a == hintA && b == hintB) || (a == hintB && b == hintA)
                if onHint || _o2f3(M, f, a, b, pu, pw) == 0 || _near_line3(M, f, a, b, pu, pw)
                    continue
                end
            end
            nn = _new_sf3!(M, v, a, b, f)
            M.sfN[][nn, 1] = n
            if n != 0
                for j in 1:3
                    M.sfN[][n, j] == k && (M.sfN[][n, j] = nn)
                end
            end
            M.tmpSF[a] = nn
            push!(M.newSF, nn)
        end
    end
    for nn in M.newSF
        b = M.sfV[][nn, 3]
        t = M.tmpSF[b]
        if t != 0
            M.sfN[][nn, 2] = t
            M.sfN[][t, 3] = nn
        end
    end
    for nn in M.newSF
        M.tmpSF[M.sfV[][nn, 2]] = 0
    end
    for k in M.cavSF
        M.sfAlive[k] = 0
    end
    isempty(M.newSF) && error("internal error: a facet insertion made no subfacet")
    M.facSF[f] = M.newSF[1]
    return nothing
end

# Set the subfacet number on the faces of the cells M.newCells. A new
# subfacet that is not a face of the 3-D triangulation goes to queue B as
# a forced entry. Call the function after a facet insertion.
function _link_subfacets3!(M::NamedTuple)
    for nn in M.newSF
        a = M.sfV[][nn, 1]; b = M.sfV[][nn, 2]; d = M.sfV[][nn, 3]
        found = false
        for c in M.newCells
            ha = 0; hb = 0; hd = 0
            for i in 1:4
                v = M.cellV[][c, i]
                v == a && (ha = i)
                v == b && (hb = i)
                v == d && (hd = i)
            end
            (ha != 0 && hb != 0 && hd != 0) || continue
            li = 10 - ha - hb - hd
            M.cellF[][c, li] = nn
            n = M.cellN[][c, li]
            n != 0 && (M.cellF[][n, _slot_of3(M, n, c)] = nn)
            found = true
            break
        end
        found || _qpush!(M.qB, 2 * nn + 1)
    end
    return nothing
end

# ── Input geometry ───────────────────────────────────────────────────────

# Return the two axes that lie in a face with the normal axis a.
@inline _plane_axes(a::Int) = a == 1 ? (2, 3) : (a == 2 ? (1, 3) : (1, 2))

# Build the vertices, the facets, and the segments of an axis-aligned
# box. The box has the corner lo and n cells along each axis. A cell has
# the width delta[a] along the axis a. The numbers lo and delta must be
# multiples of one power of two, so that every lattice coordinate is
# exact. The array labels holds the label of each face in the order xmin,
# xmax, ymin, ymax, zmin, zmax. A vertex gets the lowest label of the
# faces that it lies on. The function does not insert the vertices into
# the 3-D triangulation. It returns the opposite corner of the box.
function _build_box_plc3!(M::NamedTuple, lo::NTuple{3,Float64}, delta::NTuple{3,Float64},
                         n::NTuple{3,Int}, labels::NTuple{6,Int})
    hi = (lo[1] + n[1] * delta[1], lo[2] + n[2] * delta[2], lo[3] + n[3] * delta[3])
    s = [collect(0.0:Float64(n[a])) for a in 1:3]
    lattice = Dict{NTuple{3,Int},Int}()
    # The label of a boundary lattice point is the lowest label of its faces.
    function label_of(idx::NTuple{3,Int})
        best = typemax(Int)
        for a in 1:3
            if idx[a] == 0
                best = min(best, labels[2 * (a - 1) + 1])
            end
            if idx[a] == n[a]
                best = min(best, labels[2 * (a - 1) + 2])
            end
        end
        return best
    end
    function add_point!(idx::NTuple{3,Int}, x::Float64, y::Float64, z::Float64)
        lattice[idx] = _add_vertex3!(M, x, y, z, label_of(idx), 0)
    end
    # Corners.
    for iz in 0:1, iy in 0:1, ix in 0:1
        idx = (ix * n[1], iy * n[2], iz * n[3])
        add_point!(idx, ix == 0 ? lo[1] : hi[1], iy == 0 ? lo[2] : hi[2], iz == 0 ? lo[3] : hi[3])
    end
    # Points inside the 12 edges. Point k is lo + k * delta.
    for a in 1:3
        b, c = _plane_axes(a)
        for p in 0:1, q in 0:1
            n[a] >= 2 || continue
            cnt = n[a] - 1
            xs = zeros(Float64, cnt); ys = zeros(Float64, cnt); zs = zeros(Float64, cnt)
            p0 = [0.0, 0.0, 0.0]
            p0[a] = lo[a]
            p0[b] = p == 0 ? lo[b] : hi[b]
            p0[c] = q == 0 ? lo[c] : hi[c]
            p1 = copy(p0)
            p1[a] = lo[a] + delta[a]
            line_points_kernel(p0[1], p0[2], p0[3], p1[1], p1[2], p1[3], s[a], 2, n[a], 0, xs, ys, zs)
            for k in 1:cnt
                idx = [0, 0, 0]
                idx[a] = k; idx[b] = p * n[b]; idx[c] = q * n[c]
                add_point!((idx[1], idx[2], idx[3]), xs[k], ys[k], zs[k])
            end
        end
    end
    # Interior grid points of the 6 faces, and the facets.
    for a in 1:3, side in 0:1
        b, c = _plane_axes(a)
        n1 = n[b]; n2 = n[c]
        if n1 >= 2 && n2 >= 2
            cnt = (n1 - 1) * (n2 - 1)
            xs = zeros(Float64, cnt); ys = zeros(Float64, cnt); zs = zeros(Float64, cnt)
            o = [lo[1], lo[2], lo[3]]
            o[a] = side == 0 ? lo[a] : hi[a]
            u = [0.0, 0.0, 0.0]; w = [0.0, 0.0, 0.0]
            u[b] = delta[b]; w[c] = delta[c]
            rect_grid_points_kernel(o[1], o[2], o[3], u[1], u[2], u[3], w[1], w[2], w[3], s[b], s[c], n1, n2, 0, xs, ys, zs)
            for jl in 1:n2 - 1, il in 1:n1 - 1
                idx = [0, 0, 0]
                idx[a] = side * n[a]; idx[b] = il; idx[c] = jl
                k = (jl - 1) * (n1 - 1) + il
                add_point!((idx[1], idx[2], idx[3]), xs[k], ys[k], zs[k])
            end
        end
        fid = _add_facet3!(M, labels[2 * (a - 1) + side + 1], a)
        tris = NTuple{3,Int}[]
        vid(i, j) = begin
            idx = [0, 0, 0]
            idx[a] = side * n[a]; idx[b] = i; idx[c] = j
            lattice[(idx[1], idx[2], idx[3])]
        end
        for j in 0:n2 - 1, i in 0:n1 - 1
            push!(tris, (vid(i, j), vid(i + 1, j), vid(i + 1, j + 1)))
            push!(tris, (vid(i, j), vid(i + 1, j + 1), vid(i, j + 1)))
        end
        _facet_build3!(M, fid, tris)
    end
    return hi
end

# Build an icosphere with the given number of subdivisions. Return the
# direction vectors of the vertices (n-by-3) and the triangles.
function _icosphere3(level::Int)
    t = (1.0 + sqrt(5.0)) / 2.0
    dirs = Vector{NTuple{3,Float64}}()
    for (x, y, z) in ((-1.0, t, 0.0), (1.0, t, 0.0), (-1.0, -t, 0.0), (1.0, -t, 0.0),
                      (0.0, -1.0, t), (0.0, 1.0, t), (0.0, -1.0, -t), (0.0, 1.0, -t),
                      (t, 0.0, -1.0), (t, 0.0, 1.0), (-t, 0.0, -1.0), (-t, 0.0, 1.0))
        len = sqrt(x * x + y * y + z * z)
        push!(dirs, (x / len, y / len, z / len))
    end
    tris = NTuple{3,Int}[(1, 12, 6), (1, 6, 2), (1, 2, 8), (1, 8, 11), (1, 11, 12),
                         (2, 6, 10), (6, 12, 5), (12, 11, 3), (11, 8, 7), (8, 2, 9),
                         (4, 10, 5), (4, 5, 3), (4, 3, 7), (4, 7, 9), (4, 9, 10),
                         (5, 10, 6), (3, 5, 12), (7, 3, 11), (9, 7, 8), (10, 9, 2)]
    for _ in 1:level
        cache = Dict{Tuple{Int,Int},Int}()
        function midpoint(i::Int, j::Int)
            key = i < j ? (i, j) : (j, i)
            haskey(cache, key) && return cache[key]
            p = dirs[i]; q = dirs[j]
            mx = p[1] + q[1]; my = p[2] + q[2]; mz = p[3] + q[3]
            len = sqrt(mx * mx + my * my + mz * mz)
            push!(dirs, (mx / len, my / len, mz / len))
            cache[key] = length(dirs)
            return length(dirs)
        end
        newTris = NTuple{3,Int}[]
        for (a, b, c) in tris
            ab = midpoint(a, b); bc = midpoint(b, c); ca = midpoint(c, a)
            push!(newTris, (a, ab, ca)); push!(newTris, (b, bc, ab))
            push!(newTris, (c, ca, bc)); push!(newTris, (ab, bc, ca))
        end
        tris = newTris
    end
    out = zeros(Float64, length(dirs), 3)
    for (i, d) in enumerate(dirs)
        out[i, 1] = d[1]; out[i, 2] = d[2]; out[i, 3] = d[3]
    end
    return out, tris
end

# Return the number of icosphere subdivisions for a sphere of radius r and
# the target edge length h. The edge of an icosphere with L subdivisions
# is about 1.0515 * r / 2^L.
function _sphere_level3(r::Float64, h::Float64)
    level = 1
    while level < 6 && 1.0515 * r / 2^level > h
        level += 1
    end
    return level
end

# Build the vertices, the facets, and the segments of a sphere as an
# icosphere. Each triangle is one facet with the given label.
function _build_sphere_plc3!(M::NamedTuple, cx::Float64, cy::Float64, cz::Float64, r::Float64, level::Int, label::Int)
    dirs, tris = _icosphere3(level)
    n = size(dirs, 1)
    # The icosphere is exactly symmetric. Four of its vertices can lie on
    # one circle in one plane. A rounding error on a tilted facet then
    # makes a flat cell between the two ways to triangulate the circle.
    # A small jitter of the directions removes such circles. The
    # jitter is 0.2 percent of an edge. It is the same for every run.
    rng = MersenneTwister(1000 + label)
    jitter = 0.002 * 1.0515 / 2^level
    rx = dirs[:, 1] .+ jitter .* (2.0 .* rand(rng, n) .- 1.0)
    ry = dirs[:, 2] .+ jitter .* (2.0 .* rand(rng, n) .- 1.0)
    rz = dirs[:, 3] .+ jitter .* (2.0 .* rand(rng, n) .- 1.0)
    xs = zeros(Float64, n); ys = zeros(Float64, n); zs = zeros(Float64, n)
    sphere_project_kernel(cx, cy, cz, r, n, rx, ry, rz, 0, xs, ys, zs)
    ids = zeros(Int, n)
    for k in 1:n
        qx, qy, qz = _quantize_point3(M, xs[k], ys[k], zs[k])
        ids[k] = _add_vertex3!(M, qx, qy, qz, label, 0)
    end
    for (a, b, c) in tris
        va = ids[a]; vb = ids[b]; vc = ids[c]
        ux = M.vx[vb] - M.vx[va]; uy = M.vy[vb] - M.vy[va]; uz = M.vz[vb] - M.vz[va]
        wx = M.vx[vc] - M.vx[va]; wy = M.vy[vc] - M.vy[va]; wz = M.vz[vc] - M.vz[va]
        nx = abs(uy * wz - uz * wy); ny = abs(uz * wx - ux * wz); nz = abs(ux * wy - uy * wx)
        drop = (nx >= ny && nx >= nz) ? 1 : (ny >= nz ? 2 : 3)
        fid = _add_facet3!(M, label, drop)
        _facet_build3!(M, fid, [(va, vb, vc)])
    end
    return nothing
end

# Link the opposite faces and the box edges of the outer box that mirror
# each other. The outer box facets have the numbers 1 to 6. Each
# periodic axis links the two faces that are normal to it. Segments on
# the box edges that map onto each other form a group.
function _setup_periodic3!(M::NamedTuple, lo::NTuple{3,Float64}, hi::NTuple{3,Float64}, per::NTuple{3,Bool})
    for a in 1:3
        per[a] || continue
        f0 = 2 * (a - 1) + 1
        f1 = f0 + 1
        M.facMate[f0] = f1; M.facMate[f1] = f0
        M.facMateAxis[f0] = a; M.facMateAxis[f1] = a
        M.facMateCoord[f0] = hi[a]; M.facMateCoord[f1] = lo[a]
    end
    keyOf = Dict{NTuple{5,Float64},Int}()
    empty!(M.groups)
    for s in 1:M.nS[]
        ax = M.segAxis[s]
        ax == 0 && continue
        (M.segFac[][s, 1] in 1:6 && M.segFac[][s, 2] in 1:6) || continue
        a = M.segV[][s, 1]; b = M.segV[][s, 2]
        pa = (M.vx[a], M.vy[a], M.vz[a])
        pb = (M.vx[b], M.vy[b], M.vz[b])
        b1, c1 = _plane_axes(ax)
        canon(d, val) = (per[d] && val == hi[d]) ? lo[d] : val
        key = (Float64(ax), pa[ax], pb[ax], canon(b1, pa[b1]), canon(c1, pa[c1]))
        g = get(keyOf, key, 0)
        if g == 0
            push!(M.groups, Int[])
            g = length(M.groups)
            keyOf[key] = g
        end
        push!(M.groups[g], s)
        M.segGroup[s] = g
    end
    return nothing
end

# ── Scanning for work ────────────────────────────────────────────────────

# Examine the cells ids[1:n]. Queue each segment that an apex vertex
# encroaches on, and each subfacet that an apex vertex encroaches on.
# When M.doBad is true, queue each cell in the domain that is too large
# or has a shape that is too poor.
function _scan3!(M::NamedTuple, ids::Vector{Int}, n::Int)
    nPair = 0
    for k in 1:n
        c = ids[k]
        for (la, lb) in _tet_local_edges()
            s = _seg_of_edge3(M, M.cellV[][c, la], M.cellV[][c, lb])
            s == 0 && continue
            for l in 1:4
                (l == la || l == lb) && continue
                nPair += 1
                _ws3!(M, nPair)
                t = M.cellV[][c, l]
                M.wSeg[nPair] = s
                M.wTx[nPair] = M.vx[t]; M.wTy[nPair] = M.vy[t]; M.wTz[nPair] = M.vz[t]
            end
        end
    end
    if nPair > 0
        segment_encroach_3d_kernel(M.segV[], M.vx, M.vy, M.vz, M.wSeg, M.wTx, M.wTy, M.wTz, nPair, M.wFlag)
        for k in 1:nPair
            M.wFlag[k] == 1 && _qpush!(M.qA, 2 * M.wSeg[k])
        end
    end
    nPair = 0
    for k in 1:n
        c = ids[k]
        for i in 1:4
            s = _live_sf3(M, c, i)
            s == 0 && continue
            nPair += 1
            _ws3!(M, nPair)
            t = M.cellV[][c, i]
            M.wSeg[nPair] = s
            M.wTx[nPair] = M.vx[t]; M.wTy[nPair] = M.vy[t]; M.wTz[nPair] = M.vz[t]
        end
    end
    if nPair > 0
        subfacet_encroach_kernel(M.sfV[], M.vx, M.vy, M.vz, M.wSeg, M.wTx, M.wTy, M.wTz, nPair, M.wFlag)
        for k in 1:nPair
            M.wFlag[k] == 1 && _qpush!(M.qB, 2 * M.wSeg[k])
        end
    end
    if M.doBad[]
        nIn = 0
        for k in 1:n
            c = ids[k]
            if M.cellAlive[c] == 1 && M.cellIn[c] == 1
                nIn += 1
                _ws3!(M, nIn)
                M.wCell[nIn] = c
            end
        end
        nIn > 0 && _bad_cells3!(M, nIn, true)
    end
    return nothing
end

# Compute the quality of the cells M.wCell[1:n]. Queue the bad cells
# when push is true. Return the number of bad cells.
function _bad_cells3!(M::NamedTuple, n::Int, push::Bool)
    cell_circumcenter_3d_kernel(M.cellV[], M.vx, M.vy, M.vz, M.wCell, n, M.ccx, M.ccy, M.ccz, M.ccR)
    size_field_3d_kernel(M.ccx, M.ccy, M.ccz, n, M.obsKind[], M.obsCx[], M.obsCy[], M.obsCz[],
                         M.obsHx[], M.obsHy[], M.obsHz[], M.nObs[], M.hMax[], M.hObs[], M.grading[], M.hAt)
    cell_quality_3d_kernel(M.cellV[], M.vx, M.vy, M.vz, M.wCell, n, M.ccR, M.hAt, M.ratio, M.sizeRatio)
    nBad = 0
    for k in 1:n
        if M.ratio[k] > M.ratioBound[] || M.sizeRatio[k] > 1.0
            nBad += 1
            push && _qpush!(M.qC, M.wCell[k])
        end
    end
    return nBad
end

# Return true when a vertex encroaches on segment s. Test the vertices
# of the cells around the segment.
function _seg_encroached3(M::NamedTuple, s::Int)
    a = M.segV[][s, 1]
    b = M.segV[][s, 2]
    _star3!(M, a)
    nPair = 0
    for c in M.star
        hasB = false
        for i in 1:4
            M.cellV[][c, i] == b && (hasB = true)
        end
        hasB || continue
        for i in 1:4
            t = M.cellV[][c, i]
            (t == a || t == b) && continue
            nPair += 1
            _ws3!(M, nPair)
            M.wSeg[nPair] = s
            M.wTx[nPair] = M.vx[t]; M.wTy[nPair] = M.vy[t]; M.wTz[nPair] = M.vz[t]
        end
    end
    nPair == 0 && return false
    segment_encroach_3d_kernel(M.segV[], M.vx, M.vy, M.vz, M.wSeg, M.wTx, M.wTy, M.wTz, nPair, M.wFlag)
    for k in 1:nPair
        M.wFlag[k] == 1 && return true
    end
    return false
end

# Return true when a vertex encroaches on subfacet k. The cell c holds
# the face i of the subfacet. Test the apex vertex of each side.
function _sf_encroached3(M::NamedTuple, k::Int, c::Int, i::Int)
    _ws3!(M, 2)
    t = M.cellV[][c, i]
    M.wSeg[1] = k
    M.wTx[1] = M.vx[t]; M.wTy[1] = M.vy[t]; M.wTz[1] = M.vz[t]
    nPair = 1
    n = M.cellN[][c, i]
    if n != 0
        t2 = M.cellV[][n, _slot_of3(M, n, c)]
        nPair = 2
        M.wSeg[2] = k
        M.wTx[2] = M.vx[t2]; M.wTy[2] = M.vy[t2]; M.wTz[2] = M.vz[t2]
    end
    subfacet_encroach_kernel(M.sfV[], M.vx, M.vy, M.vz, M.wSeg, M.wTx, M.wTy, M.wTz, nPair, M.wFlag)
    return M.wFlag[1] == 1 || (nPair == 2 && M.wFlag[2] == 1)
end

# Queue the segments and subfacets that the last insertion deleted. A
# forced entry always splits. A subfacet that is still alive after the
# facet insertion is missing from the 3-D triangulation.
function _queue_lost3!(M::NamedTuple)
    for t in unique(M.lostSeg)
        M.segAlive[t] == 1 && _qpush!(M.qA, 2 * t + 1)
    end
    for k in unique(M.lostSF)
        M.sfAlive[k] == 1 && _qpush!(M.qB, 2 * k + 1)
    end
    empty!(M.lostSeg)
    empty!(M.lostSF)
    return nothing
end

# ── Splitting segments ───────────────────────────────────────────────────

# Return the midpoint of segment s.
function _segment_mid3(M::NamedTuple, s::Int)
    _ws3!(M, 1)
    M.wSeg[1] = s
    segment_midpoint_kernel(M.segV[], M.vx, M.vy, M.vz, M.wSeg, 1, M.wTx, M.wTy, M.wTz)
    return M.wTx[1], M.wTy[1], M.wTz[1]
end

# Return the first cell near the segment s for a point walk.
function _near_segment3(M::NamedTuple, s::Int)
    a = M.segV[][s, 1]
    c = _find_edge3(M, a, M.segV[][s, 2])
    return c != 0 ? c : M.vCell[a]
end

# Insert the point (px, py, pz) on segment s. The first half keeps the
# number s. The second half gets a new number. Return that number. The
# caller must have tested the insertion with _dry3! first.
function _split_segment_apply3!(M::NamedTuple, s::Int, px::Float64, py::Float64, pz::Float64)
    a = M.segV[][s, 1]
    b = M.segV[][s, 2]
    f1 = M.segFac[][s, 1]
    f2 = M.segFac[][s, 2]
    label = typemax(Int)
    f1 != 0 && (label = min(label, M.facLabel[f1]))
    f2 != 0 && (label = min(label, M.facLabel[f2]))
    start = _near_segment3(M, s)
    _dry3!(M, px, py, pz, start, 2, s, f1, f2) == 0 ||
        error("internal error: a segment split failed after a successful test")
    v = _add_vertex3!(M, px, py, pz, label, 1)
    _commit3!(M, v)
    sB = _add_segment3!(M, v, b, f1, f2)
    delete!(M.segOf, a < b ? (a, b) : (b, a))
    M.segV[][s, 2] = v
    M.segOf[a < v ? (a, v) : (v, a)] = s
    for f in (f1, f2)
        f == 0 && continue
        pu = _proj_u3(M, f, px, py, pz)
        pw = _proj_w3(M, f, px, py, pz)
        code, id = _facet_locate3(M, f, pu, pw, M.facSF[f], a, b)
        code == 1 || error("internal error: a segment split point is not in its facet (locate code $code, id $id, facet $f)")
        _facet_insert3!(M, f, v, id, a, b)
        _link_subfacets3!(M)
    end
    _queue_lost3!(M)
    _scan3!(M, M.newCells, length(M.newCells))
    return sB
end

# Process one entry of queue A. The entry is 2*s + force. A forced entry
# always splits. An entry that is not forced splits only when a vertex
# still encroaches on the segment, or when the segment is not an edge.
# The members of the periodic group of the segment split together. A
# split is never rejected. This gives queue A the priority that makes the
# whole refinement end. Each split queues the items that it deletes.
function _process_segment3!(M::NamedTuple, item::Int)
    s = div(item, 2)
    force = mod(item, 2)
    M.segAlive[s] == 1 || return nothing
    present = _find_edge3(M, M.segV[][s, 1], M.segV[][s, 2]) != 0
    if present && force == 0 && !_seg_encroached3(M, s)
        return nothing
    end
    members = Int[s]
    g = M.segGroup[s]
    if g != 0
        for t in M.groups[g]
            (t != s && M.segAlive[t] == 1) && push!(members, t)
        end
    end
    px, py, pz = _segment_mid3(M, s)
    px, py, pz = _quantize_point3(M, px, py, pz)
    ax = M.segAxis[s]
    pts = NTuple{3,Float64}[(px, py, pz)]
    for t in members[2:end]
        a = M.segV[][t, 1]
        qx = ax == 1 ? px : M.vx[a]
        qy = ax == 2 ? py : M.vy[a]
        qz = ax == 3 ? pz : M.vz[a]
        push!(pts, _quantize_point3(M, qx, qy, qz))
    end
    seconds = Int[]
    for (idx, t) in enumerate(members)
        q = pts[idx]
        push!(seconds, _split_segment_apply3!(M, t, q[1], q[2], q[3]))
    end
    if g != 0
        push!(M.groups, seconds)
        ng = length(M.groups)
        for sb in seconds
            M.segGroup[sb] = ng
        end
    end
    return nothing
end

# ── Splitting subfacets ──────────────────────────────────────────────────

# Insert the point (px, py, pz) into facet f. The subfacet startSF holds
# the point or lies next to it. The caller must have tested the
# insertion with _dry3! first.
function _facet_vertex_apply3!(M::NamedTuple, f::Int, px::Float64, py::Float64, pz::Float64, startSF::Int)
    start = M.vCell[M.sfV[][startSF, 1]]
    _dry3!(M, px, py, pz, start, 3, 0, f, 0) == 0 ||
        error("internal error: a facet insertion failed after a successful test")
    v = _add_vertex3!(M, px, py, pz, M.facLabel[f], 0)
    _commit3!(M, v)
    _facet_insert3!(M, f, v, startSF, 0, 0)
    _link_subfacets3!(M)
    _queue_lost3!(M)
    _scan3!(M, M.newCells, length(M.newCells))
    return v
end

# Process one entry of queue B. The entry is 2*k + force. A forced entry
# always splits. An entry that is not forced splits only when a vertex
# still encroaches on the subfacet, or when the subfacet is not a face.
function _process_subfacet3!(M::NamedTuple, item::Int)
    k = div(item, 2)
    force = mod(item, 2)
    M.sfAlive[k] == 1 || return nothing
    a = M.sfV[][k, 1]; b = M.sfV[][k, 2]; d = M.sfV[][k, 3]
    c, i = _find_face3(M, a, b, d)
    if c != 0
        M.cellF[][c, i] = k
        n = M.cellN[][c, i]
        n != 0 && (M.cellF[][n, _slot_of3(M, n, c)] = k)
        if force == 0 && !_sf_encroached3(M, k, c, i)
            return nothing
        end
    end
    f = M.sfFac[k]
    _ws3!(M, 1)
    M.wSeg[1] = k
    subfacet_circumcenter_kernel(M.sfV[], M.vx, M.vy, M.vz, M.wSeg, 1, M.ccx, M.ccy, M.ccz, M.ccR)
    px = M.ccx[1]; py = M.ccy[1]; pz = M.ccz[1]
    isfinite(px + py + pz) || return nothing
    px, py, pz = _quantize_point3(M, px, py, pz)
    # Locate the point in the facet. A boundary edge on the way is a
    # segment that the point encroaches on.
    code, id = _facet_locate3(M, f, _proj_u3(M, f, px, py, pz), _proj_w3(M, f, px, py, pz), k, 0, 0)
    if code == 2
        id != 0 && _qpush!(M.qA, 2 * id + 1)
        _qpush!(M.qB, item)
        return nothing
    elseif code == 3
        return nothing
    end
    startSF = id
    f2 = M.facMate[f]
    qx = px; qy = py; qz = pz
    startSF2 = 0
    if f2 != 0
        ax = M.facMateAxis[f]
        qx = ax == 1 ? M.facMateCoord[f] : px
        qy = ax == 2 ? M.facMateCoord[f] : py
        qz = ax == 3 ? M.facMateCoord[f] : pz
        qx, qy, qz = _quantize_point3(M, qx, qy, qz)
        code2, id2 = _facet_locate3(M, f2, _proj_u3(M, f2, qx, qy, qz), _proj_w3(M, f2, qx, qy, qz), M.facSF[f2], 0, 0)
        if code2 == 2
            id2 != 0 && _qpush!(M.qA, 2 * id2 + 1)
            _qpush!(M.qB, item)
            return nothing
        elseif code2 == 3
            return nothing
        end
        startSF2 = id2
    end
    blockSeg = Int[]
    blockSF = Int[]
    reject = false
    if _dry3!(M, px, py, pz, M.vCell[M.sfV[][startSF, 1]], 3, 0, f, 0) != 0
        reject = true
        append!(blockSeg, M.violSeg); append!(blockSF, M.violSF)
    end
    if f2 != 0 && _dry3!(M, qx, qy, qz, M.vCell[M.sfV[][startSF2, 1]], 3, 0, f2, 0) != 0
        reject = true
        append!(blockSeg, M.violSeg); append!(blockSF, M.violSF)
    end
    if reject
        for t in unique(blockSeg)
            _qpush!(M.qA, 2 * t + 1)
        end
        for q in unique(blockSF)
            q != k && _qpush!(M.qB, 2 * q + 1)
        end
        _qpush!(M.qB, item)
        return nothing
    end
    _facet_vertex_apply3!(M, f, px, py, pz, startSF)
    f2 != 0 && _facet_vertex_apply3!(M, f2, qx, qy, qz, startSF2)
    return nothing
end

# ── Splitting cells ──────────────────────────────────────────────────────

# Process one cell of queue C. Insert its circumcenter when the cell is
# still bad. When the insertion is not allowed, queue the segments and
# subfacets that block it and put the cell back at the end of queue C.
function _process_cell3!(M::NamedTuple, c::Int)
    (M.cellAlive[c] == 1 && M.cellIn[c] == 1) || return nothing
    _ws3!(M, 1)
    M.wCell[1] = c
    _bad_cells3!(M, 1, false) == 0 && return nothing
    px = M.ccx[1]; py = M.ccy[1]; pz = M.ccz[1]
    isfinite(px + py + pz) || return nothing
    px, py, pz = _quantize_point3(M, px, py, pz)
    if _dry3!(M, px, py, pz, c, 1, 0, 0, 0) == 0
        v = _add_vertex3!(M, px, py, pz, 0, 0)
        _commit3!(M, v)
        _scan3!(M, M.newCells, length(M.newCells))
    else
        for t in unique(M.violSeg)
            _qpush!(M.qA, 2 * t + 1)
        end
        for q in unique(M.violSF)
            _qpush!(M.qB, 2 * q + 1)
        end
        _qpush!(M.qC, c)
    end
    return nothing
end

# ── Queues and checks ────────────────────────────────────────────────────

# Process queue A and queue B until both are empty. Queue A has priority.
function _run_boundary_queues3!(M::NamedTuple)
    iter = 0
    while true
        if !_qempty(M.qA)
            _process_segment3!(M, _qpop!(M.qA))
        elseif !_qempty(M.qB)
            _process_subfacet3!(M, _qpop!(M.qB))
        else
            break
        end
        iter += 1
        iter > M.iterLimit && error("recovery of the boundary does not end. Check the obstacles and mesh_size.")
    end
    return nothing
end

# Process all three queues until all are empty. A has priority over B,
# and B has priority over C.
function _run_refinement3!(M::NamedTuple)
    iter = 0
    while true
        if !_qempty(M.qA)
            _process_segment3!(M, _qpop!(M.qA))
        elseif !_qempty(M.qB)
            _process_subfacet3!(M, _qpop!(M.qB))
        elseif !_qempty(M.qC)
            _process_cell3!(M, _qpop!(M.qC))
        else
            break
        end
        iter += 1
        iter > M.iterLimit && error("refinement does not end. Check the obstacles, mesh_size, and radius_edge_bound.")
    end
    return nothing
end

# Check that each live segment is an edge and each live subfacet is a
# face. Set the subfacet numbers on the faces. Queue each item that is
# missing as a forced entry. Return the number of missing items.
function _mark_present3!(M::NamedTuple)
    nMissing = 0
    for k in 1:M.nSF[]
        M.sfAlive[k] == 1 || continue
        c, i = _find_face3(M, M.sfV[][k, 1], M.sfV[][k, 2], M.sfV[][k, 3])
        if c == 0
            _qpush!(M.qB, 2 * k + 1)
            nMissing += 1
        else
            M.cellF[][c, i] = k
            n = M.cellN[][c, i]
            n != 0 && (M.cellF[][n, _slot_of3(M, n, c)] = k)
        end
    end
    for s in 1:M.nS[]
        M.segAlive[s] == 1 || continue
        if _find_edge3(M, M.segV[][s, 1], M.segV[][s, 2]) == 0
            _qpush!(M.qA, 2 * s + 1)
            nMissing += 1
        end
    end
    return nMissing
end

# Make every segment an edge and every subfacet a face. Remove all
# encroachment. The domain flags are not used yet.
function _recover3!(M::NamedTuple)
    M.doBad[] = false
    rounds = 0
    while true
        rounds += 1
        rounds > 1000 && error("recovery of the boundary does not end")
        nMissing = _mark_present3!(M)
        ids = Int[c for c in 1:M.nC[] if M.cellAlive[c] == 1]
        _scan3!(M, ids, length(ids))
        queued = !_qempty(M.qA) || !_qempty(M.qB)
        _run_boundary_queues3!(M)
        (nMissing == 0 && !queued) && break
    end
    return nothing
end

# ── Domain flags ─────────────────────────────────────────────────────────

# Spread the flag 0 from the cells on the stack to their neighbors. The
# spread does not cross a face that carries a subfacet.
function _flood3!(M::NamedTuple, stack::Vector{Int})
    while !isempty(stack)
        c = pop!(stack)
        for i in 1:4
            n = M.cellN[][c, i]
            (n == 0 || _live_sf3(M, c, i) != 0) && continue
            if M.cellIn[n] == 1
                M.cellIn[n] = 0
                push!(stack, n)
            end
        end
    end
    return nothing
end

# Set cellIn to 1 for each cell in the domain and to 0 for each cell
# outside it. A cell is outside when the flood from the bounding
# tetrahedron or from an obstacle center reaches it.
function _classify3!(M::NamedTuple)
    for c in 1:M.nC[]
        M.cellIn[c] = M.cellAlive[c]
    end
    stack = Int[]
    for c in 1:M.nC[]
        M.cellAlive[c] == 1 || continue
        if M.cellV[][c, 1] <= 4 || M.cellV[][c, 2] <= 4 || M.cellV[][c, 3] <= 4 || M.cellV[][c, 4] <= 4
            M.cellIn[c] = 0
            push!(stack, c)
        end
    end
    _flood3!(M, stack)
    for k in 1:M.nObs[]
        code, id = _locate3(M, M.obsCx[][k], M.obsCy[][k], M.obsCz[][k], M.vCell[5], false)
        c = code == 3 ? M.vCell[id] : id
        if M.cellIn[c] == 1
            M.cellIn[c] = 0
            push!(stack, c)
            _flood3!(M, stack)
        end
    end
    return nothing
end

# ── Options and obstacles ────────────────────────────────────────────────

# Check the mesh_options record. Set the defaults. Return a NamedTuple.
function _resolve_mesh_options3(raw)
    raw isa AbstractDict || error("'mesh_options' must be a mapping")
    known = ("mesh_size", "mesh_size_obstacle", "grading", "radius_edge_bound", "smoothing_sweeps", "max_vertices", "min_dihedral", "max_dihedral")
    for k in keys(raw)
        k in known || error("mesh_options: unknown key '$k'. Valid keys: $(join(known, ", "))")
    end
    haskey(raw, "mesh_size") || error("mesh_options: 'mesh_size' is required")
    num(key, default) = begin
        v = get(raw, key, default)
        v isa Real || error("mesh_options['$key'] must be a number, got $(repr(v))")
        Float64(v)
    end
    meshSize = num("mesh_size", nothing)
    meshSize > 0 || error("mesh_options['mesh_size'] must be greater than 0")
    meshSizeObs = num("mesh_size_obstacle", meshSize)
    meshSizeObs > 0 || error("mesh_options['mesh_size_obstacle'] must be greater than 0")
    grading = num("grading", 0.3)
    grading >= 0 || error("mesh_options['grading'] must be 0 or more")
    bound = num("radius_edge_bound", 2.0)
    bound >= 2.0 || error("mesh_options['radius_edge_bound'] must be 2.0 or more")
    sweeps = num("smoothing_sweeps", 3)
    (sweeps >= 0 && sweeps == floor(sweeps)) || error("mesh_options['smoothing_sweeps'] must be a whole number, 0 or more")
    maxVertices = num("max_vertices", 1000000)
    (maxVertices >= 100 && maxVertices == floor(maxVertices)) || error("mesh_options['max_vertices'] must be a whole number, 100 or more")
    minDihedral = num("min_dihedral", 10.0)
    (minDihedral >= 0 && minDihedral <= 60) || error("mesh_options['min_dihedral'] must be in [0, 60]")
    maxDihedral = num("max_dihedral", 0.0)
    (maxDihedral == 0 || (maxDihedral >= 90 && maxDihedral < 180)) || error("mesh_options['max_dihedral'] must be 0 (off) or in [90, 180)")
    return (meshSize = meshSize, meshSizeObstacle = meshSizeObs, grading = grading,
            radiusEdgeBound = bound, smoothingSweeps = Int(sweeps), maxVertices = Int(maxVertices),
            minDihedral = minDihedral, maxDihedral = maxDihedral)
end

# Check each obstacle record. Raise an error for a bad field, for an
# obstacle that crosses or touches the box, and for two obstacles whose
# bounding boxes overlap. Warn when a gap is smaller than hObs.
function _validate_obstacles3(logBuf::IOBuffer, obstacles, box::Tuple, hObs::Float64)
    isempty(obstacles) && return
    lo = Vector{Vector{Float64}}()
    hi = Vector{Vector{Float64}}()
    for (i, obs) in enumerate(obstacles)
        idx = i - 1
        shape = get(obs, "shape", nothing)
        if !(shape in _obstacle_shapes_3d())
            error("obstacle $idx: shape must be one of $(_obstacle_shapes_3d()) for a 3-D mesh, got '$shape'")
        end
        c = get(obs, "center", nothing)
        (c === nothing || length(c) != 3) && error("obstacle $idx: 'center' must have 3 entries")
        c = Float64.(c)
        if shape == "sphere"
            r = get(obs, "radius", nothing)
            (r === nothing || r <= 0) && error("obstacle $idx: 'radius' must be a positive number")
            half = fill(Float64(r), 3)
        else
            he = get(obs, "half_extents", nothing)
            (he === nothing || length(he) != 3 || any(<=(0), he)) && error("obstacle $idx: 'half_extents' must have 3 positive entries")
            half = Float64.(he)
        end
        for d in 1:3
            if c[d] - half[d] <= 0 || c[d] + half[d] >= box[d]
                error("obstacle $idx: crosses or touches the box boundary along axis $(_dirs()[d])")
            end
            gap = min(c[d] - half[d], box[d] - c[d] - half[d])
            gap < hObs && _log(logBuf, "  [warn] obstacle $idx: the gap to the box along axis $(_dirs()[d]) is smaller than mesh_size_obstacle")
        end
        push!(lo, c .- half)
        push!(hi, c .+ half)
    end
    n = length(lo)
    for i in 1:n, j in (i + 1):n
        separation = maximum(max(lo[j][d] - hi[i][d], lo[i][d] - hi[j][d]) for d in 1:3)
        if separation <= 0
            error("obstacles $(i-1) and $(j-1): the bounding boxes overlap")
        elseif separation < hObs
            _log(logBuf, "  [warn] obstacles $(i-1) and $(j-1): the gap is smaller than mesh_size_obstacle")
        end
    end
end

# Build the record of one obstacle for the HDF5 file.
function _obstacle_record3(obs, marker::Int)
    center = Tuple(Float64.(obs["center"]))
    if obs["shape"] == "sphere"
        return (shape=obs["shape"], center=center, radius=Float64(obs["radius"]), marker=marker)
    else
        return (shape=obs["shape"], center=center, half_extents=Tuple(Float64.(obs["half_extents"])), marker=marker)
    end
end

# ── The 3-D mesh builder ─────────────────────────────────────────────────

# Collect the cells in the domain. Renumber the vertices that these cells
# use. Return the vertex array, the cell array, and the vertex labels.
function _compact3(M::NamedTuple)
    used = zeros(Int, M.nV[])
    nCells = 0
    for c in 1:M.nC[]
        (M.cellAlive[c] == 1 && M.cellIn[c] == 1) || continue
        nCells += 1
        for k in 1:4
            used[M.cellV[][c, k]] = 1
        end
    end
    (used[1] + used[2] + used[3] + used[4]) == 0 ||
        error("internal error: a cell in the domain uses a corner of the bounding tetrahedron")
    newIdx = zeros(Int, M.nV[])
    nPts = 0
    for v in 1:M.nV[]
        if used[v] == 1
            nPts += 1
            newIdx[v] = nPts
        end
    end
    pts = zeros(Float64, nPts, 3)
    labels = zeros(Int32, nPts)
    for v in 1:M.nV[]
        if used[v] == 1
            pts[newIdx[v], 1] = M.vx[v]
            pts[newIdx[v], 2] = M.vy[v]
            pts[newIdx[v], 3] = M.vz[v]
            labels[newIdx[v]] = M.vLabel[v]
        end
    end
    cells = zeros(Int, nCells, 4)
    k = 0
    for c in 1:M.nC[]
        (M.cellAlive[c] == 1 && M.cellIn[c] == 1) || continue
        k += 1
        for j in 1:4
            cells[k, j] = newIdx[M.cellV[][c, j]]
        end
    end
    return pts, cells, labels
end

# Return the largest radius-edge ratio of all cells.
function _max_ratio3(pts::Matrix{Float64}, cells::Matrix{Int})
    m = size(cells, 1)
    vx = pts[:, 1]; vy = pts[:, 2]; vz = pts[:, 3]
    ids = collect(1:m)
    ccx = zeros(Float64, m); ccy = zeros(Float64, m); ccz = zeros(Float64, m); ccR = zeros(Float64, m)
    hAt = ones(Float64, m); ratio = zeros(Float64, m); sizeRatio = zeros(Float64, m)
    cell_circumcenter_3d_kernel(cells, vx, vy, vz, ids, m, ccx, ccy, ccz, ccR)
    cell_quality_3d_kernel(cells, vx, vy, vz, ids, m, ccR, hAt, ratio, sizeRatio)
    return maximum(ratio)
end

# Return the smallest and the largest dihedral angle of all cells, in
# degrees. Also return the number of cells with a dihedral angle below
# goalRad. The goal is in radians.
function _dihedral_stats3(pts::Matrix{Float64}, cells::Matrix{Int}, goalRad::Float64, goalMaxRad::Float64 = Inf)
    m = size(cells, 1)
    minD = zeros(Float64, m)
    maxD = zeros(Float64, m)
    tet_dihedral_kernel(cells, pts, m, minD, maxD)
    nLow = count(<(goalRad), minD)
    nHigh = count(>(goalMaxRad), maxD)
    return rad2deg(minimum(minD)), rad2deg(maximum(maxD)), nLow, nHigh
end

# Compute four numbers for each tetrahedron: the smallest and the largest
# dihedral angle in radians, the signed volume, and the radius-edge ratio.
# The function calls the kernels that the writer uses everywhere.
function _cell_numbers3(pts::Matrix{Float64}, cells::Matrix{Int})
    n = size(cells, 1)
    minD = zeros(Float64, n)
    maxD = zeros(Float64, n)
    vol = zeros(Float64, n)
    tet_dihedral_kernel(cells, pts, n, minD, maxD)
    signed_measure_3d_kernel(cells, pts, n, vol)
    ids = collect(1:n)
    vx = pts[:, 1]; vy = pts[:, 2]; vz = pts[:, 3]
    ccx = zeros(Float64, n); ccy = zeros(Float64, n); ccz = zeros(Float64, n); ccR = zeros(Float64, n)
    cell_circumcenter_3d_kernel(cells, vx, vy, vz, ids, n, ccx, ccy, ccz, ccR)
    ratio = zeros(Float64, n)
    sizeRatio = zeros(Float64, n)
    cell_quality_3d_kernel(cells, vx, vy, vz, ids, n, ccR, ones(Float64, n), ratio, sizeRatio)
    return minD, maxD, vol, ratio
end

# Pack the three vertex numbers of a face into one integer. The integer
# does not depend on the order of the numbers. Each number must be below 2^21.
function _face_key3(a::Int, b::Int, c::Int)
    lo = min(a, min(b, c))
    hi = max(a, max(b, c))
    mid = a + b + c - lo - hi
    return (lo << 42) | (mid << 21) | hi
end

# Return the neighbor array of a cell list. Entry [c, f] is the cell on the
# other side of local face f of cell c. It is 0 for a boundary face.
function _cell_neighbors3(cells::Matrix{Int})
    n = size(cells, 1)
    keys = Vector{Int}(undef, 4 * n)
    for c in 1:n, f in 1:4
        fv = _face3(f)
        keys[4 * (c - 1) + f] = _face_key3(cells[c, fv[1]], cells[c, fv[2]], cells[c, fv[3]])
    end
    perm = sortperm(keys)
    nbr = zeros(Int, n, 4)
    j = 1
    while j < 4 * n
        if keys[perm[j]] == keys[perm[j + 1]]
            p = perm[j]; q = perm[j + 1]
            cp = div(p - 1, 4) + 1; fp = mod(p - 1, 4) + 1
            cq = div(q - 1, 4) + 1; fq = mod(q - 1, 4) + 1
            nbr[cp, fp] = cq
            nbr[cq, fq] = cp
            j += 2
        else
            j += 1
        end
    end
    return nbr
end

# Walk around the edge between local vertices i and j of cell c. Return
# the cells around the edge when they form a closed ring of exactly 3
# cells. Otherwise return an empty list. A boundary face or a longer ring
# gives an empty list.
function _edge_ring3(cells::Matrix{Int}, nbr::Matrix{Int}, c::Int, i::Int, j::Int)
    u = cells[c, i]; w = cells[c, j]
    keep = 0; leave = 0
    for t in 1:4
        (t == i || t == j) && continue
        keep == 0 ? (keep = cells[c, t]) : (leave = cells[c, t])
    end
    ring = Int[c]
    cur = c
    for step in 1:3
        # leave the cell through the face opposite the vertex `leave`
        fl = 0
        for t in 1:4
            cells[cur, t] == leave && (fl = t)
        end
        nxt = nbr[cur, fl]
        nxt == 0 && return Int[]
        if nxt == c
            return step == 3 ? ring : Int[]
        end
        step == 3 && return Int[]
        push!(ring, nxt)
        # the new cell holds the vertex `keep` and one new vertex z
        z = 0
        for t in 1:4
            v = cells[nxt, t]
            (v != u && v != w && v != keep) && (z = v)
        end
        leave = keep
        keep = z
        cur = nxt
    end
    return Int[]
end

# Return the score of each cell for the flip pass. Without a goal for the
# largest angle (goalMaxRad is Inf), the score is the smallest dihedral
# angle. With a goal, the score is the smaller of two margins. The first
# margin is the smallest angle minus its goal. The second margin is the
# goal of the largest angle minus the largest angle. A score below 0
# means that a goal is not met. When goalRad is 0, only the second margin
# counts.
function _angle_score3(minD::Vector{Float64}, maxD::Vector{Float64}, goalRad::Float64, goalMaxRad::Float64)
    isfinite(goalMaxRad) || return minD
    n = length(minD)
    score = zeros(Float64, n)
    for c in 1:n
        hi = goalMaxRad - maxD[c]
        score[c] = goalRad > 0 ? min(minD[c] - goalRad, hi) : hi
    end
    return score
end

# Make one pass of flips on a cell list. A flip changes the cells around
# one face (2 cells become 3) or around one edge (3 cells become 2). It
# does not move a vertex. The pass tries each flip that touches a cell
# with a dihedral angle below the goal.
# The pass keeps a flip only when all of these are true:
#   1. The new cells have a positive volume and fill the same space.
#   2. The smallest dihedral angle is at least gainMin higher.
#   3. The radius-edge ratio stays within the larger of the old ratio
#      and ratioBound.
# Flips that share a cell do not run in the same pass. Return the new
# cell list, the number of candidates and the number of flips done.
# The vector why counts the failed candidates by test. Entry 1 is for
# test 1. Entry 2 is for test 3. Entry 3 is for test 2.
function _flip_pass3(pts::Matrix{Float64}, cells::Matrix{Int}, goalRad::Float64, ratioBound::Float64, gainMin::Float64,
                     why::Vector{Int} = zeros(Int, 4), goalMaxRad::Float64 = Inf)
    nC = size(cells, 1)
    size(pts, 1) < 2^21 || return cells, 0, 0
    minD, maxD, vol, ratio = _cell_numbers3(pts, cells)
    bad = Int[c for c in 1:nC if minD[c] < goalRad || maxD[c] > goalMaxRad]
    isempty(bad) && return cells, 0, 0
    nbr = _cell_neighbors3(cells)

    # ---- List the candidates. Kind 1 is a 2-3 flip. Kind 3 is a 3-2 flip.
    newRows = Vector{NTuple{4,Int}}()
    newStart = Int[]; newCount = Int[]
    oldIdx = Int[]; oldStart = Int[]; oldCount = Int[]
    kind = Int[]
    addCand(k, tets, olds) = begin
        push!(kind, k)
        push!(newStart, length(newRows) + 1); push!(newCount, length(tets))
        append!(newRows, tets)
        push!(oldStart, length(oldIdx) + 1); push!(oldCount, length(olds))
        append!(oldIdx, olds)
    end
    for c in bad
        for f in 1:4
            d = nbr[c, f]
            d == 0 && continue
            fv = _face3(f)
            a = cells[c, f]
            p = cells[c, fv[1]]; q = cells[c, fv[2]]; r = cells[c, fv[3]]
            b = 0
            for t in 1:4
                v = cells[d, t]
                (v != p && v != q && v != r) && (b = v)
            end
            addCand(1, [(p, q, a, b), (q, r, a, b), (r, p, a, b)], [c, d])
        end
        for (i, j) in _tet_local_edges()
            ring = _edge_ring3(cells, nbr, c, i, j)
            isempty(ring) && continue
            u = cells[c, i]; w = cells[c, j]
            rv = Int[]
            for cc in ring, t in 1:4
                v = cells[cc, t]
                (v != u && v != w && !(v in rv)) && push!(rv, v)
            end
            length(rv) == 3 || continue
            addCand(3, [(rv[1], rv[2], rv[3], u), (rv[1], rv[2], rv[3], w)], ring)
        end
    end
    nG = length(kind)
    nG == 0 && return cells, 0, 0
    nNew = length(newRows)
    cand = Matrix{Int}(undef, nNew, 4)
    for k in 1:nNew
        cand[k, 1], cand[k, 2], cand[k, 3], cand[k, 4] = newRows[k]
    end

    # ---- Give each new cell a positive orientation. A candidate can fill
    # a different space than the old cells. Then the sum of the positive
    # volumes differs from the old sum, and the volume test below fails it.
    cvol = zeros(Float64, nNew)
    signed_measure_3d_kernel(cand, pts, nNew, cvol)
    for k in 1:nNew
        if cvol[k] < 0
            cand[k, 3], cand[k, 4] = cand[k, 4], cand[k, 3]
        end
    end
    cD, cDmax, cvol, cratio = _cell_numbers3(pts, cand)
    # The score of a cell is its smallest dihedral angle. With a goal for
    # the largest angle, the score is the smaller margin to the two goals.
    scoreOld = _angle_score3(minD, maxD, goalRad, goalMaxRad)
    scoreNew = _angle_score3(cD, cDmax, goalRad, goalMaxRad)

    # ---- Reduce the numbers of the new and of the old cells.
    idNew = collect(1:nNew)
    lo = zeros(Float64, nG); hi = zeros(Float64, nG); sm = zeros(Float64, nG)
    newVolMin = zeros(Float64, nG); newVolSum = zeros(Float64, nG)
    group_extrema_kernel(cvol, idNew, newStart, newCount, nG, newVolMin, hi, newVolSum)
    newDMin = zeros(Float64, nG)
    group_extrema_kernel(scoreNew, idNew, newStart, newCount, nG, newDMin, hi, sm)
    newRatioMax = zeros(Float64, nG)
    group_extrema_kernel(cratio, idNew, newStart, newCount, nG, lo, newRatioMax, sm)
    oldVolSum = zeros(Float64, nG)
    group_extrema_kernel(vol, oldIdx, oldStart, oldCount, nG, lo, hi, oldVolSum)
    oldDMin = zeros(Float64, nG)
    group_extrema_kernel(scoreOld, oldIdx, oldStart, oldCount, nG, oldDMin, hi, sm)
    oldRatioMax = zeros(Float64, nG)
    group_extrema_kernel(ratio, oldIdx, oldStart, oldCount, nG, lo, oldRatioMax, sm)

    # ---- Accept the good candidates. Sort them by gain.
    gain = fill(-Inf, nG)
    for g in 1:nG
        if !(newVolMin[g] > 1e-10 * oldVolSum[g] && abs(newVolSum[g] - oldVolSum[g]) <= 1e-9 * oldVolSum[g])
            why[1] += 1
            continue
        end
        if newRatioMax[g] > max(oldRatioMax[g], ratioBound)
            why[2] += 1
            continue
        end
        gn = newDMin[g] - oldDMin[g]
        if gn >= gainMin
            gain[g] = gn
        else
            why[3] += 1
        end
    end
    order = sortperm(gain; rev = true)

    # ---- Do the flips that do not share a cell.
    locked = falses(nC)
    dead = falses(nC)
    out = copy(cells)
    extra = Vector{NTuple{4,Int}}()
    nDone = 0
    for g in order
        gain[g] == -Inf && break
        olds = oldIdx[oldStart[g]:oldStart[g] + oldCount[g] - 1]
        any(locked[o] for o in olds) && continue
        for o in olds
            locked[o] = true
        end
        rows = newStart[g]:newStart[g] + newCount[g] - 1
        for (n, o) in enumerate(olds)
            if n <= length(rows)
                out[o, :] = cand[rows[n], :]
            else
                dead[o] = true
            end
        end
        if kind[g] == 1
            k = rows[3]
            push!(extra, (cand[k, 1], cand[k, 2], cand[k, 3], cand[k, 4]))
        end
        nDone += 1
    end
    keepRows = Int[c for c in 1:nC if !dead[c]]
    out = out[keepRows, :]
    if !isempty(extra)
        add = Matrix{Int}(undef, length(extra), 4)
        for (k, t) in enumerate(extra)
            add[k, 1], add[k, 2], add[k, 3], add[k, 4] = t
        end
        out = vcat(out, add)
    end
    return out, nG, nDone
end

# Build the tetrahedral mesh of a box with obstacles. Return the vertex
# array, the cell array, the edge array, the vertex markers, and a
# NamedTuple of statistics. The function needs no package.
#
# L: box lengths. per: periodic flags.
# obsKind, obsCx, ..., obsHz: flat obstacle arrays. Kind 1 is a sphere.
#   Kind 2 is a box.
# hMax, hObs, grading, ratioBound, sweeps, maxVertices: the mesh options
# minDihedral: goal for the smallest dihedral angle, in degrees. The flip
#   pass runs when it is above 0.
# maxDihedral: goal for the largest dihedral angle, in degrees. 0 is off.
#   The flip pass runs when it is above 0.
function _triangulate_domain3(logBuf::IOBuffer, L::NTuple{3,Float64}, per::NTuple{3,Bool},
                              obsKind::Vector{Int}, obsCx::Vector{Float64}, obsCy::Vector{Float64}, obsCz::Vector{Float64},
                              obsHx::Vector{Float64}, obsHy::Vector{Float64}, obsHz::Vector{Float64},
                              hMax::Float64, hObs::Float64, grading::Float64, ratioBound::Float64,
                              sweeps::Int, maxVertices::Int, minDihedral::Float64 = 0.0, maxDihedral::Float64 = 0.0)
    nObs = length(obsKind)
    hObsUsed = nObs == 0 ? hMax : hObs
    M = _new_mesher3(4096, 16384, 1024, 1024, maxVertices)
    M.obsKind[] = obsKind; M.obsCx[] = obsCx; M.obsCy[] = obsCy; M.obsCz[] = obsCz
    M.obsHx[] = obsHx; M.obsHy[] = obsHy; M.obsHz[] = obsHz
    M.nObs[] = nObs
    M.hMax[] = hMax; M.hObs[] = hObsUsed; M.grading[] = grading; M.ratioBound[] = ratioBound
    # The lattice quantum and the fine quantum are powers of two.
    e = ceil(Int, log2(max(L[1], L[2], L[3])))
    qLat = 2.0^(e - 38)
    M.qFine[] = 2.0^(e - 46)
    # ---- Build the input geometry. The lattice of each box has an exact
    # cell width. The width is a multiple of qLat. The box size changes by
    # less than 1e-11 of its length.
    nBox = (max(1, ceil(Int, L[1] / hMax)), max(1, ceil(Int, L[2] / hMax)), max(1, ceil(Int, L[3] / hMax)))
    dBox = (round(L[1] / nBox[1] / qLat) * qLat, round(L[2] / nBox[2] / qLat) * qLat, round(L[3] / nBox[3] / qLat) * qLat)
    Lq = (nBox[1] * dBox[1], nBox[2] * dBox[2], nBox[3] * dBox[3])
    c1 = _start_bounding_tetra3!(M, 0.0, 0.0, 0.0, Lq[1], Lq[2], Lq[3])
    _build_box_plc3!(M, (0.0, 0.0, 0.0), dBox, nBox, (1, 2, 3, 4, 5, 6))
    for i in 1:nObs
        label = 100 + (i - 1)
        if obsKind[i] == 1
            level = _sphere_level3(obsHx[i], hObsUsed)
            _build_sphere_plc3!(M, obsCx[i], obsCy[i], obsCz[i], obsHx[i], level, label)
        else
            nn = (max(1, ceil(Int, 2.0 * obsHx[i] / hObsUsed)), max(1, ceil(Int, 2.0 * obsHy[i] / hObsUsed)),
                  max(1, ceil(Int, 2.0 * obsHz[i] / hObsUsed)))
            lo = (round((obsCx[i] - obsHx[i]) / qLat) * qLat, round((obsCy[i] - obsHy[i]) / qLat) * qLat,
                  round((obsCz[i] - obsHz[i]) / qLat) * qLat)
            dd = (round(2.0 * obsHx[i] / nn[1] / qLat) * qLat, round(2.0 * obsHy[i] / nn[2] / qLat) * qLat,
                  round(2.0 * obsHz[i] / nn[3] / qLat) * qLat)
            _build_box_plc3!(M, lo, dd, nn, (label, label, label, label, label, label))
        end
    end
    _setup_periodic3!(M, (0.0, 0.0, 0.0), Lq, per)

    # ---- Insert the input vertices.
    start = c1
    for v in 5:M.nV[]
        _insert3!(M, v, start, 0, 0, 0, 0) == 0 || error("internal error: plain insertion failed")
        start = M.vCell[v]
    end

    # ---- Make every segment an edge and every subfacet a face.
    _recover3!(M)

    # ---- Set the domain flags. Then refine the cells in the domain.
    _classify3!(M)
    M.doBad[] = true
    ids = Int[c for c in 1:M.nC[] if M.cellAlive[c] == 1 && M.cellIn[c] == 1]
    isempty(ids) && error("the domain has no cell. Check the box and the obstacles.")
    _scan3!(M, ids, length(ids))
    _run_refinement3!(M)
    # A subfacet that an insertion deletes returns later. Between the two
    # events, the flags of the cells next to it can be wrong. Repeat the
    # steps below until a pass inserts no vertex. First, make every
    # segment an edge and every subfacet a face. Then set the flags again
    # from the finished surface. Then refine the cells that the new flags
    # put in the domain.
    rounds = 0
    while true
        rounds += 1
        rounds > 50 && error("internal error: the domain flags do not become stable")
        inner = 0
        while _mark_present3!(M) > 0
            inner += 1
            inner > 100 && error("internal error: a segment or subfacet is missing after refinement")
            _run_refinement3!(M)
        end
        _classify3!(M)
        nBefore = M.nV[]
        ids = Int[c for c in 1:M.nC[] if M.cellAlive[c] == 1 && M.cellIn[c] == 1]
        _scan3!(M, ids, length(ids))
        _run_refinement3!(M)
        M.nV[] == nBefore && break
    end
    exactShare = M.insCount[1] == 0 ? 0.0 : M.insCount[2] / M.insCount[1]

    # ---- Remove the cells outside the domain.
    pts, cells, labels = _compact3(M)
    nPts = size(pts, 1)
    nCells = size(cells, 1)
    meas = zeros(Float64, nCells)
    signed_measure_3d_kernel(cells, pts, nCells, meas)
    nNeg = zeros(Int, 1)
    count_negative_kernel(meas, nCells, nNeg)
    if nNeg[1] > 0
        _log(logBuf, "  [warn] tetra: $(nNeg[1]) cells had a negative volume. The writer reordered them.")
        fix_orientation_kernel(meas, nCells, 4, cells)
    end

    # ---- Flip the cells around slivers. Each pass flips faces and edges
    # of the cells below the goal. The vertices stay where they are.
    reportRad = deg2rad(minDihedral > 0 ? minDihedral : 10.0)
    goalMaxRad = maxDihedral > 0 ? deg2rad(maxDihedral) : Inf
    nBelowStart = _dihedral_stats3(pts, cells, reportRad)[3]
    nAboveStart = _dihedral_stats3(pts, cells, reportRad, goalMaxRad)[4]
    flipsDone = 0
    if minDihedral > 0 || maxDihedral > 0
        for pass in 1:10
            why = zeros(Int, 4)
            cells, nTried, nDone = _flip_pass3(pts, cells, deg2rad(minDihedral), ratioBound, deg2rad(1.0), why, goalMaxRad)
            _log(logBuf, "  flip pass $pass: $nTried candidates, $nDone flips. Rejected: $(why[1]) by volume, $(why[2]) by radius-edge ratio, $(why[3]) by angle gain")
            flipsDone += nDone
            nDone == 0 && break
        end
        nCells = size(cells, 1)
        meas = zeros(Float64, nCells)
    end
    edges, c2e = _build_edge_registry(cells, 3)
    nEdges = size(edges, 1)

    # ---- Smooth the interior vertices. A sweep stays only when no cell
    # flips and the smallest dihedral angle does not shrink.
    deg = zeros(Int, nPts)
    for e in 1:nEdges
        deg[edges[e, 1]] += 1
        deg[edges[e, 2]] += 1
    end
    nbrOff = ones(Int, nPts + 1)
    for v in 1:nPts
        nbrOff[v + 1] = nbrOff[v] + deg[v]
    end
    nbrIdx = zeros(Int, nbrOff[nPts + 1] - 1)
    cursor = copy(nbrOff[1:nPts])
    for e in 1:nEdges
        u = edges[e, 1]; w = edges[e, 2]
        nbrIdx[cursor[u]] = w; cursor[u] += 1
        nbrIdx[cursor[w]] = u; cursor[w] += 1
    end
    avgPts = zeros(Float64, nPts, 3)
    newPts = zeros(Float64, nPts, 3)
    v2cOff, v2cIdx = _vertex_to_cell(nPts, cells)
    minDeg, maxDeg, nLow, nHigh = _dihedral_stats3(pts, cells, reportRad, goalMaxRad)
    sweepsDone = 0
    for _ in 1:sweeps
        smooth_jacobi_3d_kernel(pts, nbrOff, nbrIdx, nPts, avgPts)
        # Try a full half step. Halve the step when the sweep makes the
        # smallest dihedral angle worse or flips a cell.
        accepted = false
        omega = 0.5
        for _try in 1:3
            smooth_accept_3d_kernel(pts, avgPts, labels, v2cOff, v2cIdx, cells, nPts, omega, newPts)
            signed_measure_3d_kernel(cells, newPts, nCells, meas)
            nNeg[1] = 0
            count_negative_kernel(meas, nCells, nNeg)
            if nNeg[1] == 0
                minNew, maxNew, nLowNew, nHighNew = _dihedral_stats3(newPts, cells, reportRad, goalMaxRad)
                # With a goal for the largest angle, the largest angle may not grow.
                if minNew >= minDeg && (maxDihedral == 0 || maxNew <= maxDeg)
                    pts, newPts = newPts, pts
                    minDeg, maxDeg, nLow, nHigh = minNew, maxNew, nLowNew, nHighNew
                    accepted = true
                    break
                end
            end
            omega *= 0.5
        end
        accepted || break
        sweepsDone += 1
    end
    maxRatio = _max_ratio3(pts, cells)

    info = (nVertices = nPts, nCells = nCells, nEdges = nEdges, minDihedralDeg = minDeg,
            maxDihedralDeg = maxDeg, nBelowGoal = nLow, nBelowStart = nBelowStart,
            nAboveGoal = nHigh, nAboveStart = nAboveStart,
            goalDeg = rad2deg(reportRad), flipsDone = flipsDone, maxRatio = maxRatio,
            sweepsDone = sweepsDone, exactShare = exactShare, boxSize = Lq)
    return pts, cells, edges, c2e, labels, info
end

# =====================================================================
# builders, config check, periodicity -- driver code, no package needed
# =====================================================================

# ── Helpers of the structured builders ───────────────────────────────────

# Turn a diagonal name into the integer mode code the kernels take.
function _diagonal_mode_2d(name::AbstractString)
    name == "backslash" && return 0
    name == "slash" && return 1
    name == "alternate" && return 2
    error("diagonal must be one of (\"backslash\", \"slash\", \"alternate\") for a 2-D mesh, got '$name'")
end

function _diagonal_mode_3d(name::AbstractString)
    (name == "main" || name == "fixed") && return 0
    name == "flip_x" && return 1
    name == "flip_y" && return 2
    name == "flip_z" && return 3
    name == "alternate_layers" && return 4
    error("diagonal must be one of (\"main\", \"flip_x\", \"flip_y\", \"flip_z\", " *
          "\"alternate_layers\") for a 3-D mesh, got '$name'")
end

# 2-D only: gridsplit and equilateral both need edge markers.
# Structured grid indices, recovered from a point's flat 0-based
# index, tell directly whether an edge lies on one box side.
function _edge_markers_grid_2d(edges::Matrix{Int}, nx::Int, ny::Int)
    nCols = nx + 1
    edgeMk = zeros(Int32, size(edges, 1))
    for r in 1:size(edges, 1)
        # edges holds 1-based Julia row indices; divrem needs the
        # underlying 0-based flat grid index
        u, v = edges[r, 1] - 1, edges[r, 2] - 1
        ju, iu = divrem(u, nCols)
        jv, iv = divrem(v, nCols)
        marker = 0
        if ju == 0 && jv == 0
            marker = 1
        elseif iu == nx && iv == nx
            marker = 2
        elseif ju == ny && jv == ny
            marker = 3
        elseif iu == 0 && iv == 0
            marker = 4
        end
        edgeMk[r] = marker
    end
    return edgeMk
end

# Fill nx+1 breakpoints, scaled to [0, L], from cfg's stretching keys
# for axis d (1-based). Falls back to uniform spacing when cfg carries
# no stretching keys at all (the equilateral case).
function _axis_positions(logBuf::IOBuffer, cfg::Dict, d::Int, n::Int, L::Float64)
    pos = zeros(Float64, n + 1)
    ratio = haskey(cfg, "stretching_ratio") ? cfg["stretching_ratio"][d] : 1.0
    symmetricBool = haskey(cfg, "stretching_symmetric") ? cfg["stretching_symmetric"][d] : false
    symmetric = symmetricBool ? 1 : 0
    if symmetricBool && mod(n, 2) != 0
        _log(logBuf, "  [warn] stretching_symmetric on axis $(_dirs()[d]): n_divisions=$n is odd -- " *
             "the two mirrored halves will have slightly different cell counts, giving a " *
             "small kink at the center; use an even value for a perfectly symmetric grid")
    end
    stretch_positions_kernel(n, ratio, symmetric, pos)
    for i in 1:n + 1
        pos[i] *= L
    end
    return pos
end

function _check_cell_orientation(logBuf::IOBuffer, signedMeasure::Vector{Float64}, meshType::AbstractString)
    nPos = count(>(0), signedMeasure)
    nNeg = count(<(0), signedMeasure)
    if nPos > 0 && nNeg > 0
        _log(logBuf, "  [warn] $meshType: cells are not consistently oriented ($nPos positive, $nNeg negative signed measure)")
    end
end

# ── Builders ─────────────────────────────────────────────────────────────
#
# A builder takes the checked config. It returns one NamedTuple, made by
# _builder_result. The fields are:
#   pts             vertex array with dim columns
#   cells           cell array, 1-based, positive orientation
#   ptMk            Int32 vertex markers
#   edgeMk          Int32 edge markers (2-D), or nothing
#   edges, c2e      edge array and cell-to-edge array, or nothing. A
#                   builder that has them already (the Delaunay types
#                   need them for smoothing) passes them on. When they
#                   are nothing, main builds them once.
#   obstacles       vector of obstacle records, can be empty
#   boxSize         exact box lengths that the vertices use, for the
#                   periodic search
#   info            NamedTuple of statistics, or nothing
#   edgeMarkerRule  :given (edgeMk is ready), :grid (main computes the
#                   edge markers of a 2-D grid mesh once the edge array
#                   exists), or :none (3-D, no edge markers)
#   edgeMarkerArgs  extra arguments of the rule, a Tuple
function _builder_result(pts, cells, ptMk, edgeMk, edges, c2e, obstacles, boxSize, info,
                         edgeMarkerRule::Symbol, edgeMarkerArgs::Tuple)
    return (pts = pts, cells = cells, ptMk = ptMk, edgeMk = edgeMk, edges = edges, c2e = c2e,
            obstacles = obstacles, boxSize = boxSize, info = info,
            edgeMarkerRule = edgeMarkerRule, edgeMarkerArgs = edgeMarkerArgs)
end

function _build_gridsplit(logBuf::IOBuffer, cfg::Dict)
    Lx, Ly = cfg["box_length"]
    nx, ny = cfg["n_divisions"]
    mode = _diagonal_mode_2d(cfg["diagonal"])
    px = _axis_positions(logBuf, cfg, 1, nx, Lx)
    py = _axis_positions(logBuf, cfg, 2, ny, Ly)

    nPts = (nx + 1) * (ny + 1)
    pts = zeros(Float64, nPts, 2)
    grid_points_2d_kernel(nx, ny, px, py, pts)

    nCells = 2 * nx * ny
    cells0 = zeros(Int, nCells, 3)
    grid_triangles_2d_kernel(nx, ny, mode, cells0)

    meas = zeros(Float64, nCells)
    signed_measure_2d_kernel(cells0, pts, nCells, meas)
    fix_orientation_kernel(meas, nCells, 3, cells0)
    # Measure again, so the check sees the cells after the fix.
    signed_measure_2d_kernel(cells0, pts, nCells, meas)
    _check_cell_orientation(logBuf, meas, "gridsplit")

    ptMk = zeros(Int32, nPts)
    point_markers_grid_2d_kernel(nx, ny, ptMk)

    return _builder_result(pts, cells0, ptMk, nothing, nothing, nothing, NamedTuple[],
                           cfg["box_length"], nothing, :grid, (nx, ny))
end

function _build_equilateral(logBuf::IOBuffer, cfg::Dict)
    Lx, Ly = cfg["box_length"]
    nx, ny = cfg["n_divisions"]
    dx = Lx / nx
    dy = Ly / ny
    dyEq = dx * sqrt(3.0) / 2.0
    if abs(dy - dyEq) / dyEq > 0.02
        _log(logBuf, "  [info] equilateral: the interior triangles are equilateral at dy = " *
             "$(round(dyEq, digits = 5)) (Ly/ny = $(round(dy, digits = 5)))")
    end
    if mod(ny, 2) != 0 && cfg["periodic"][2]
        _log(logBuf, "  [warn] equilateral: periodic_y requested with an odd n_divisions[2]; " *
             "the top and bottom rows will not line up exactly")
    end

    nPts = (div(ny, 2) + 1) * (nx + 1) + div(ny + 1, 2) * (nx + 2)
    pts = zeros(Float64, nPts, 2)
    equilateral_points_kernel(nx, ny, Lx, Ly, pts)

    nCells = ny * (2 * nx + 1)
    cells0 = zeros(Int, nCells, 3)
    equilateral_triangles_kernel(nx, ny, cells0)

    meas = zeros(Float64, nCells)
    signed_measure_2d_kernel(cells0, pts, nCells, meas)
    fix_orientation_kernel(meas, nCells, 3, cells0)
    # Measure again, so the check sees the cells after the fix.
    signed_measure_2d_kernel(cells0, pts, nCells, meas)
    _check_cell_orientation(logBuf, meas, "equilateral")

    # The rows have different lengths, so the markers come from the
    # coordinates, as for the Delaunay mesh.
    edges, c2e = _build_edge_registry(cells0, 2)
    nEdges = size(edges, 1)
    edgeCount = zeros(Int, nEdges)
    for ci in 1:nCells, li in 1:3
        edgeCount[c2e[ci, li] + 1] += 1
    end
    isBnd = zeros(Int, nEdges)
    for e in 1:nEdges
        edgeCount[e] == 1 && (isBnd[e] = 1)
    end
    tol = 1.0e-9 * max(Lx, Ly)
    noKind = zeros(Int, 0)
    noVal = zeros(Float64, 0)
    ptMk = zeros(Int32, nPts)
    point_markers_2d_domain_kernel(pts, nPts, Lx, Ly, tol, noKind, noVal, noVal, noVal, noVal, 0, ptMk)
    edgeMk = zeros(Int32, nEdges)
    edge_markers_2d_domain_kernel(pts, edges, isBnd, nEdges, Lx, Ly, tol, noKind, noVal, noVal, noVal, noVal, 0, edgeMk)

    return _builder_result(pts, cells0, ptMk, edgeMk, edges, c2e, NamedTuple[],
                           cfg["box_length"], nothing, :given, ())
end

function _build_kuhn(logBuf::IOBuffer, cfg::Dict)
    Lx, Ly, Lz = cfg["box_length"]
    nx, ny, nz = cfg["n_divisions"]
    mode = _diagonal_mode_3d(cfg["diagonal"])
    px = _axis_positions(logBuf, cfg, 1, nx, Lx)
    py = _axis_positions(logBuf, cfg, 2, ny, Ly)
    pz = _axis_positions(logBuf, cfg, 3, nz, Lz)

    nPts = (nx + 1) * (ny + 1) * (nz + 1)
    pts = zeros(Float64, nPts, 3)
    grid_points_3d_kernel(nx, ny, nz, px, py, pz, pts)

    nCells = 6 * nx * ny * nz
    cells0 = zeros(Int, nCells, 4)
    kuhn_tets_kernel(nx, ny, nz, mode, _kuhn_perm_axis(), _kuhn_perm_parity(), cells0)

    meas = zeros(Float64, nCells)
    signed_measure_3d_kernel(cells0, pts, nCells, meas)
    fix_orientation_kernel(meas, nCells, 4, cells0)
    # Measure again, so the check sees the cells after the fix.
    signed_measure_3d_kernel(cells0, pts, nCells, meas)
    _check_cell_orientation(logBuf, meas, "kuhn")

    ptMk = zeros(Int32, nPts)
    point_markers_grid_3d_kernel(nx, ny, nz, ptMk)

    return _builder_result(pts, cells0, ptMk, nothing, nothing, nothing, NamedTuple[],
                           cfg["box_length"], nothing, :none, ())
end

# Build one triangle mesh of a box with obstacles from a config.
function _build_triangle(logBuf::IOBuffer, cfg::Dict)
    Lx, Ly = cfg["box_length"]
    per = cfg["periodic"]
    obstacles = cfg["obstacles"]
    opts = cfg["mesh_options"]
    hObs = isempty(obstacles) ? opts.meshSize : opts.meshSizeObstacle
    _validate_obstacles2(logBuf, obstacles, cfg["box_length"], hObs)

    nObs = length(obstacles)
    obsKind = zeros(Int, nObs)
    obsCx = zeros(Float64, nObs)
    obsCy = zeros(Float64, nObs)
    obsHx = zeros(Float64, nObs)
    obsHy = zeros(Float64, nObs)
    obstacleMeta = NamedTuple[]
    for (i, obs) in enumerate(obstacles)
        c = Float64.(obs["center"])
        obsCx[i] = c[1]
        obsCy[i] = c[2]
        if obs["shape"] == "circle"
            obsKind[i] = 1
            obsHx[i] = Float64(obs["radius"])
            obsHy[i] = Float64(obs["radius"])
        else
            obsKind[i] = 2
            he = Float64.(obs["half_extents"])
            obsHx[i] = he[1]
            obsHy[i] = he[2]
        end
        push!(obstacleMeta, _obstacle_record2(obs, 10 + (i - 1)))
    end

    pts, cells, edges, c2e, ptMk, edgeMk, info = _triangulate_domain2(logBuf,
        Lx, Ly, per[1], per[2], obsKind, obsCx, obsCy, obsHx, obsHy,
        opts.meshSize, hObs, opts.grading, opts.minAngle,
        opts.smoothingSweeps, opts.maxVertices)
    _log(logBuf, "  triangle: $(info.nVertices) vertices, $(info.nCells) cells")
    _log(logBuf, "  smallest cell angle: $(round(info.minAngleDeg, digits = 2)) degrees")
    _log(logBuf, "  smoothing sweeps kept: $(info.sweepsDone) of $(opts.smoothingSweeps)")
    return _builder_result(pts, cells, ptMk, edgeMk, edges, c2e, obstacleMeta,
                           cfg["box_length"], info, :given, ())
end

# Build one tetrahedral mesh of a box with obstacles from a config.
function _build_tetra(logBuf::IOBuffer, cfg::Dict)
    L = cfg["box_length"]
    per = cfg["periodic"]
    obstacles = cfg["obstacles"]
    opts = cfg["mesh_options"]
    hObs = isempty(obstacles) ? opts.meshSize : opts.meshSizeObstacle
    _validate_obstacles3(logBuf, obstacles, L, hObs)

    nObs = length(obstacles)
    obsKind = zeros(Int, nObs)
    obsCx = zeros(Float64, nObs); obsCy = zeros(Float64, nObs); obsCz = zeros(Float64, nObs)
    obsHx = zeros(Float64, nObs); obsHy = zeros(Float64, nObs); obsHz = zeros(Float64, nObs)
    obstacleMeta = NamedTuple[]
    for (i, obs) in enumerate(obstacles)
        c = Float64.(obs["center"])
        obsCx[i] = c[1]; obsCy[i] = c[2]; obsCz[i] = c[3]
        if obs["shape"] == "sphere"
            obsKind[i] = 1
            obsHx[i] = Float64(obs["radius"]); obsHy[i] = obsHx[i]; obsHz[i] = obsHx[i]
        else
            obsKind[i] = 2
            he = Float64.(obs["half_extents"])
            obsHx[i] = he[1]; obsHy[i] = he[2]; obsHz[i] = he[3]
        end
        push!(obstacleMeta, _obstacle_record3(obs, 100 + (i - 1)))
    end

    pts, cells, edges, c2e, ptMk, info = _triangulate_domain3(logBuf,
        L, per, obsKind, obsCx, obsCy, obsCz, obsHx, obsHy, obsHz,
        opts.meshSize, hObs, opts.grading, opts.radiusEdgeBound,
        opts.smoothingSweeps, opts.maxVertices, opts.minDihedral, opts.maxDihedral)
    _log(logBuf, "  tetra: $(info.nVertices) vertices, $(info.nCells) cells")
    _log(logBuf, "  dihedral angles: smallest $(round(info.minDihedralDeg, digits = 2)), largest $(round(info.maxDihedralDeg, digits = 2)) degrees")
    _log(logBuf, "  cells with a dihedral angle below $(info.goalDeg) degrees: $(info.nBelowStart) before the flips, $(info.nBelowGoal) at the end ($(info.flipsDone) flips)")
    if opts.minDihedral > 0 && info.nBelowGoal > 0
        _log(logBuf, "  [warn] tetra: $(info.nBelowGoal) cells have a dihedral angle below the goal of $(opts.minDihedral) degrees. The flips could not remove them.")
    end
    if opts.maxDihedral > 0
        _log(logBuf, "  cells with a dihedral angle above $(opts.maxDihedral) degrees: $(info.nAboveStart) before the flips, $(info.nAboveGoal) at the end")
        if info.nAboveGoal > 0
            _log(logBuf, "  [warn] tetra: $(info.nAboveGoal) cells have a dihedral angle above the goal of $(opts.maxDihedral) degrees. The flips could not remove them.")
        end
    end
    _log(logBuf, "  largest radius-edge ratio after smoothing: $(round(info.maxRatio, digits = 3))")
    _log(logBuf, "  smoothing sweeps kept: $(info.sweepsDone) of $(opts.smoothingSweeps)")
    _log(logBuf, "  exact insphere calls: $(round(100 * info.exactShare, digits = 1)) percent")
    return _builder_result(pts, cells, ptMk, nothing, edges, c2e, obstacleMeta,
                           info.boxSize, info, :none, ())
end

# Call the builder of the mesh type in a checked config. The call goes
# through invokelatest on purpose. Without it, Julia infers all five
# builders, and so both Delaunay meshers, the first time that main runs,
# for every mesh type. With it, Julia compiles only the builder that
# runs. The result has the type Any here. main asserts the types.
function _build_mesh(logBuf::IOBuffer, cfg::Dict)
    mt = cfg["mesh_type"]
    mt == "gridsplit" && return Base.invokelatest(_build_gridsplit, logBuf, cfg)
    mt == "equilateral" && return Base.invokelatest(_build_equilateral, logBuf, cfg)
    mt == "kuhn" && return Base.invokelatest(_build_kuhn, logBuf, cfg)
    mt == "triangle" && return Base.invokelatest(_build_triangle, logBuf, cfg)
    mt == "tetra" && return Base.invokelatest(_build_tetra, logBuf, cfg)
    error("mesh_type must be one of $(join(_mesh_types(), ", ")), got '$mt'")
end

# Compute the edge markers of a 2-D grid mesh from the edge array.
function _grid_edge_markers(res, edges::Matrix{Int})
    res.edgeMk !== nothing && return res.edgeMk
    if res.edgeMarkerRule == :grid
        nx, ny = res.edgeMarkerArgs
        return _edge_markers_grid_2d(edges, nx, ny)
    end
    return nothing
end

# ── Hierarchy: nested refinement ────────────────────────────────────────
#
# Level 0 is the coarsest mesh. Level k comes from k red refinements.
# Every cell splits into 4 (2-D) or 8 (3-D) children. Every edge gets one
# new vertex at its midpoint. The fine vertices 1 to nV are the coarse
# vertices. The midpoint of coarse edge e is fine vertex nV + e. Let c be
# 4 in 2-D and 8 in 3-D. The children of coarse cell k are the fine rows
# c (k - 1) + 1 to c k. No flip and no smoothing run on a fine level,
# because either one would break the nesting.

# Return the parents of the fine vertices of one refinement. A coarse
# vertex has its own number twice. A midpoint has the two end vertices
# of its coarse edge. All numbers are 1-based.
function _vertex_parents(edges::Matrix{Int}, nV::Int)
    nE = size(edges, 1)
    vp = zeros(Int, nV + nE, 2)
    for v in 1:nV
        vp[v, 1] = v
        vp[v, 2] = v
    end
    for e in 1:nE
        vp[nV + e, 1] = edges[e, 1]
        vp[nV + e, 2] = edges[e, 2]
    end
    return vp
end

# Return the parent cell of each child. The cell numbers are 1-based.
function _cell_parents(nCells::Int, children::Int)
    cp = zeros(Int, children * nCells)
    for k in 1:nCells, j in 1:children
        cp[children * (k - 1) + j] = k
    end
    return cp
end

# Refine a 2-D result record one time. The record must hold edges, c2e
# and edgeMk. Return the fine record, with its own edges, c2e and edgeMk,
# and a NamedTuple with the two parent arrays. A midpoint takes the
# marker of its coarse edge. A fine edge that is half of a coarse edge
# takes the marker of that edge. Every other fine edge has marker 0.
function _refine_mesh2(res)
    pts = res.pts::Matrix{Float64}
    cells = res.cells::Matrix{Int}
    ptMk = res.ptMk::Vector{Int32}
    edges = res.edges::Matrix{Int}
    c2e = res.c2e::Matrix{Int}
    edgeMk = res.edgeMk::Vector{Int32}
    nV = size(pts, 1)
    nC = size(cells, 1)
    nE = size(edges, 1)
    nF = nV + nE

    ptsF = zeros(Float64, nF, 2)
    ptsF[1:nV, :] = pts
    edge_midpoint_kernel(pts, edges, nV, nE, 2, ptsF)

    cellsF = zeros(Int, 4 * nC, 3)
    refine_triangles_2d_kernel(cells, c2e, nV, nC, cellsF)
    # The children of a positive cell are positive. Keep the check, so
    # that a cell with a wrong orientation cannot pass.
    meas = zeros(Float64, 4 * nC)
    signed_measure_2d_kernel(cellsF, hcat(ptsF, zeros(Float64, nF)), 4 * nC, meas)
    fix_orientation_kernel(meas, 4 * nC, 3, cellsF)

    ptMkF = zeros(Int32, nF)
    ptMkF[1:nV] = ptMk
    for e in 1:nE
        ptMkF[nV + e] = edgeMk[e]
    end

    edgesF, c2eF = _build_edge_registry(cellsF, 2)
    edgeMkF = zeros(Int32, size(edgesF, 1))
    for f in 1:size(edgesF, 1)
        u = edgesF[f, 1]
        w = edgesF[f, 2]
        # A half edge joins one midpoint to one end of its coarse edge.
        if u > nV && w <= nV
            mid, other = u, w
        elseif w > nV && u <= nV
            mid, other = w, u
        else
            continue
        end
        ec = mid - nV
        if other == edges[ec, 1] || other == edges[ec, 2]
            edgeMkF[f] = edgeMk[ec]
        end
    end

    resF = merge(res, (pts = ptsF, cells = cellsF, ptMk = ptMkF, edgeMk = edgeMkF,
                       edges = edgesF, c2e = c2eF, info = nothing))
    parents = (vertexParents = _vertex_parents(edges, nV), cellParent = _cell_parents(nC, 4))
    return resF, parents
end

# Return the table of the 8 children of a tetrahedron for each of the 3
# diagonals of the octahedron. The entry (dg, j, i) is the code of vertex
# i of child j. The codes 1 to 4 are the vertices of the parent. The codes
# 5 to 10 are the midpoints of the local edges (1,2), (1,3), (1,4),
# (2,3), (2,4) and (3,4). The children 1 to 4 are the corner children.
# The children 5 to 8 cut the octahedron along its diagonal. The three
# diagonals join the midpoints 5-10, 6-9 and 7-8.
function _bey_table()
    corners = [1 5 6 7; 2 5 8 9; 3 6 8 10; 4 7 9 10]
    inner = [5 10 6 7; 5 10 7 9; 5 10 9 8; 5 10 8 6;
             6 9 5 7; 6 9 7 10; 6 9 10 8; 6 9 8 5;
             7 8 5 6; 7 8 6 10; 7 8 10 9; 7 8 9 5]
    bey = zeros(Int, 3, 8, 4)
    for dg in 1:3, j in 1:4, i in 1:4
        bey[dg, j, i] = corners[j, i]
        bey[dg, 4 + j, i] = inner[4 * (dg - 1) + j, i]
    end
    return bey
end

# Return the marker of a boundary face with the vertices a, b and c. A
# face with three vertices of one obstacle (marker 100 or more) has the
# marker of the obstacle. A face on a box side has the number of the
# side, and the lowest number wins. The function returns 0 when no rule
# applies.
function _boundary_face_marker3(pts::Matrix{Float64}, ptMk::Vector{Int32}, box, a::Int, b::Int, c::Int)
    ma = Int(ptMk[a])
    if ma >= 100 && ma == ptMk[b] && ma == ptMk[c]
        return ma
    end
    tol = 1e-9 * max(box[1], box[2], box[3])
    for d in 1:3
        if abs(pts[a, d]) < tol && abs(pts[b, d]) < tol && abs(pts[c, d]) < tol
            return 2 * d - 1
        end
        if abs(pts[a, d] - box[d]) < tol && abs(pts[b, d] - box[d]) < tol && abs(pts[c, d] - box[d]) < tol
            return 2 * d
        end
    end
    return 0
end

# Refine a 3-D result record one time. The record must hold edges and c2e.
# Return the fine record and a NamedTuple with the two parent arrays. The
# argument fixedDiagonal is 0 for the max-min rule. The values 1, 2 and 3
# force one diagonal for every cell. The marker of a midpoint is the lowest
# marker of the boundary faces that contain its edge. A boundary face has
# one cell.
function _refine_mesh3(res, fixedDiagonal::Int = 0)
    pts = res.pts::Matrix{Float64}
    cells = res.cells::Matrix{Int}
    ptMk = res.ptMk::Vector{Int32}
    edges = res.edges::Matrix{Int}
    c2e = res.c2e::Matrix{Int}
    nV = size(pts, 1)
    nC = size(cells, 1)
    nE = size(edges, 1)
    nF = nV + nE
    0 <= fixedDiagonal <= 3 || error("fixedDiagonal must be 0, 1, 2 or 3")

    ptsF = zeros(Float64, nF, 3)
    ptsF[1:nV, :] = pts
    edge_midpoint_kernel(pts, edges, nV, nE, 3, ptsF)

    # Choose the diagonal of each cell. The midpoints are known, so the
    # shape of every candidate child can be measured.
    bey = _bey_table()
    dgSel = fill(max(fixedDiagonal, 1), nC)
    if fixedDiagonal == 0
        cand = zeros(Int, 8 * nC, 4)
        dgTry = zeros(Int, nC)
        eta = [zeros(Float64, 8 * nC) for _ in 1:3]
        for dg in 1:3
            fill!(dgTry, dg)
            refine_tets_3d_kernel(cells, c2e, bey, dgTry, nV, nC, cand)
            tet_shape_kernel(cand, ptsF, 8 * nC, eta[dg])
        end
        diagonal_choice_kernel(eta[1], eta[2], eta[3], nC, dgSel)
    end

    cellsF = zeros(Int, 8 * nC, 4)
    refine_tets_3d_kernel(cells, c2e, bey, dgSel, nV, nC, cellsF)
    meas = zeros(Float64, 8 * nC)
    signed_measure_3d_kernel(cellsF, ptsF, 8 * nC, meas)
    fix_orientation_kernel(meas, 8 * nC, 4, cellsF)

    # Markers of the midpoints.
    nbr = _cell_neighbors3(cells)
    localEdge = zeros(Int, 4, 4)
    for (li, (p, q)) in enumerate(_tet_local_edges())
        localEdge[p, q] = li
        localEdge[q, p] = li
    end
    ptMkF = zeros(Int32, nF)
    ptMkF[1:nV] = ptMk
    nNoMarker = 0
    for c in 1:nC, f in 1:4
        nbr[c, f] == 0 || continue
        fv = _face3(f)
        mk = _boundary_face_marker3(pts, ptMk, res.boxSize, cells[c, fv[1]], cells[c, fv[2]], cells[c, fv[3]])
        if mk == 0
            nNoMarker += 1
            continue
        end
        for (p, q) in ((fv[1], fv[2]), (fv[1], fv[3]), (fv[2], fv[3]))
            m = nV + c2e[c, localEdge[p, q]] + 1
            cur = ptMkF[m]
            ptMkF[m] = cur == 0 ? Int32(mk) : Int32(min(Int(cur), mk))
        end
    end
    nNoMarker == 0 || error("hierarchy: $nNoMarker boundary faces have no marker. The marker rule found neither an obstacle nor a box side.")

    edgesF, c2eF = _build_edge_registry(cellsF, 3)
    resF = merge(res, (pts = ptsF, cells = cellsF, ptMk = ptMkF, edges = edgesF, c2e = c2eF, info = nothing))
    parents = (vertexParents = _vertex_parents(edges, nV), cellParent = _cell_parents(nC, 8))
    return resF, parents
end

# Return the exact volume of the obstacles and the share that the mesh
# captures. The share is the volume that the mesh removes from the box,
# divided by the exact obstacle volume. A rectangle or a box gives 1. A
# circle or a sphere gives less than 1, because the mesh holds a polygon
# or a polyhedron inside the shape. The record fields follow
# _obstacle_record2 and _obstacle_record3.
function _obstacle_capture(obstacles, box::Tuple, dim::Int, totalVolume::Float64)
    exact = 0.0
    for obs in obstacles
        if obs.shape == "circle"
            exact += pi * obs.radius^2
        elseif obs.shape == "sphere"
            exact += 4.0 / 3.0 * pi * obs.radius^3
        else
            exact += 2.0^dim * prod(obs.half_extents)
        end
    end
    boxVolume = prod(box[1:dim])
    return exact, (boxVolume - totalVolume) / exact
end

# Check the config value of the key hierarchy. Return nothing when the
# key is absent. Otherwise return a NamedTuple with the fields levels,
# boundary and maxCells.
function _resolve_hierarchy_options(raw)
    raw === nothing && return nothing
    raw isa AbstractDict || error("'hierarchy' must be a mapping")
    known = ("levels", "boundary", "max_cells")
    for k in keys(raw)
        k in known || error("hierarchy: unknown key '$k'. Valid keys: $(join(known, ", "))")
    end
    haskey(raw, "levels") || error("hierarchy: 'levels' is required")
    whole(v) = v isa Integer && !(v isa Bool)
    levels = raw["levels"]
    (whole(levels) && levels >= 1) || error("hierarchy['levels'] must be a whole number, 1 or more, got $(repr(levels))")
    boundary = get(raw, "boundary", "frozen")
    if boundary == "projected"
        error("hierarchy['boundary'] 'projected' is not available yet. Use 'frozen'.")
    end
    boundary == "frozen" || error("hierarchy['boundary'] must be 'frozen', got $(repr(boundary))")
    maxCells = get(raw, "max_cells", 2000000)
    (whole(maxCells) && maxCells >= 1) || error("hierarchy['max_cells'] must be a whole number, 1 or more, got $(repr(maxCells))")
    return (levels = Int(levels), boundary = String(boundary), maxCells = Int(maxCells))
end

# ── Config check ─────────────────────────────────────────────────────────
#
# This function needs no package. It takes the Dict that the YAML reader
# gives. It sets the defaults, checks every key against the mesh type,
# and changes the values to Tuples. It returns the same Dict.

# Return a short text with the valid values of mesh_type.
function _mesh_type_list()
    return join(_mesh_types(), ", ")
end

function _check_config(cfg)
    cfg isa AbstractDict || error("the config must be a YAML mapping")
    cfg isa Dict || (cfg = Dict{Any,Any}(cfg))

    # ---- Mesh type.
    mt = get(cfg, "mesh_type", nothing)
    if mt === nothing
        if haskey(cfg, "mesh_options")
            error("'mesh_options' is present but 'mesh_type' is not set. " *
                  "Set mesh_type to one of $(_mesh_type_list()).")
        end
        mt = "gridsplit"
    end
    if !(mt isa AbstractString) || !(mt in _mesh_types())
        error("mesh_type must be one of $(_mesh_type_list()), got '$mt'")
    end
    dim = (mt == "kuhn" || mt == "tetra") ? 3 : 2
    structured = mt in _structured_types()

    # ---- Keys that no type accepts, and keys that one family rejects.
    haskey(cfg, "gmsh_options") &&
        error("'gmsh_options' is not valid for mesh_type='$mt'. This file does not use Gmsh.")
    if structured
        haskey(cfg, "mesh_options") &&
            error("'mesh_options' is not valid for mesh_type='$mt'")
        obs = get(cfg, "obstacles", nothing)
        if obs !== nothing && !isempty(obs)
            error("'obstacles' is not valid for mesh_type='$mt'. A structured grid has no obstacles.")
        end
        if mt == "equilateral"
            for key in ("diagonal", "stretching_ratio", "stretching_symmetric")
                haskey(cfg, key) && error("'$key' is not valid for mesh_type='$mt'")
            end
        end
    else
        for key in ("n_divisions", "diagonal", "stretching_ratio", "stretching_symmetric")
            haskey(cfg, key) && error("'$key' is not valid for mesh_type='$mt'")
        end
        haskey(cfg, "mesh_options") ||
            error("'mesh_options' is required for mesh_type='$mt'. It needs at least 'mesh_size'.")
    end

    # ---- Required keys and defaults.
    haskey(cfg, "box_length") || error("'box_length' is required")
    haskey(cfg, "periodic") || (cfg["periodic"] = fill(false, dim))
    (!haskey(cfg, "obstacles") || cfg["obstacles"] === nothing) && (cfg["obstacles"] = [])
    if structured
        haskey(cfg, "n_divisions") || error("'n_divisions' is required for mesh_type='$mt'")
    end
    keysToMeasure = structured ? ("box_length", "periodic", "n_divisions") : ("box_length", "periodic")
    for key in keysToMeasure
        if length(cfg[key]) != dim
            error("'$key' needs $dim entries for mesh_type='$mt'")
        end
    end
    if mt == "gridsplit" || mt == "kuhn"
        haskey(cfg, "diagonal") || (cfg["diagonal"] = mt == "gridsplit" ? "backslash" : "main")
        haskey(cfg, "stretching_ratio") || (cfg["stretching_ratio"] = fill(1.0, dim))
        haskey(cfg, "stretching_symmetric") || (cfg["stretching_symmetric"] = fill(false, dim))
        if length(cfg["stretching_ratio"]) != dim || length(cfg["stretching_symmetric"]) != dim
            error("'stretching_ratio' and 'stretching_symmetric' need $dim entries")
        end
    end

    # ---- Convert the values.
    cfg["box_length"] = Tuple(Float64.(cfg["box_length"]))
    cfg["periodic"] = Tuple(Bool.(cfg["periodic"]))
    if structured
        cfg["n_divisions"] = Tuple(Int.(cfg["n_divisions"]))
    elseif mt == "triangle"
        cfg["mesh_options"] = _resolve_mesh_options2(cfg["mesh_options"])
    else
        cfg["mesh_options"] = _resolve_mesh_options3(cfg["mesh_options"])
    end
    cfg["hierarchy"] = _resolve_hierarchy_options(get(cfg, "hierarchy", nothing))
    cfg["mesh_type"] = mt
    cfg["dim"] = dim
    return cfg
end

# ── Periodicity ───────────────────────────────────────────────────────────
#
# These functions use coordinates, not indices. They find the pairs of
# boundary vertices that lie on opposite box sides or faces. The vertex
# array has 3 columns for every mesh type. A 2-D mesh has a zero z column.

# Find the pairs of vertices on opposite sides of the box. Return a Dict
# from the axis name to an array of 0-based vertex pairs. The function
# puts the vertices of the high side into a hash table. The key is the
# pair of the other two coordinates, rounded to the match tolerance.
# Then the function looks up each vertex of the low side in the 9 nearest
# keys. The time is proportional to the number of vertices. The result
# equals the result of a search over all vertices. Two points that are
# closer than matchTol have keys that differ by 1 at most.
function _periodic_pairs(logBuf::IOBuffer, pts::Matrix{Float64}, box::Tuple, periodic::Tuple; bndTol = 1e-9, matchTol = 1e-7)
    pairs = Dict{String,Matrix{Int}}()
    for d in 1:length(periodic)
        periodic[d] || continue
        L, dname = box[d], _dirs()[d]
        other = [k for k in 1:size(pts, 2) if k != d]
        lo = [i for i in 1:size(pts, 1) if abs(pts[i, d]) < bndTol]
        hi = [i for i in 1:size(pts, 1) if abs(pts[i, d] - L) < bndTol]
        if isempty(lo) || isempty(hi)
            _log(logBuf, "  [warn] periodic_$dname: no boundary vertices found")
            pairs[dname] = zeros(Int, 0, 2)
            continue
        end
        length(lo) != length(hi) && _log(logBuf, "  [warn] periodic_$dname: $(length(lo)) lo vs $(length(hi)) hi vertices")
        table = Dict{Tuple{Int,Int},Vector{Int}}()
        for vj in hi
            key = (round(Int, pts[vj, other[1]] / matchTol), round(Int, pts[vj, other[2]] / matchTol))
            push!(get!(table, key, Int[]), vj)
        end
        matched = Tuple{Int,Int}[]
        nMiss = 0
        for vi in lo
            k1 = round(Int, pts[vi, other[1]] / matchTol)
            k2 = round(Int, pts[vi, other[2]] / matchTol)
            bestJ, bestD = 0, Inf
            for d1 in -1:1, d2 in -1:1
                for vj in get(table, (k1 + d1, k2 + d2), Int[])
                    dist = hypot(pts[vi, other[1]] - pts[vj, other[1]], pts[vi, other[2]] - pts[vj, other[2]])
                    dist < bestD && (bestD = dist; bestJ = vj)
                end
            end
            bestD < matchTol ? push!(matched, (vi - 1, bestJ - 1)) : (nMiss += 1)
        end
        pairs[dname] = isempty(matched) ? zeros(Int, 0, 2) :
            reduce(vcat, [reshape([a, b], 1, 2) for (a, b) in matched])
        status = "$(length(matched)) pairs" * (nMiss > 0 ? ", $nMiss unmatched" : "")
        _log(logBuf, "  periodic_$dname: $status")
    end
    return pairs
end
function _periodic_groups(perPairs::Dict{String,Matrix{Int}})
    parent = Dict{Int,Int}()
    function find(x::Int)
        haskey(parent, x) || (parent[x] = x)
        root = x
        while parent[root] != root
            root = parent[root]
        end
        while parent[x] != root
            parent[x], x = root, parent[x]
        end
        return root
    end
    for arr in values(perPairs), r in 1:size(arr, 1)
        ra, rb = find(arr[r, 1]), find(arr[r, 2])
        ra != rb && (parent[rb] = ra)
    end
    groups = Dict{Int,Vector{Int}}()
    for v in keys(parent)
        push!(get!(groups, find(v), Int[]), v)
    end
    return groups
end

function _all_periodic_pairs(perPairs::Dict{String,Matrix{Int}})
    groups = _periodic_groups(perPairs)
    allPairs = Tuple{Int,Int}[]
    for members in values(groups)
        length(members) < 2 && continue
        s = sort(members)
        for i in 1:length(s)-1, j in (i+1):length(s)
            push!(allPairs, (s[i], s[j]))
        end
    end
    isempty(allPairs) && return zeros(Int, 0, 2)
    out = reduce(vcat, [reshape([a, b], 1, 2) for (a, b) in allPairs])
    return out[sortperm(1:size(out, 1), by = i -> (out[i, 1], out[i, 2])), :]
end

function _enforce_periodicity_vertex_field(vertexField::Vector{Float64}, perPairs::Dict{String,Matrix{Int}})
    groups = _periodic_groups(perPairs)
    isempty(groups) && return vertexField
    merged = copy(vertexField)
    for members in values(groups)
        total = sum(vertexField[v + 1] for v in members)
        for v in members
            merged[v + 1] = total
        end
    end
    return merged
end

# =====================================================================
# package-free part: ends
# =====================================================================

# ── Configuration ────────────────────────────────────────────────────────

# Read one YAML config. Set the defaults. Check the keys.
function load_config(path::AbstractString)
    return _check_config(YAML.load_file(path))
end

# ── HDF5 and XDMF output ──────────────────────────────────────────────────

# Write one dataset. Skip an empty array.
function _write_ds(loc, name::AbstractString, data)
    length(data) == 0 && return
    # Julia arrays are column-major and HDF5 datasets are row-major, so
    # HDF5.jl reverses the axes of what it writes. Reverse them here once
    # more. The dataset then has the shape (N, k) that XDMF, h5py and
    # ParaView expect, and _read_ds undoes it for the Julia kernels.
    out = ndims(data) >= 2 ? permutedims(data, ndims(data):-1:1) : data
    d = create_dataset(loc, name, datatype(out), dataspace(out); chunk = size(out), deflate = 4)
    write(d, out)
end

# Read one dataset back as the Julia array that _write_ds got. A matrix
# has the shape (N, k) again. The axes are permuted once, at load time.
function _read_ds(loc, name::AbstractString)
    data = read(loc[name])
    return ndims(data) >= 2 ? permutedims(data, ndims(data):-1:1) : data
end

# Read one mesh file that either writer made. Return a NamedTuple with the
# arrays that the kernels use. A matrix has the shape (N, k). Every number
# of a vertex, cell or edge is 1-based, as in the kernels. The file keeps
# them 0-based. The fields are:
# dim, meshType, boxLength, schemaVersion: the scalar attributes.
# quality: a Dict of the quality attributes.
# pts: the points, (N, 3).
# cells: the vertices of each cell, (nCells, dim+1).
# edges: the vertices of each edge, (nEdges, 2).
# c2e: the edges of each cell, (nCells, 3 or 6).
# cellVol, vertexVol: the volumes. cvn: the vertex normals, (nCells, dim+1, 3).
# v2cOff, v2cIdx: the cells of vertex v are v2cIdx[v2cOff[v]:v2cOff[v+1]-1].
# ptMk: the vertex markers. edgeMk: the edge markers of a 2-D mesh, else nothing.
# periodic: one Bool per axis. perPairs: a Dict from x, y, z to a pair matrix.
# allPairs: all pairs, (nPairs, 2). obstacles: a list of NamedTuples.
# hierarchy: nothing, or the fields level, nLevels, refinement, boundary,
# childrenPerCell, vertexParents and cellParent. The last two are nothing
# at level 0. Level 0 is the coarsest.
# The function stops with an error when schema_version is not 1.
function load_mesh(path::AbstractString)
    h5open(path, "r") do f
        a = attrs(f)
        haskey(a, "schema_version") || error("$path has no schema_version attribute. It is not a mesh file of this writer.")
        version = Int(a["schema_version"])
        version == 1 || error("$path has schema_version $version. load_mesh reads version 1 only.")
        dim = Int(a["dimension"])
        quality = Dict{String,Float64}()
        for key in ("min_angle_deg", "max_angle_deg", "min_dihedral_deg", "max_dihedral_deg", "max_radius_edge_ratio")
            haskey(a, key) && (quality[key] = Float64(a[key]))
        end
        rd(g, name) = _read_ds(g, name)
        pts = rd(f, "points")
        cells = Int.(rd(f, "cell_to_vertex")) .+ 1
        edges = Int.(rd(f, "vertex_to_vertex")) .+ 1
        c2e = Int.(rd(f, "cell_to_edge")) .+ 1
        v2cOff = Int.(rd(f["vertex_to_cell"], "offsets")) .+ 1
        v2cIdx = Int.(rd(f["vertex_to_cell"], "indices")) .+ 1
        edgeMk = haskey(f, "boundary/edge_markers") ? rd(f["boundary"], "edge_markers") : nothing

        p = f["periodicity"]
        pa = attrs(p)
        periodic = Tuple(haskey(pa, "periodic_$d") && pa["periodic_$d"] != 0 for d in _dirs()[1:dim])
        perPairs = Dict{String,Matrix{Int}}()
        for d in _dirs()[1:dim]
            haskey(p, "$(d)_pairs") && (perPairs[d] = Int.(rd(p, "$(d)_pairs")) .+ 1)
        end
        allPairs = haskey(p, "all_pairs") ? Int.(rd(p, "all_pairs")) .+ 1 : zeros(Int, 0, 2)

        obstacles = NamedTuple[]
        if haskey(f, "obstacles")
            o = f["obstacles"]
            for i in 0:Int(attrs(o)["count"]) - 1
                og = o["obstacle_$i"]
                oa = attrs(og)
                shape = String(oa["shape"])
                extent = haskey(oa, "radius") ? Float64(oa["radius"]) : Float64.(oa["half_extents"])
                push!(obstacles, (shape = shape, marker = Int(oa["marker"]),
                                  center = Float64.(oa["center"]), extent = extent))
            end
        end
        hierarchy = nothing
        if haskey(f, "hierarchy")
            h = f["hierarchy"]
            ha = attrs(h)
            level = Int(ha["level"])
            hierarchy = (level = level, nLevels = Int(ha["n_levels"]), refinement = String(ha["refinement"]),
                         boundary = String(ha["boundary"]), childrenPerCell = Int(ha["children_per_cell"]),
                         vertexParents = level > 0 ? Int.(rd(h, "fine_vertex_parents")) .+ 1 : nothing,
                         cellParent = level > 0 ? Int.(rd(h, "cell_parent")) .+ 1 : nothing)
        end
        return (dim = dim, meshType = String(a["mesh_type"]), boxLength = Float64.(a["box_length"]),
                schemaVersion = version, quality = quality,
                pts = pts, cells = cells, edges = edges, c2e = c2e,
                cellVol = rd(f, "cell_vol"), cvn = rd(f, "cell_vertex_normal"), vertexVol = rd(f, "vertex_vol"),
                v2cOff = v2cOff, v2cIdx = v2cIdx, ptMk = rd(f["boundary"], "point_markers"), edgeMk = edgeMk,
                periodic = periodic, perPairs = perPairs, allPairs = allPairs, obstacles = obstacles,
                hierarchy = hierarchy)
    end
end

# Read all levels of a hierarchy that main wrote. The argument is the path
# of the YAML config, as for main. The function returns a vector of the
# load_mesh results, with level 0 (the coarsest) first. It stops with an
# error when a level file is missing or does not belong to the hierarchy.
function load_hierarchy(yamlPath::AbstractString)
    stem = replace(yamlPath, r"\.ya?ml$" => "")
    file0 = "$(stem)_L0.h5"
    isfile(file0) || error("$file0 not found. Run writer.jl on a config with the key hierarchy first.")
    m0 = load_mesh(file0)
    m0.hierarchy === nothing && error("$file0 has no hierarchy group.")
    nLevels = m0.hierarchy.nLevels
    levels = NamedTuple[m0]
    for k in 1:nLevels - 1
        file = "$(stem)_L$k.h5"
        isfile(file) || error("$file not found. The hierarchy has $nLevels levels.")
        m = load_mesh(file)
        if m.hierarchy === nothing || m.hierarchy.level != k || m.hierarchy.nLevels != nLevels
            error("$file does not belong to the hierarchy of $file0.")
        end
        push!(levels, m)
    end
    return levels
end

# Return the value of the HDF5 attribute mesh_type. It names the cell
# type, as writer.py does: triangle for every 2-D mesh and tetra for
# every 3-D mesh. The config key mesh_type records how the mesh was made.
function _h5_mesh_type(meshType::AbstractString)
    (meshType == "kuhn" || meshType == "tetra") && return "tetra"
    return "triangle"
end

# Return the name of the top XDMF grid. It is the same for every mesh
# type, as in writer.py.
function _xdmf_grid_name(meshType::AbstractString)
    return "MeshCollection"
end

# Format a number with 6 significant digits, for log lines.
function _g6(x::Float64)
    return string(round(x; sigdigits = 6))
end

# Compute the quality figures of the mesh. Return a dictionary of
# attribute names and values. A 2-D mesh gets the smallest and the largest
# triangle angle. A 3-D mesh gets the smallest and the largest dihedral
# angle and the largest radius-edge ratio. The angles are in degrees.
function _cell_quality(cells::Matrix{Int}, pts3d::Matrix{Float64}, dim::Int)
    nCells = size(cells, 1)
    q = Dict{String,Float64}()
    if dim == 2
        minAng = zeros(Float64, nCells)
        maxAng = zeros(Float64, nCells)
        triangle_angle_kernel(cells, pts3d, nCells, minAng, maxAng)
        q["min_angle_deg"] = minimum(minAng)
        q["max_angle_deg"] = maximum(maxAng)
    else
        # The dihedral kernel works in radians. The radius-edge ratio comes
        # from the circumcenter kernel and the quality kernel.
        minD = zeros(Float64, nCells)
        maxD = zeros(Float64, nCells)
        tet_dihedral_kernel(cells, pts3d, nCells, minD, maxD)
        ids = collect(1:nCells)
        vx = pts3d[:, 1]; vy = pts3d[:, 2]; vz = pts3d[:, 3]
        ccx = zeros(Float64, nCells); ccy = zeros(Float64, nCells); ccz = zeros(Float64, nCells)
        ccR = zeros(Float64, nCells)
        cell_circumcenter_3d_kernel(cells, vx, vy, vz, ids, nCells, ccx, ccy, ccz, ccR)
        hAt = ones(Float64, nCells)
        ratio = zeros(Float64, nCells)
        sizeRatio = zeros(Float64, nCells)
        cell_quality_3d_kernel(cells, vx, vy, vz, ids, nCells, ccR, hAt, ratio, sizeRatio)
        q["min_dihedral_deg"] = rad2deg(minimum(minD))
        q["max_dihedral_deg"] = rad2deg(maximum(maxD))
        q["max_radius_edge_ratio"] = maximum(ratio)
    end
    return q
end

# Write one mesh with all its arrays to one HDF5 file.
function write_hdf5(logBuf::IOBuffer, h5Path, meshType, pts, cells, edges, c2e, cellVol, cvn, vertexVol,
                    v2cOff, v2cIdx, ptMk, edgeMk, quality::Dict{String,Float64},
                    periodic::Tuple, perPairs::Dict, allPairs, box::Tuple, dim::Int,
                    obstacles::Vector, hier)
    h5open(h5Path, "w") do f
        attrs(f)["schema_version"] = 1
        attrs(f)["mesh_type"] = _h5_mesh_type(meshType)
        attrs(f)["dimension"] = dim
        attrs(f)["box_length"] = collect(box)
        attrs(f)["n_vertices"] = size(pts, 1)
        attrs(f)["n_cells"] = size(cells, 1)
        attrs(f)["n_edges"] = size(edges, 1)
        for key in sort(collect(keys(quality)))
            attrs(f)[key] = quality[key]
        end

        _write_ds(f, "points", pts)
        _write_ds(f, "cell_to_vertex", cells)
        _write_ds(f, "cell_to_edge", c2e)
        _write_ds(f, "cell_vol", cellVol)
        _write_ds(f, "cell_vertex_normal", cvn)
        _write_ds(f, "vertex_to_vertex", edges)

        g = create_group(f, "vertex_to_cell")
        _write_ds(g, "offsets", v2cOff)
        _write_ds(g, "indices", v2cIdx)
        _write_ds(f, "vertex_vol", vertexVol)

        b = create_group(f, "boundary")
        _write_ds(b, "point_markers", ptMk)
        edgeMk !== nothing && _write_ds(b, "edge_markers", edgeMk)

        p = create_group(f, "periodicity")
        for (d, dname) in enumerate(_dirs()[1:dim])
            attrs(p)["periodic_$dname"] = periodic[d]
        end
        for (dname, arr) in perPairs
            attrs(p)["periodic_$dname"] = true
            _write_ds(p, "$(dname)_pairs", arr)
        end
        _write_ds(p, "all_pairs", allPairs)

        if hier !== nothing
            h = create_group(f, "hierarchy")
            attrs(h)["level"] = hier.level
            attrs(h)["n_levels"] = hier.nLevels
            attrs(h)["refinement"] = "red"
            attrs(h)["boundary"] = hier.boundary
            attrs(h)["children_per_cell"] = hier.children
            if hier.level > 0
                _write_ds(h, "fine_vertex_parents", hier.vertexParents .- 1)
                _write_ds(h, "cell_parent", hier.cellParent .- 1)
            end
        end

        if !isempty(obstacles)
            o = create_group(f, "obstacles")
            attrs(o)["count"] = length(obstacles)
            for (i, obs) in enumerate(obstacles)
                og = create_group(o, "obstacle_$(i-1)")
                attrs(og)["shape"] = obs.shape
                attrs(og)["marker"] = obs.marker
                attrs(og)["center"] = collect(obs.center)
                if obs.shape == "circle" || obs.shape == "sphere"
                    attrs(og)["radius"] = obs.radius
                else
                    attrs(og)["half_extents"] = collect(obs.half_extents)
                end
            end
        end
    end
    _log(logBuf, "")
    _log(logBuf, "  HDF5 written -> $h5Path")
    _log(logBuf, "    vertices          : $(size(pts,1))")
    _log(logBuf, "    cells             : $(size(cells,1))")
    _log(logBuf, "    edges             : $(size(edges,1))")
    if haskey(quality, "min_angle_deg")
        _log(logBuf, "    quality           : angle min=$(_g6(quality["min_angle_deg"])) max=$(_g6(quality["max_angle_deg"])) deg")
    else
        _log(logBuf, "    quality           : dihedral min=$(_g6(quality["min_dihedral_deg"])) max=$(_g6(quality["max_dihedral_deg"])) deg, " *
                     "radius-edge max=$(_g6(quality["max_radius_edge_ratio"]))")
    end
end

# Write one XDMF file that points at the mesh data in one HDF5 file.
function write_xdmf(logBuf::IOBuffer, xdmfPath, h5Path, nPts::Int, nCells::Int, nEdges::Int, dim::Int, meshType::AbstractString)
    h5name = basename(h5Path)
    tt = dim == 2 ? "Triangle" : "Tetrahedron"
    nn = dim == 2 ? 3 : 4
    geo = "$nPts 3"
    gridName = _xdmf_grid_name(meshType)
    function di(dims, dset; dtype = "Int", prec = nothing)
        pr = prec === nothing ? "" : " Precision=\"$prec\""
        return "<DataItem Format=\"HDF\" DataType=\"$dtype\"$pr Dimensions=\"$dims\">$h5name:$dset</DataItem>"
    end
    cellGrid = """
        <Grid Name="cells" GridType="Uniform">
          <Topology TopologyType="$tt" NumberOfElements="$nCells">
            $(di("$nCells $nn", "/cell_to_vertex"))
          </Topology>
          <Geometry GeometryType="XYZ">
            $(di(geo, "/points"; dtype = "Float", prec = 8))
          </Geometry>
          <Attribute Name="BoundaryMarker" AttributeType="Scalar" Center="Node">
            $(di("$nPts", "/boundary/point_markers"))
          </Attribute>
          <Attribute Name="VertexVolume" AttributeType="Scalar" Center="Node">
            $(di("$nPts", "/vertex_vol"; dtype = "Float", prec = 8))
          </Attribute>
          <Attribute Name="CellVolume" AttributeType="Scalar" Center="Cell">
            $(di("$nCells", "/cell_vol"; dtype = "Float", prec = 8))
          </Attribute>
        </Grid>"""
    edgeGrid = """
        <Grid Name="edges" GridType="Uniform">
          <Topology TopologyType="Polyline" NumberOfElements="$nEdges"
                    NodesPerElement="2">
            $(di("$nEdges 2", "/vertex_to_vertex"))
          </Topology>
          <Geometry GeometryType="XYZ">
            $(di(geo, "/points"; dtype = "Float", prec = 8))
          </Geometry>
        </Grid>"""
    xdmf = """<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf Version="3.0">
  <Domain>
    <Grid Name="$gridName" GridType="Collection"
          CollectionType="Spatial">
$cellGrid
<!-- $edgeGrid -->
    </Grid>
  </Domain>
</Xdmf>
"""
    write(xdmfPath, xdmf)
    _log(logBuf, "  XDMF written -> $xdmfPath")
end

# ── Provenance ────────────────────────────────────────────────────────────

# Make one string safe to place inside an XML comment.
function _xml_comment_safe(text::AbstractString)
    safe = text
    while occursin("--", safe)
        safe = replace(safe, "--" => "- -")
    end
    endswith(safe, "-") && (safe *= " ")
    return safe
end

# Append three record blocks to the end of one XDMF file.
function _append_provenance(logBuf::IOBuffer, xdmfPath, yamlPath, yamlText, terminalText, scriptText)
    block(label, content) = "\n\n<!-- ===== $label =====\n$(_xml_comment_safe(content))\n===== end $label ===== -->\n"
    open(xdmfPath, "a") do f
        write(f, "\n")
        write(f, block("YAML CONFIG: $(basename(yamlPath))", yamlText))
        write(f, block("TERMINAL OUTPUT", terminalText))
        write(f, block("WRITER.JL SOURCE", scriptText))
    end
    _log(logBuf, "  Provenance appended -> $xdmfPath")
end

# ── Entry point ───────────────────────────────────────────────────────────

# Build one mesh from one YAML config. Write the output files. Every mesh
# type goes through the same steps after the builder, in _write_level.
function main(yamlPath::AbstractString)
    logBuf = IOBuffer()
    cfg = load_config(yamlPath)
    _log(logBuf, "Building '$(cfg["mesh_type"])' mesh …")
    res = _build_mesh(logBuf, cfg)
    stem = replace(yamlPath, r"\.ya?ml$" => "")
    if cfg["hierarchy"] === nothing
        _write_level(logBuf, yamlPath, stem, cfg, res, nothing)
    else
        _write_hierarchy(logBuf, yamlPath, stem, cfg, res)
    end
    return nothing
end

# Write the levels of a hierarchy. Level 0 is the mesh of the builder.
# Each next level is the red refinement of the one before. The files are
# stem_L0, stem_L1, and so on. The function stops with an error before
# the refinement when the finest level would have more cells than
# max_cells.
function _write_hierarchy(logBuf::IOBuffer, yamlPath::AbstractString, stem::AbstractString, cfg::Dict, res)
    H = cfg["hierarchy"]
    dim = cfg["dim"]
    children = dim == 2 ? 4 : 8
    nCoarse = size(res.cells, 1)
    estimate = Float64(nCoarse) * Float64(children)^(H.levels - 1)
    if estimate > H.maxCells
        error("hierarchy: the finest level of $(H.levels) levels would have about $(round(Int, estimate)) cells. " *
              "The limit max_cells is $(H.maxCells). The coarsest mesh has $nCoarse cells. " *
              "Use fewer levels, a coarser mesh, or a larger max_cells.")
    end

    # Level 0 needs the edge array and the edge markers before it refines.
    if res.edges === nothing
        edges, c2e = _build_edge_registry(res.cells::Matrix{Int}, dim)
        res = merge(res, (edges = edges, c2e = c2e))
    end
    res = merge(res, (edgeMk = _grid_edge_markers(res, res.edges::Matrix{Int}),))

    parents = nothing
    for k in 0:H.levels - 1
        buf = k == 0 ? logBuf : IOBuffer()
        if k > 0
            _log(buf, "Refining level $(k - 1) to make level $k …")
            res, parents = dim == 2 ? _refine_mesh2(res) : _refine_mesh3(res)
        end
        _log(buf, "Level $k of $(H.levels - 1): $(size(res.pts, 1)) vertices, $(size(res.cells, 1)) cells")
        hier = (level = k, nLevels = H.levels, boundary = H.boundary, children = children,
                vertexParents = parents === nothing ? nothing : parents.vertexParents,
                cellParent = parents === nothing ? nothing : parents.cellParent)
        _write_level(buf, yamlPath, "$(stem)_L$k", cfg, res, hier)
    end
    return nothing
end

# Run the steps after the builder on one result record. The function
# computes the edges, the volumes, the normals, the quality, the
# vertex-to-cell lists and the periodic pairs. Then it writes the HDF5
# file and the XDMF file. The path stem has no extension. At the end it
# appends the provenance record. yamlPath names the config for it.
# hier is nothing for a mesh without hierarchy. Otherwise it is the
# NamedTuple of one level. A level after level 0 carries a pointer to
# level 0 instead of the writer source.
function _write_level(logBuf::IOBuffer, yamlPath::AbstractString, stem::AbstractString, cfg::Dict, res, hier)
    mt, dim, per, box = cfg["mesh_type"], cfg["dim"], cfg["periodic"], cfg["box_length"]
    cells = res.cells::Matrix{Int}
    resPts = res.pts::Matrix{Float64}
    ptMk = res.ptMk::Vector{Int32}
    nPts = size(resPts, 1)
    nCells = size(cells, 1)

    # Points and periodicity work in 3 columns throughout. A 2-D mesh
    # carries a zero z column, so the same writers serve both.
    pts3d = dim == 2 ? hcat(resPts, zeros(Float64, nPts)) : resPts

    # The edge registry is built once. The Delaunay builders have it
    # already, because their smoothing needs it.
    _log(logBuf, "Computing cell-to-edge connectivity …")
    if res.edges === nothing
        edges, c2e = _build_edge_registry(cells, dim)
    else
        edges, c2e = res.edges::Matrix{Int}, res.c2e::Matrix{Int}
    end
    edgeMk = _grid_edge_markers(res, edges)

    _log(logBuf, "Computing cell volumes …")
    cellVol = zeros(Float64, nCells)
    if dim == 2
        signed_measure_2d_kernel(cells, pts3d, nCells, cellVol)
    else
        signed_measure_3d_kernel(cells, pts3d, nCells, cellVol)
    end
    for c in 1:nCells
        cellVol[c] = abs(cellVol[c])
    end
    if !isempty(res.obstacles)
        exactObs, share = _obstacle_capture(res.obstacles, box, dim, sum(cellVol))
        _log(logBuf, "  obstacle volume captured by the mesh: $(round(100 * share, digits = 2)) percent of $(round(exactObs, sigdigits = 6))")
        if hier !== nothing && hier.level == 0 && share < 0.95
            _log(logBuf, "  [warn] hierarchy: the coarsest mesh captures $(round(100 * share, digits = 1)) percent of the obstacle volume. " *
                         "Every level keeps this boundary. Use a smaller mesh size on the obstacles.")
        end
    end

    _log(logBuf, "Computing per-cell vertex normals …")
    cvn = zeros(Float64, nCells, dim + 1, 3)
    if dim == 2
        cell_vertex_normal_2d_kernel(cells, pts3d, nCells, cvn)
    else
        cell_vertex_normal_3d_kernel(cells, pts3d, nCells, cvn)
    end

    _log(logBuf, "Computing cell quality …")
    quality = _cell_quality(cells, pts3d, dim)

    _log(logBuf, "Computing vertex volumes …")
    vertexVol = zeros(Float64, nPts)
    vertex_volume_kernel(cells, cellVol, nCells, dim + 1, vertexVol)

    _log(logBuf, "Computing vertex-to-cell connectivity …")
    v2cOff, v2cIdx = _vertex_to_cell(nPts, cells)

    _log(logBuf, "Detecting periodic vertex pairs …")
    perPairs = _periodic_pairs(logBuf, pts3d, res.boxSize, per)
    allPairs = _all_periodic_pairs(perPairs)

    _log(logBuf, "Enforcing periodicity on vertex volumes …")
    vertexVol = _enforce_periodicity_vertex_field(vertexVol, perPairs)

    # Every array above holds 1-based row numbers. Only the stored
    # connectivity becomes 0-based, here, to match the XDMF file.
    cellsOut = cells .- 1
    edgesOut = edges .- 1

    hdf5Path = stem * ".h5"
    xdmfPath = stem * ".xdmf"
    write_hdf5(logBuf, hdf5Path, mt, pts3d, cellsOut, edgesOut, c2e, cellVol, cvn, vertexVol,
               v2cOff, v2cIdx, ptMk, edgeMk, quality, per, perPairs, allPairs, box, dim,
               res.obstacles, hier)
    write_xdmf(logBuf, xdmfPath, hdf5Path, nPts, nCells, size(edges, 1), dim, mt)

    terminalText = String(take!(logBuf))
    if hier !== nothing && hier.level > 0
        scriptText = "The writer source is in $(replace(basename(stem), r"_L[0-9]+$" => "_L0")).xdmf."
    else
        scriptText = read(@__FILE__, String)
    end
    _append_provenance(logBuf, xdmfPath, yamlPath, read(yamlPath, String), terminalText, scriptText)
end

# Run main on one YAML file, or on every YAML file in one folder.
function _process_yaml_or_dir(pathStr::AbstractString)
    if isdir(pathStr)
        yamlFiles = sort(filter(f -> endswith(f, ".yaml"), readdir(pathStr, join = true)))
        if isempty(yamlFiles)
            println(stderr, "No .yaml files found in directory: $pathStr")
            exit(1)
        end
        println("Found $(length(yamlFiles)) YAML config(s) in $pathStr\n")
        failures = String[]
        for (i, yf) in enumerate(yamlFiles)
            println("[$i/$(length(yamlFiles))] $(basename(yf))")
            try
                main(yf)
            catch e
                println("  [error] $(basename(yf)): $(sprint(showerror, e))")
                push!(failures, basename(yf))
            end
            println()
        end
        if !isempty(failures)
            println(stderr, "Completed with $(length(failures)) failure(s): $(join(failures, ", "))")
            exit(1)
        end
        return
    end
    if !isfile(pathStr)
        println(stderr, "Path not found: $pathStr")
        exit(1)
    end
    main(pathStr)
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 1
        println(stderr, "Usage: julia $(PROGRAM_FILE) config.yaml | config_dir/")
        exit(1)
    end
    _process_yaml_or_dir(ARGS[1])
end
