#!/usr/bin/env julia
# Check a mesh hierarchy of writer.jl with a small geometric multigrid
# solver. The script reads the files that writer.jl wrote, the same way a
# solver does. It uses only load_hierarchy and the arrays in the files.
#
# For each hierarchy the script does these checks:
#   1. The prolongation matrix P, built from fine_vertex_parents, copies
#      every linear function of the coarse level to the fine level.
#   2. On a periodic mesh, a periodic function of the coarse level is
#      still periodic after P. The two periodic images of a fine vertex
#      must get the same value.
#   3. The Galerkin product P' * A_fine * P equals the stiffness matrix
#      A_coarse, which the script assembles on the coarse mesh. This holds
#      only when the coarse space lies inside the fine space, so it tests
#      the nesting, the vertex order and the cell order in one number.
#   4. A fine vertex on the boundary has no free coarse vertex as parent.
#      Without this, the zero values on the boundary would not nest.
#   5. A V(2,2)-cycle solves a Poisson problem. The script prints the
#      convergence factor for 2, 3, ... levels. A factor of 0.7 or more
#      counts as a failure. A broken hierarchy gives a factor near 1.
# The problem has a zero value at every vertex on an obstacle and on a box
# side without periodicity. The vertices that are periodic images of each
# other are one unknown. A mesh that is periodic in every direction and has
# no obstacle has no zero value. The script then solves -Laplace(u) + u = 1
# with a lumped mass matrix, because the pure Laplace matrix is singular. The smoother is a plain
# Gauss-Seidel sweep, so the factor depends on the mesh quality. It can
# grow a little with the number of levels on unstructured 3-D meshes.
# Usage: julia --project=. tools/hierarchy_multigrid_check.jl [NAME ...]
# NAME is one of the names that _configs() returns. Without a name, the script runs all.
using SparseArrays, LinearAlgebra, Random
include(joinpath(@__DIR__, "..", "writer.jl"))

# ── Kernels ──────────────────────────────────────────────────────────────

# p1_gradients_kernel(cells, pts, nc, dim, grad, vol)
#
# Compute the gradients of the linear basis functions of each cell, and
# the cell volume. A 2-D cell has 3 vertices, a 3-D cell has 4.
#
# cells: vertex numbers of each cell, one row per cell
# pts: vertex coordinates, one row per vertex, 3 columns
# nc: number of cells
# dim: 2 or 3
# grad: output, nc by nloc * dim. Entry (c, (a - 1) * dim + d) is
#   component d of the gradient of basis function a in cell c
# vol: output, volume of each cell
function p1_gradients_kernel(cells, pts, nc, dim, grad, vol)
    for c = 1:nc
        v0 = cells[c, 1]
        v1 = cells[c, 2]
        v2 = cells[c, 3]
        ax = pts[v1, 1] - pts[v0, 1]
        ay = pts[v1, 2] - pts[v0, 2]
        bx = pts[v2, 1] - pts[v0, 1]
        by = pts[v2, 2] - pts[v0, 2]
        if dim == 2
            # Edge vectors a and b from vertex 1. The determinant is twice the area.
            det2 = ax * by - ay * bx
            vol[c] = abs(det2) / 2.0
            grad[c, 3] = (by) / det2
            grad[c, 4] = (-bx) / det2
            grad[c, 5] = (-ay) / det2
            grad[c, 6] = (ax) / det2
            grad[c, 1] = -grad[c, 3] - grad[c, 5]
            grad[c, 2] = -grad[c, 4] - grad[c, 6]
        else
            # The gradient of a basis function is a cross product of two edge
            # vectors divided by six times the signed volume. The two edges
            # leave vertex 1 and do not end in the vertex of the basis function.
            v3 = cells[c, 4]
            az = pts[v1, 3] - pts[v0, 3]
            bz = pts[v2, 3] - pts[v0, 3]
            cx = pts[v3, 1] - pts[v0, 1]
            cy = pts[v3, 2] - pts[v0, 2]
            cz = pts[v3, 3] - pts[v0, 3]
            det6 = ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) + az * (bx * cy - by * cx)
            vol[c] = abs(det6) / 6.0
            # grad of basis 2: (b x c) / det6
            grad[c, 4] = (by * cz - bz * cy) / det6
            grad[c, 5] = (bz * cx - bx * cz) / det6
            grad[c, 6] = (bx * cy - by * cx) / det6
            # grad of basis 3: (c x a) / det6
            grad[c, 7] = (cy * az - cz * ay) / det6
            grad[c, 8] = (cz * ax - cx * az) / det6
            grad[c, 9] = (cx * ay - cy * ax) / det6
            # grad of basis 4: (a x b) / det6
            grad[c, 10] = (ay * bz - az * by) / det6
            grad[c, 11] = (az * bx - ax * bz) / det6
            grad[c, 12] = (ax * by - ay * bx) / det6
            # The basis functions add up to 1, so their gradients add up to 0.
            grad[c, 1] = -grad[c, 4] - grad[c, 7] - grad[c, 10]
            grad[c, 2] = -grad[c, 5] - grad[c, 8] - grad[c, 11]
            grad[c, 3] = -grad[c, 6] - grad[c, 9] - grad[c, 12]
        end
    end
    return nothing
end

# p1_stiffness_kernel(grad, vol, cells, nc, nloc, dim, ii, jj, vv)
#
# Compute the entries of the stiffness matrix of the Laplace operator.
# Each cell gives nloc * nloc entries. Entries with the same row and
# column are added later, when the matrix is built.
#
# grad: gradients from p1_gradients_kernel
# vol: cell volumes
# cells: vertex numbers of each cell
# nc: number of cells
# nloc: number of vertices in a cell
# dim: 2 or 3
# ii, jj, vv: outputs, row, column and value of each entry
function p1_stiffness_kernel(grad, vol, cells, nc, nloc, dim, ii, jj, vv)
    # one entry per (cell, a, b), the entries are independent
    for idx = 1:nc * nloc * nloc
        c = div(idx - 1, nloc * nloc) + 1
        a = div(mod(idx - 1, nloc * nloc), nloc) + 1
        b = mod(idx - 1, nloc) + 1
        s = 0.0
        for d = 1:dim
            s += grad[c, (a - 1) * dim + d] * grad[c, (b - 1) * dim + d]
        end
        ii[idx] = cells[c, a]
        jj[idx] = cells[c, b]
        vv[idx] = s * vol[c]
    end
    return nothing
end

# prolongation_kernel(parents, nf, ii, jj, vv)
#
# Compute the entries of the prolongation matrix. A fine vertex that is a
# coarse vertex has the same parent twice. A fine vertex at the middle of
# a coarse edge has the two ends of the edge as parents. Each parent gets
# the weight 0.5, and the two weights add up to 1 when the parents are equal.
#
# parents: the array fine_vertex_parents, nf by 2
# nf: number of fine vertices
# ii, jj, vv: outputs, 2 entries for each fine vertex
function prolongation_kernel(parents, nf, ii, jj, vv)
    for idx = 1:2 * nf
        v = div(idx - 1, 2) + 1
        k = mod(idx - 1, 2) + 1
        ii[idx] = v
        jj[idx] = parents[v, k]
        vv[idx] = 0.5
    end
    return nothing
end

# gauss_seidel_kernel(colptr, rowval, nzval, n, b, x, forward)
#
# Do one Gauss-Seidel sweep for the system A x = b. A is symmetric, so the
# compressed columns of A are also its compressed rows.
#
# colptr, rowval, nzval: the compressed columns of A
# n: number of unknowns
# b: right-hand side
# x: current solution, changed in place
# forward: 1 to sweep from the first unknown to the last, 0 to sweep back
function gauss_seidel_kernel(colptr, rowval, nzval, n, b, x, forward)
    # Each unknown uses the new values of the unknowns before it, so the
    # order of the loop matters.
    for k = 1:n
        i = forward * k + (1 - forward) * (n + 1 - k)
        s = b[i]
        dg = 1.0
        for p = colptr[i]:colptr[i + 1] - 1
            j = rowval[p]
            if j == i
                dg = nzval[p]
            else
                s -= nzval[p] * x[j]
            end
        end
        x[i] = s / dg
    end
    return nothing
end

# ── Driver ───────────────────────────────────────────────────────────────

# Assemble the stiffness matrix of a mesh and the load vector for a source
# term of 1. Return the matrix as a sparse matrix.
function assemble(m)
    nc = size(m.cells, 1)
    nloc = m.dim + 1
    grad = zeros(Float64, nc, nloc * m.dim)
    vol = zeros(Float64, nc)
    p1_gradients_kernel(m.cells, m.pts, nc, m.dim, grad, vol)
    ii = zeros(Int, nc * nloc * nloc)
    jj = zeros(Int, nc * nloc * nloc)
    vv = zeros(Float64, nc * nloc * nloc)
    p1_stiffness_kernel(grad, vol, m.cells, nc, nloc, m.dim, ii, jj, vv)
    nv = size(m.pts, 1)
    A = sparse(ii, jj, vv, nv, nv)
    b = zeros(Float64, nv)
    for c in 1:nc, a in 1:nloc
        b[m.cells[c, a]] += vol[c] / nloc
    end
    return A, b, vol
end

# Build the prolongation matrix from a fine level to the level above it.
function prolongation(fine, nCoarse)
    nf = size(fine.pts, 1)
    ii = zeros(Int, 2 * nf)
    jj = zeros(Int, 2 * nf)
    vv = zeros(Float64, 2 * nf)
    prolongation_kernel(fine.hierarchy.vertexParents, nf, ii, jj, vv)
    return sparse(ii, jj, vv, nf, nCoarse)
end

# Apply one V-cycle on level k of the list lev. Each entry holds the
# matrix A of the free vertices, the matrix P to the next finer level and
# the factorization of the coarsest matrix. The number nu is the number of
# forward sweeps before the correction and of backward sweeps after it.
function vcycle!(lev, k, b, x, nu)
    A = lev[k].A
    if k == 1
        x .= lev[1].fact \ b
        return nothing
    end
    for _ in 1:nu
        gauss_seidel_kernel(A.colptr, A.rowval, A.nzval, size(A, 1), b, x, 1)
    end
    r = b - A * x
    P = lev[k].P
    ec = zeros(Float64, size(P, 2))
    vcycle!(lev, k - 1, P' * r, ec, nu)
    x .+= P * ec
    for _ in 1:nu
        gauss_seidel_kernel(A.colptr, A.rowval, A.nzval, size(A, 1), b, x, 0)
    end
    return nothing
end

# Group the vertices of one level into classes. Two vertices are in the same
# class when they are periodic images of each other. A vertex without a
# partner is a class of its own. The classes are numbered in the order of
# their lowest vertex, so vertex 1 is always in class 1. Return the class
# of each vertex and the number of classes.
function vertex_classes(m)
    n = size(m.pts, 1)
    root = collect(1:n)
    function find(x)
        while root[x] != x
            root[x] = root[root[x]]
            x = root[x]
        end
        return x
    end
    for r in 1:size(m.allPairs, 1)
        ra = find(m.allPairs[r, 1])
        rb = find(m.allPairs[r, 2])
        ra != rb && (root[max(ra, rb)] = min(ra, rb))
    end
    cls = zeros(Int, n)
    nCls = 0
    for v in 1:n
        if find(v) == v
            nCls += 1
            cls[v] = nCls
        else
            cls[v] = cls[find(v)]
        end
    end
    return cls, nCls
end

# Find the classes with a zero value (Dirichlet). A class has one when a
# vertex lies on a box side with no periodicity. It also has one when
# the vertex lies on an obstacle or on a step. These vertices have a
# marker above 2 times the dimension. Return a vector of Bool, one
# entry for each class. All entries are false when the mesh is periodic
# in every direction and has no obstacle or step.
function dirichlet_classes(m, cls, nCls)
    isD = fill(false, nCls)
    tol = 1e-9 * maximum(m.boxLength)
    for v in 1:size(m.pts, 1)
        hit = m.ptMk[v] > 2 * m.dim
        for d in 1:m.dim
            m.periodic[d] && continue
            hit = hit || abs(m.pts[v, d]) < tol || abs(m.pts[v, d] - m.boxLength[d]) < tol
        end
        hit && (isD[cls[v]] = true)
    end
    return isD
end

# Run the checks on one hierarchy. The argument yamlPath is the config
# of a hierarchy that writer.jl has already written. Return a vector of
# failed check names. A periodic direction is allowed. Vertices that are
# periodic images of each other are one unknown.
function check_hierarchy(name, yamlPath)
    levels = load_hierarchy(yamlPath)
    nL = length(levels)
    failed = String[]
    per = levels[1].periodic
    println("== $name: $(nL) levels, $(levels[1].dim)-D, periodic $(join(Int.(collect(per)), "")), cells ", join([size(l.cells, 1) for l in levels], " / "))
    cls = Vector{Int}[]
    nCls = Int[]
    isD = Vector{Bool}[]
    R = SparseMatrixCSC{Float64,Int}[]
    for l in levels
        c, n = vertex_classes(l)
        push!(cls, c)
        push!(nCls, n)
        push!(isD, dirichlet_classes(l, c, n))
        # R maps the value of each class to the value of each vertex.
        push!(R, sparse(1:length(c), c, 1.0, length(c), n))
    end
    # Sel maps the values at the vertices to the values of the classes. It takes the lowest vertex of each class.
    rep = [zeros(Int, nCls[k]) for k in 1:nL]
    for k in 1:nL, v in 1:length(cls[k])
        rep[k][cls[k][v]] == 0 && (rep[k][cls[k][v]] = v)
    end
    Sel = [sparse(1:nCls[k], rep[k], 1.0, nCls[k], size(levels[k].pts, 1)) for k in 1:nL]
    println("   classes of vertices: ", join(["$(nCls[k]) of $(size(levels[k].pts, 1))" for k in 1:nL], " / "))
    free = [findall(.!isD[k]) for k in 1:nL]
    Afull = [assemble(l)[1] for l in levels]
    bfull = [assemble(l)[2] for l in levels]
    Pfull = [k == 1 ? nothing : prolongation(levels[k], size(levels[k - 1].pts, 1)) for k in 1:nL]
    # Pc maps the values of the classes of level k-1 to the values of the classes of level k.
    Pc = [k == 1 ? nothing : Sel[k] * Pfull[k] * R[k - 1] for k in 1:nL]

    # Check 1: linear functions go from the coarse to the fine level unchanged.
    rng = MersenneTwister(7)
    for k in 2:nL
        g = randn(rng, 3)
        uc = levels[k - 1].pts * g .+ 0.3
        uf = levels[k].pts * g .+ 0.3
        err = maximum(abs.(Pfull[k] * uc - uf))
        rowsum = maximum(abs.(Pfull[k] * ones(size(Pfull[k], 2)) .- 1.0))
        println("   linear exactness  level $(k - 2) -> $(k - 1): error $(round(err, sigdigits = 3)), row sums $(round(rowsum, sigdigits = 3))")
        (err < 1e-12 && rowsum < 1e-14) || push!(failed, "$name: prolongation is not exact for linear functions at level $(k - 1)")
    end

    # Check 2: a periodic function of the coarse level stays periodic on the fine level.
    # This needs two periodic images on the fine level to have parents that are periodic images too.
    if any(per)
        for k in 2:nL
            uc = R[k - 1] * randn(rng, nCls[k - 1])
            uf = Pfull[k] * uc
            err = maximum(abs.(R[k] * (Sel[k] * uf) - uf))
            println("   periodic images   level $(k - 2) -> $(k - 1): largest jump between images $(round(err, sigdigits = 3))")
            err < 1e-13 || push!(failed, "$name: a periodic function is not periodic after the prolongation at level $(k - 1)")
        end
    end

    # Check 3: the Galerkin product equals the stiffness matrix of the coarse mesh.
    # The matrices are those of the unknowns of the classes.
    for k in 2:nL
        Af = R[k]' * Afull[k] * R[k]
        Ac = R[k - 1]' * Afull[k - 1] * R[k - 1]
        G = Pc[k]' * Af * Pc[k]
        err = norm(G - Ac, Inf) / norm(Ac, Inf)
        println("   Galerkin product  level $(k - 1) -> $(k - 2): relative difference $(round(err, sigdigits = 3))")
        err < 1e-12 || push!(failed, "$name: Galerkin product differs from the coarse matrix at level $(k - 2)")
    end

    # Check 4: a fine vertex on the boundary has no coarse parent that is free.
    # Without this, the zero values on the boundary would not nest.
    for k in 2:nL
        blk = Pc[k][isD[k], free[k - 1]]
        err = isempty(blk) ? 0.0 : maximum(abs.(blk))
        println("   boundary parents  level $(k - 2) -> $(k - 1): largest weight of a free coarse vertex on a fine boundary vertex $(round(err, sigdigits = 3))")
        err < 1e-14 || push!(failed, "$name: a fine boundary vertex has a free coarse parent at level $(k - 1)")
    end

    # Check 5: V-cycles. Level 0 stays the coarsest level. The finest level changes.
    for top in 2:nL
        lev = NamedTuple[]
        for k in 1:top
            I = free[k]
            A = R[k]' * Afull[k] * R[k]
            # a mesh with no zero value: add the lumped mass matrix, so the matrix is not singular
            any(isD[k]) || (A = A + R[k]' * spdiagm(0 => bfull[k]) * R[k])
            A = A[I, I]
            P = k == 1 ? nothing : Pc[k][free[k], free[k - 1]]
            push!(lev, (A = A, P = P, fact = k == 1 ? cholesky(Symmetric(A)) : nothing))
        end
        A = lev[top].A
        b = (R[top]' * bfull[top])[free[top]]
        exact = cholesky(Symmetric(A)) \ b
        # V(2,2)-cycle: 2 forward sweeps, correction, 2 backward sweeps
        x = zeros(Float64, length(b))
        res = [norm(b - A * x)]
        for _ in 1:12
            vcycle!(lev, top, b, x, 2)
            push!(res, norm(b - A * x))
        end
        rate = (res[end] / res[3])^(1 / (length(res) - 3))
        err = norm(x - exact) / norm(exact)
        println("   V(2,2)-cycle, $(top) levels, $(length(b)) unknowns: factor $(round(rate, digits = 3)), error after 12 cycles $(round(err, sigdigits = 3))")
        rate < 0.7 || push!(failed, "$name: V-cycle with $top levels is too slow (factor $(round(rate, digits = 3)))")
    end
    return failed
end

# Return the test hierarchies. Each entry holds a name and the config text.
# A function replaces a constant, so that the file has no global value.
function _configs()
    return [
    ("triangle_circle", """
mesh_type: triangle
box_length: [2.0, 1.0]
obstacles:
  - shape: circle
    center: [0.7, 0.5]
    radius: 0.2
mesh_options:
  mesh_size: 0.4
  mesh_size_obstacle: 0.2
hierarchy:
  levels: 5
"""),
    ("gridsplit_step", """
mesh_type: gridsplit
box_length: [4.0, 2.0]
n_divisions: [8, 4]
diagonal: slash
step:
  length: 1.0
  height: 1.0
hierarchy:
  levels: 4
"""),
    ("gridsplit_stretched", """
mesh_type: gridsplit
box_length: [1.0, 1.0]
n_divisions: [4, 3]
diagonal: slash
stretching_ratio: [1.3, 1.5]
stretching_symmetric: [false, true]
hierarchy:
  levels: 5
"""),
    ("kuhn_box", """
mesh_type: kuhn
box_length: [1.0, 1.0, 1.0]
n_divisions: [2, 2, 2]
diagonal: main
hierarchy:
  levels: 4
"""),
    ("triangle_periodic", """
mesh_type: triangle
box_length: [3.0, 1.0]
periodic: [true, true]
obstacles:
  - shape: circle
    center: [1.0, 0.5]
    radius: 0.2
  - shape: rectangle
    center: [2.0, 0.5]
    half_extents: [0.25, 0.15]
mesh_options:
  mesh_size: 0.4
  mesh_size_obstacle: 0.2
hierarchy:
  levels: 4
"""),
    ("gridsplit_periodic_x", """
mesh_type: gridsplit
box_length: [1.0, 1.0]
n_divisions: [4, 3]
diagonal: slash
periodic: [true, false]
hierarchy:
  levels: 5
"""),
    ("kuhn_periodic", """
mesh_type: kuhn
box_length: [1.0, 1.0, 1.0]
n_divisions: [2, 2, 2]
periodic: [true, false, true]
diagonal: main
hierarchy:
  levels: 4
"""),
    ("kuhn_fully_periodic", """
mesh_type: kuhn
box_length: [1.0, 1.0, 1.0]
n_divisions: [2, 2, 2]
periodic: [true, true, true]
diagonal: main
hierarchy:
  levels: 4
"""),
    ("tetra_periodic", """
mesh_type: tetra
box_length: [2.0, 1.0, 1.0]
periodic: [false, true, true]
obstacles:
  - shape: sphere
    center: [0.5, 0.5, 0.5]
    radius: 0.2
  - shape: box
    center: [1.4, 0.5, 0.5]
    half_extents: [0.15, 0.2, 0.15]
mesh_options:
  mesh_size: 0.5
  mesh_size_obstacle: 0.3
  max_dihedral: 150
hierarchy:
  levels: 3
"""),
    ("tetra_sphere", """
mesh_type: tetra
box_length: [1.0, 1.0, 1.0]
obstacles:
  - shape: sphere
    center: [0.5, 0.5, 0.5]
    radius: 0.2
mesh_options:
  mesh_size: 0.5
  mesh_size_obstacle: 0.25
  max_dihedral: 150
hierarchy:
  levels: 4
"""),
    ]
end

function run_all()
    dir = mktempdir()
    failed = String[]
    configs = _configs()
    names = isempty(ARGS) ? [c[1] for c in configs] : ARGS
    for (name, text) in configs
        name in names || continue
        yamlPath = joinpath(dir, name * ".yaml")
        write(yamlPath, text)
        redirect_stdout(devnull) do
            main(yamlPath)
        end
        append!(failed, check_hierarchy(name, yamlPath))
    end
    println()
    if isempty(failed)
        println("All checks passed.")
    else
        println("Failed checks:")
        foreach(f -> println("  ", f), failed)
        exit(1)
    end
    return nothing
end

if abspath(PROGRAM_FILE) == (@__FILE__)
    run_all()
end
