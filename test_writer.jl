#!/usr/bin/env julia
#=
Test script for writer.jl (the merged writer for all five mesh types).

Usage:
    julia test_writer.jl [path/to/writer.jl]

The script reads the writer file. It cuts out the package-free part and
evaluates that part in the module WCore. The tests of the predicates, the
mesh invariants, the helper functions, the config check, and the names
run in WCore. They need no YAML or HDF5 package. The end-to-end tests
load the whole file in a second module. They run only when the YAML and
HDF5 packages are present.
The environment variable N_RANDOM sets the number of random cases. The
default is 25 for 2-D and 12 for 3-D.
The environment variable WRITER_REF_DIR can name a folder that holds
reference outputs (see tools/README.md). The regression test then
compares the new files with them, bit by bit.
=#

using Random
using LinearAlgebra

writerPath = abspath(length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "writer.jl"))
text = read(writerPath, String)
i0 = first(findfirst("# package-free part: begins", text))
i1 = last(findfirst("# package-free part: ends", text))
WCore = Module(:WCore)
Base.include_string(WCore, "using LinearAlgebra\nusing Random\n" * text[i0:i1])

nFail = 0
nPass = 0
# Count one check. Print a line for each failed check.
function check(name, ok)
    global nFail, nPass
    if ok
        nPass += 1
    else
        nFail += 1
        println("  FAIL: ", name)
    end
end

# Return true when f raises an error.
function raises(f)
    try
        f()
    catch
        return true
    end
    return false
end

# Return the message of the error that f raises, or nothing.
function errmsg(f)
    try
        f()
    catch e
        return sprint(showerror, e)
    end
    return nothing
end

# ======================================================================
# 2-D mesher tests
# ======================================================================

# ── Predicate tests ──────────────────────────────────────────────────────

# Exact sign of the orientation determinant, with rational numbers.
function orient_oracle(ax, ay, bx, by, cx, cy)
    r(x) = Rational{BigInt}(x)
    d = (r(ax) - r(cx)) * (r(by) - r(cy)) - (r(ay) - r(cy)) * (r(bx) - r(cx))
    return Int(sign(d))
end

# Exact sign of the incircle determinant, with rational numbers.
function incircle_oracle(ax, ay, bx, by, cx, cy, dx, dy)
    r(x) = Rational{BigInt}(x)
    adx = r(ax) - r(dx); ady = r(ay) - r(dy)
    bdx = r(bx) - r(dx); bdy = r(by) - r(dy)
    cdx = r(cx) - r(dx); cdy = r(cy) - r(dy)
    al = adx^2 + ady^2; bl = bdx^2 + bdy^2; cl = cdx^2 + cdy^2
    d = al * (bdx * cdy - cdx * bdy) + bl * (cdx * ady - adx * cdy) + cl * (adx * bdy - bdx * ady)
    return Int(sign(d))
end

# Run the predicate tests.
function predicate_tests2()
    println("Predicate tests (2-D)")
    rng = MersenneTwister(7)
    bad = 0
    for _ in 1:20000
        p = rand(rng, 6)
        WCore._orient2d(p...) == orient_oracle(p...) || (bad += 1)
    end
    check("orient2d agrees with the oracle on random input", bad == 0)
    bad = 0
    for _ in 1:20000
        p = rand(rng, 8)
        WCore._incircle(p...) == incircle_oracle(p...) || (bad += 1)
    end
    check("incircle agrees with the oracle on random input", bad == 0)
    grid = [0.1 * i for i in 0:6]
    bad = 0
    nCollinear = 0
    for _ in 1:20000
        q = [(grid[rand(rng, 1:7)], grid[rand(rng, 1:7)]) for _ in 1:4]
        o = orient_oracle(q[1]..., q[2]..., q[3]...)
        o == 0 && (nCollinear += 1)
        WCore._orient2d(q[1]..., q[2]..., q[3]...) == o || (bad += 1)
        WCore._incircle(q[1]..., q[2]..., q[3]..., q[4]...) == incircle_oracle(q[1]..., q[2]..., q[3]..., q[4]...) || (bad += 1)
    end
    check("both predicates agree with the oracle on a grid", bad == 0 && nCollinear > 1000)
end
predicate_tests2()

# ── Mesh tests ───────────────────────────────────────────────────────────

# Build one mesh with the core function. Return the vertices, cells, edges,
# vertex markers, edge markers, and info (the cell-to-edge array is left out).
function build2(Lx, Ly, perX, perY, obs, h, hObs, grading; sweeps = 3, minAngle = 20.0)
    n = length(obs)
    kinds = Int[o[1] for o in obs]
    cx = Float64[o[2] for o in obs]
    cy = Float64[o[3] for o in obs]
    hx = Float64[o[4] for o in obs]
    hy = Float64[o[5] for o in obs]
    r = WCore._triangulate_domain2(IOBuffer(), Lx, Ly, perX, perY, kinds, cx, cy, hx, hy,
                                   h, hObs, grading, minAngle, sweeps, 2000000)
    return (r[1], r[2], r[3], r[5], r[6], r[7])
end

# Signed area of the domain from its boundary edges only (Green formula).
function boundary_area(pts, cells)
    area = 0.0
    count = Dict{Tuple{Int,Int},Int}()
    for c in 1:size(cells, 1), (a, b) in ((1, 2), (2, 3), (3, 1))
        u, w = cells[c, a], cells[c, b]
        key = u < w ? (u, w) : (w, u)
        count[key] = get(count, key, 0) + 1
    end
    for c in 1:size(cells, 1), (a, b) in ((1, 2), (2, 3), (3, 1))
        u, w = cells[c, a], cells[c, b]
        key = u < w ? (u, w) : (w, u)
        if count[key] == 1
            area += 0.5 * (pts[u, 1] * pts[w, 2] - pts[w, 1] * pts[u, 2])
        end
    end
    return area
end

# Run all mesh checks on one mesh.
function check_mesh2(name, Lx, Ly, perX, perY, obs, mesh; exactArea = nothing, delaunay = false)
    pts, cells, edges, ptMk, edgeMk, info = mesh
    nCells = size(cells, 1)
    nPts = size(pts, 1)
    nEdges = size(edges, 1)
    k = length(obs)
    # 1. positive areas
    meas = zeros(nCells)
    WCore.signed_measure_2d_kernel(cells, pts, nCells, meas)
    check("$name: every cell has a positive area", all(>(0), meas))
    # 2. area from the cells and from the boundary edges
    total = sum(meas)
    check("$name: cell area equals boundary area", abs(total - boundary_area(pts, cells)) < 1e-11 * Lx * Ly)
    if exactArea !== nothing
        check("$name: area equals the exact area", abs(total - exactArea) < 1e-11 * Lx * Ly)
    end
    # 3. Euler characteristic
    check("$name: Euler characteristic is 1 - k", nPts - nEdges + nCells == 1 - k)
    # 4. edge counts
    _, c2e = WCore._build_edge_registry(cells, 2)
    cnt = zeros(Int, nEdges)
    for c in 1:nCells, l in 1:3
        cnt[c2e[c, l] + 1] += 1
    end
    check("$name: each edge belongs to 1 or 2 cells", all(x -> x == 1 || x == 2, cnt))
    # 5. locally Delaunay edges, tested only when the mesh is not smoothed
    if delaunay
        ok = true
        owner = Dict{Tuple{Int,Int},Tuple{Int,Int}}()
        for c in 1:nCells, l in 1:3
            u, w = cells[c, l % 3 + 1], cells[c, (l + 1) % 3 + 1]
            key = u < w ? (u, w) : (w, u)
            if haskey(owner, key)
                c0, l0 = owner[key]
                apex = cells[c, l]
                a, b, d = cells[c0, 1], cells[c0, 2], cells[c0, 3]
                if WCore._incircle(pts[a, 1], pts[a, 2], pts[b, 1], pts[b, 2], pts[d, 1], pts[d, 2], pts[apex, 1], pts[apex, 2]) > 0
                    ok = false
                end
            else
                owner[key] = (c, l)
            end
        end
        check("$name: every interior edge is locally Delaunay", ok)
    end
    # 6. smallest angle
    check("$name: smallest angle is at least 20 degrees", info.minAngleDeg >= 20.0 - 1e-9)
    # 7. markers
    onBoundaryVertex = falses(nPts)
    for e in 1:nEdges
        if cnt[e] == 1
            onBoundaryVertex[edges[e, 1]] = true
            onBoundaryVertex[edges[e, 2]] = true
        end
    end
    check("$name: boundary vertices have a marker above 0", all(ptMk[v] > 0 for v in 1:nPts if onBoundaryVertex[v]))
    check("$name: interior vertices have marker 0", all(ptMk[v] == 0 for v in 1:nPts if !onBoundaryVertex[v]))
    check("$name: boundary edges have a marker above 0", all(edgeMk[e] > 0 for e in 1:nEdges if cnt[e] == 1))
    check("$name: interior edges have marker 0", all(edgeMk[e] == 0 for e in 1:nEdges if cnt[e] == 2))
    markersSeen = Set(Int.(ptMk))
    for i in 0:k - 1
        check("$name: marker $(10 + i) is present", (10 + i) in markersSeen)
    end
    # 8. vertex normals sum to zero
    normal = zeros(nCells, 3, 3)
    WCore.cell_vertex_normal_2d_kernel(cells, pts, nCells, normal)
    check("$name: the vertex normals of each cell sum to zero", maximum(abs.(sum(normal, dims = 2))) < 1e-12)
    # 9. periodic pairs
    for (d, per, L) in ((1, perX, Lx), (2, perY, Ly))
        per || continue
        o = 3 - d
        lo = sort([pts[v, o] for v in 1:nPts if abs(pts[v, d]) < 1e-9])
        hi = sort([pts[v, o] for v in 1:nPts if abs(pts[v, d] - L) < 1e-9])
        check("$name: periodic axis $d has matching sides", length(lo) == length(hi) && all(abs.(lo .- hi) .< 1e-9))
    end
end

println("Mesh tests (2-D)")
configs2 = [
    ("box", 1.0, 1.0, false, false, [], 0.1, 0.1, 0.3, 1.0),
    ("box periodic x", 1.0, 2.0, true, false, [], 0.08, 0.08, 0.3, 2.0),
    ("circle", 2.0, 1.0, false, false, [(1, 1.0, 0.5, 0.2, 0.2)], 0.1, 0.04, 0.3, nothing),
    ("rectangle", 2.0, 1.0, false, false, [(2, 1.0, 0.5, 0.3, 0.15)], 0.1, 0.04, 0.3, 2.0 - 0.18),
    ("mixed periodic xy", 3.0, 1.0, true, true,
        [(1, 0.8, 0.5, 0.15, 0.15), (2, 1.6, 0.5, 0.25, 0.1), (1, 2.4, 0.5, 0.12, 0.12)], 0.08, 0.03, 0.3, nothing),
    ("obstacle near a wall", 2.0, 1.0, false, false, [(1, 1.0, 0.3, 0.2, 0.2)], 0.1, 0.05, 0.3, nothing),
    ("fine obstacle grading 0", 1.0, 1.0, false, false, [(2, 0.5, 0.5, 0.1, 0.1)], 0.1, 0.02, 0.0, 1.0 - 0.04),
]
for (name, Lx, Ly, px, py, obs, h, hObs, g, exact) in configs2
    mesh = build2(Lx, Ly, px, py, obs, h, hObs, g)
    check_mesh2(name, Lx, Ly, px, py, obs, mesh; exactArea = exact)
    meshA = build2(Lx, Ly, px, py, obs, h, hObs, g)
    check("$name: two runs give the same mesh", meshA[1] == mesh[1] && meshA[2] == mesh[2])
    meshD = build2(Lx, Ly, px, py, obs, h, hObs, g; sweeps = 0)
    check_mesh2(name * " (no smoothing)", Lx, Ly, px, py, obs, meshD; delaunay = true)
end

# The circle area must stay below the disk area and near it.
meshC = build2(2.0, 1.0, false, false, [(1, 1.0, 0.5, 0.2, 0.2)], 0.1, 0.02, 0.3)
areaC = 2.0 - sum(let m = zeros(size(meshC[2], 1)); WCore.signed_measure_2d_kernel(meshC[2], meshC[1], size(meshC[2], 1), m); m end)
check("circle: the polygon area is close to the disk area", 0.0 < pi * 0.2^2 - areaC < 2e-3)
check("circle: every circle vertex lies on the circle",
      all(abs(hypot(meshC[1][v, 1] - 1.0, meshC[1][v, 2] - 0.5) - 0.2) < 1e-12 for v in 1:size(meshC[1], 1) if meshC[4][v] == 10))

# ── Random configurations ────────────────────────────────────────────────

# Build random boxes with random obstacles. Check every mesh. The
# environment variable N_RANDOM sets the number of cases. The default is 25.
function random_tests2(nCases)
    rng = MersenneTwister(2024)
    nOk = 0
    for t in 1:nCases
        Lx = 1.0 + 3.0 * rand(rng)
        Ly = 0.5 + 1.5 * rand(rng)
        h = 0.04 + 0.12 * rand(rng)
        hObs = h * (0.2 + 0.8 * rand(rng))
        g = rand(rng) < 0.2 ? 0.0 : 0.1 + 0.6 * rand(rng)
        perX = rand(rng) < 0.4
        perY = rand(rng) < 0.4
        obs = Tuple[]
        for _ in 1:rand(rng, 0:4)
            kind = rand(rng, 1:2)
            hx = 0.04 + 0.2 * rand(rng)
            hy = kind == 1 ? hx : 0.04 + 0.2 * rand(rng)
            cx = hx + 0.02 + (Lx - 2 * hx - 0.04) * rand(rng)
            cy = hy + 0.02 + (Ly - 2 * hy - 0.04) * rand(rng)
            if cx - hx <= 0.01 || cx + hx >= Lx - 0.01 || cy - hy <= 0.01 || cy + hy >= Ly - 0.01
                continue
            end
            ok = true
            for o in obs
                sep = max(o[2] - o[4] - (cx + hx), cx - hx - (o[2] + o[4]), o[3] - o[5] - (cy + hy), cy - hy - (o[3] + o[5]))
                sep <= 0.02 && (ok = false)
            end
            ok && push!(obs, (kind, cx, cy, hx, hy))
        end
        name = "random $t"
        local mesh
        try
            mesh = build2(Lx, Ly, perX, perY, obs, h, hObs, g)
        catch e
            check("$name: the mesh builds ($(sprint(showerror, e)))", false)
            continue
        end
        before = nFail
        check_mesh2(name, Lx, Ly, perX, perY, obs, mesh)
        nFail == before && (nOk += 1)
    end
    println("  random cases without a failed check: $nOk of $nCases")
end
println("Random configurations (2-D)")
random_tests2(parse(Int, get(ENV, "N_RANDOM", "25")))

# ── Size field ───────────────────────────────────────────────────────────

println("Option and obstacle checks (2-D)")
check("mesh_size is required", raises(() -> WCore._resolve_mesh_options2(Dict("grading" => 0.1))))
check("unknown mesh_options key is an error", raises(() -> WCore._resolve_mesh_options2(Dict("mesh_size" => 0.1, "foo" => 1))))
check("min_angle above 20.7 is an error", raises(() -> WCore._resolve_mesh_options2(Dict("mesh_size" => 0.1, "min_angle" => 25))))
box2 = (2.0, 1.0)
check("obstacle that crosses the box is an error",
      raises(() -> WCore._validate_obstacles2(IOBuffer(), [Dict("shape" => "circle", "center" => [0.1, 0.5], "radius" => 0.2)], box2, 0.1)))
check("overlapping bounding boxes are an error",
      raises(() -> WCore._validate_obstacles2(IOBuffer(), [Dict("shape" => "circle", "center" => [0.8, 0.5], "radius" => 0.2),
                                             Dict("shape" => "circle", "center" => [1.1, 0.5], "radius" => 0.2)], box2, 0.1)))
check("a valid obstacle passes",
      !raises(() -> WCore._validate_obstacles2(IOBuffer(), [Dict("shape" => "circle", "center" => [1.0, 0.5], "radius" => 0.2)], box2, 0.1)))

# ======================================================================
# 3-D mesher tests
# ======================================================================

# ── Predicate tests ──────────────────────────────────────────────────────

# Exact sign of the orientation determinant of (b - a, c - a, d - a).
function orient3_oracle(a, b, c, d)
    r(x) = Rational{BigInt}(x)
    m = [r(b[1]) - r(a[1]) r(b[2]) - r(a[2]) r(b[3]) - r(a[3]);
         r(c[1]) - r(a[1]) r(c[2]) - r(a[2]) r(c[3]) - r(a[3]);
         r(d[1]) - r(a[1]) r(d[2]) - r(a[2]) r(d[3]) - r(a[3])]
    dt = m[1, 1] * (m[2, 2] * m[3, 3] - m[2, 3] * m[3, 2]) - m[1, 2] * (m[2, 1] * m[3, 3] - m[2, 3] * m[3, 1]) +
         m[1, 3] * (m[2, 1] * m[3, 2] - m[2, 2] * m[3, 1])
    return Int(sign(dt))
end

# Exact sign of the insphere test: 1 inside, -1 outside, 0 on the sphere.
# The cell (a, b, c, d) must have a positive orientation. The oracle
# builds the circumcenter with rational numbers.
function insphere_oracle(a, b, c, d, e)
    r(x) = Rational{BigInt}(x)
    P = [map(r, a), map(r, b), map(r, c), map(r, d)]
    A = [2 * (P[i][k] - P[1][k]) for i in 2:4, k in 1:3]
    rhs = [sum(P[i][k]^2 for k in 1:3) - sum(P[1][k]^2 for k in 1:3) for i in 2:4]
    det3(M) = M[1, 1] * (M[2, 2] * M[3, 3] - M[2, 3] * M[3, 2]) - M[1, 2] * (M[2, 1] * M[3, 3] - M[2, 3] * M[3, 1]) +
              M[1, 3] * (M[2, 1] * M[3, 2] - M[2, 2] * M[3, 1])
    D = det3(A)
    ctr = Vector{Rational{BigInt}}(undef, 3)
    for k in 1:3
        Ak = copy(A)
        Ak[:, k] = rhs
        ctr[k] = det3(Ak) / D
    end
    R2 = sum((P[1][k] - ctr[k])^2 for k in 1:3)
    E = map(r, e)
    return Int(sign(R2 - sum((E[k] - ctr[k])^2 for k in 1:3)))
end

# Run the predicate tests.
function predicate_tests3()
    println("Predicate tests (3-D)")
    rng = MersenneTwister(11)
    bad = 0
    badI = 0
    for _ in 1:4000
        p = [rand(rng, 3) for _ in 1:5]
        o = WCore._orient3d(p[1]..., p[2]..., p[3]..., p[4]...)
        o == orient3_oracle(p[1], p[2], p[3], p[4]) || (bad += 1)
        if o < 0
            p[3], p[4] = p[4], p[3]
        end
        WCore._insphere(p[1]..., p[2]..., p[3]..., p[4]..., p[5]...) == insphere_oracle(p[1], p[2], p[3], p[4], p[5]) || (badI += 1)
    end
    check("orient3d agrees with the oracle on random input", bad == 0)
    check("insphere agrees with the oracle on random input", badI == 0)
    grid = [0.1 * i for i in 0:4]
    bad = 0
    badI = 0
    nCoplanar = 0
    nCospherical = 0
    for _ in 1:4000
        p = [[grid[rand(rng, 1:5)], grid[rand(rng, 1:5)], grid[rand(rng, 1:5)]] for _ in 1:5]
        o = orient3_oracle(p[1], p[2], p[3], p[4])
        WCore._orient3d(p[1]..., p[2]..., p[3]..., p[4]...) == o || (bad += 1)
        o == 0 && (nCoplanar += 1)
        if o != 0
            if o < 0
                p[3], p[4] = p[4], p[3]
            end
            io = insphere_oracle(p[1], p[2], p[3], p[4], p[5])
            WCore._insphere(p[1]..., p[2]..., p[3]..., p[4]..., p[5]...) == io || (badI += 1)
            io == 0 && (nCospherical += 1)
        end
    end
    check("orient3d agrees with the oracle on a lattice", bad == 0 && nCoplanar > 200)
    check("insphere agrees with the oracle on a lattice", badI == 0 && nCospherical > 50)
    # The 2-D predicates that the facet layer uses.
    bad = 0
    for _ in 1:4000
        p = rand(rng, 8)
        r(x) = Rational{BigInt}(x)
        adx = r(p[1]) - r(p[7]); ady = r(p[2]) - r(p[8]); bdx = r(p[3]) - r(p[7]); bdy = r(p[4]) - r(p[8])
        cdx = r(p[5]) - r(p[7]); cdy = r(p[6]) - r(p[8])
        det = (adx^2 + ady^2) * (bdx * cdy - cdx * bdy) + (bdx^2 + bdy^2) * (cdx * ady - adx * cdy) +
              (cdx^2 + cdy^2) * (adx * bdy - bdx * ady)
        WCore._incircle(p...) == Int(sign(det)) || (bad += 1)
    end
    check("incircle agrees with the oracle on random input", bad == 0)
end
predicate_tests3()

# ── Mesh tests ───────────────────────────────────────────────────────────

# Build one mesh with the core function. The obstacle tuples are
# (kind, cx, cy, cz, hx, hy, hz). Kind 1 is a sphere. Kind 2 is a box.
function build3(L, per, obs, h, hObs, grading; sweeps = 3, bound = 2.0, minDihedral = 0.0)
    kinds = Int[o[1] for o in obs]
    f(i) = Float64[o[i] for o in obs]
    return WCore._triangulate_domain3(IOBuffer(), L, per, kinds, f(2), f(3), f(4), f(5), f(6), f(7),
                                      h, hObs, grading, bound, sweeps, 2000000, minDihedral)
end

# Return the faces of a cell list: a Dict from a sorted vertex triple to
# the list of (cell, local face) pairs.
function face_table(cells)
    table = Dict{NTuple{3,Int},Vector{Tuple{Int,Int}}}()
    for c in 1:size(cells, 1), i in 1:4
        f = WCore._face3(i)
        t = sort([cells[c, f[1]], cells[c, f[2]], cells[c, f[3]]])
        push!(get!(table, (t[1], t[2], t[3]), Tuple{Int,Int}[]), (c, i))
    end
    return table
end

# Run all mesh checks on one mesh.
function check_mesh3(name, L, per, obs, mesh; exactVolume = nothing, delaunay = false, bound = 2.0, hObs = 0.0)
    pts, cells, edges, c2e, ptMk, info = mesh
    nCells = size(cells, 1)
    nPts = size(pts, 1)
    nEdges = size(edges, 1)
    k = length(obs)
    Lq = info.boxSize
    # 1. positive volumes
    meas = zeros(nCells)
    WCore.signed_measure_3d_kernel(cells, pts, nCells, meas)
    check("$name: every cell has a positive volume", all(>(0), meas))
    check("$name: no cell is nearly flat", minimum(meas) > 1e-12 * prod(Lq) / nCells)
    # 2. volume from the boundary faces (divergence theorem)
    table = face_table(cells)
    bvol = 0.0
    for (key, owners) in table
        length(owners) == 1 || continue
        c, i = owners[1]
        f = WCore._face3(i)
        a = pts[cells[c, f[1]], :]; b = pts[cells[c, f[2]], :]; d = pts[cells[c, f[3]], :]
        bvol += dot(a, cross(b, d)) / 6.0
    end
    total = sum(meas)
    check("$name: cell volume equals boundary volume", abs(total - bvol) < 1e-10 * prod(Lq))
    if exactVolume !== nothing
        check("$name: volume equals the exact volume", abs(total - exactVolume) < 1e-8 * prod(L))
    end
    # 3. Euler characteristic
    nFaces = length(table)
    check("$name: Euler characteristic is 1 + k", nPts - nEdges + nFaces - nCells == 1 + k)
    # 4. face counts and a closed boundary
    check("$name: each face belongs to 1 or 2 cells", all(v -> length(v) in (1, 2), values(table)))
    edgeCount = Dict{Tuple{Int,Int},Int}()
    for (key, owners) in table
        length(owners) == 1 || continue
        for (u, w) in ((key[1], key[2]), (key[1], key[3]), (key[2], key[3]))
            edgeCount[(u, w)] = get(edgeCount, (u, w), 0) + 1
        end
    end
    check("$name: the boundary surface is closed", all(==(2), values(edgeCount)))
    # 5. locally Delaunay faces, tested only when the mesh is not smoothed
    if delaunay
        ok = true
        for (key, owners) in table
            length(owners) == 2 || continue
            (c1, i1), (c2, i2) = owners
            a, b, d, e = cells[c1, 1], cells[c1, 2], cells[c1, 3], cells[c1, 4]
            apex = cells[c2, i2]
            if WCore._insphere(pts[a, :]..., pts[b, :]..., pts[d, :]..., pts[e, :]..., pts[apex, :]...) > 0
                ok = false
            end
        end
        check("$name: every interior face is locally Delaunay", ok)
        check("$name: the radius-edge ratio is within the bound", info.maxRatio <= bound * (1 + 1e-9))
        # 6. no boundary face is encroached on
        bad = 0
        for (key, owners) in table
            length(owners) == 1 || continue
            a = pts[key[1], :]; b = pts[key[2], :]; d = pts[key[3], :]
            u = b - a; v = d - a; w = cross(u, v)
            o = a + (dot(v, v) * cross(w, u) + dot(u, u) * cross(v, w)) / (2 * dot(w, w))
            R2 = sum((a - o) .^ 2)
            for q in 1:nPts
                (q == key[1] || q == key[2] || q == key[3]) && continue
                sum((pts[q, :] - o) .^ 2) < R2 * (1 - 1e-6) && (bad += 1; break)
            end
        end
        check("$name: no boundary face is encroached on", bad == 0)
    end
    # 7. markers
    onBoundary = falses(nPts)
    for (key, owners) in table
        length(owners) == 1 || continue
        for v in key
            onBoundary[v] = true
        end
    end
    check("$name: boundary vertices have a marker above 0", all(ptMk[v] > 0 for v in 1:nPts if onBoundary[v]))
    check("$name: interior vertices have marker 0", all(ptMk[v] == 0 for v in 1:nPts if !onBoundary[v]))
    tol = 1e-8 * maximum(Lq)
    okBox = true
    for v in 1:nPts
        expected = typemax(Int)
        for a in 1:3
            abs(pts[v, a]) < tol && (expected = min(expected, 2 * (a - 1) + 1))
            abs(pts[v, a] - Lq[a]) < tol && (expected = min(expected, 2 * (a - 1) + 2))
        end
        if expected != typemax(Int)
            ptMk[v] == expected || (okBox = false)
        else
            (ptMk[v] == 0 || ptMk[v] >= 100) || (okBox = false)
        end
    end
    check("$name: vertices on box faces have the box markers", okBox)
    markersSeen = Set(Int.(ptMk))
    for i in 0:k - 1
        check("$name: marker $(100 + i) is present", (100 + i) in markersSeen)
    end
    for (i, o) in enumerate(obs)
        mk = 100 + i - 1
        ids = [v for v in 1:nPts if ptMk[v] == mk]
        if o[1] == 1
            r = o[5]
            d = [norm(pts[v, :] - [o[2], o[3], o[4]]) for v in ids]
            check("$name: sphere $(i - 1) vertices lie on or just inside the sphere",
                  all(d .<= r * (1 + 1e-9)) && all(d .>= r - 2 * (max(hObs, 1e-3))^2 / r))
        else
            lo = [o[2] - o[5], o[3] - o[6], o[4] - o[7]]
            hi = [o[2] + o[5], o[3] + o[6], o[4] + o[7]]
            ok = all(v -> all(pts[v, :] .>= lo .- 1e-8) && all(pts[v, :] .<= hi .+ 1e-8) &&
                          any(abs.(pts[v, :] .- lo) .< 1e-8) || any(abs.(pts[v, :] .- hi) .< 1e-8), ids)
            check("$name: box obstacle $(i - 1) vertices lie on its surface", ok)
        end
    end
    # 8. vertex normals sum to zero
    normal = zeros(nCells, 4, 3)
    WCore.cell_vertex_normal_3d_kernel(cells, pts, nCells, normal)
    check("$name: the vertex normals of each cell sum to zero", maximum(abs.(sum(normal, dims = 2))) < 1e-12)
    # 9. periodic pairs
    for d in 1:3
        per[d] || continue
        o = [a for a in 1:3 if a != d]
        lo = sort([(pts[v, o[1]], pts[v, o[2]]) for v in 1:nPts if abs(pts[v, d]) < tol])
        hi = sort([(pts[v, o[1]], pts[v, o[2]]) for v in 1:nPts if abs(pts[v, d] - Lq[d]) < tol])
        ok = length(lo) == length(hi) && all(abs(lo[j][1] - hi[j][1]) < 1e-9 && abs(lo[j][2] - hi[j][2]) < 1e-9 for j in eachindex(lo))
        check("$name: periodic axis $d has matching faces", ok)
    end
    # 10. the registry arrays
    check("$name: cell-to-edge agrees with the edges", all(edges[c2e[c, li] + 1, :] == sort([cells[c, a], cells[c, b]])
                                                              for c in 1:min(nCells, 300) for (li, (a, b)) in enumerate(WCore._tet_local_edges())))
end

println("Mesh tests (3-D)")
configs3 = [
    ("box", (1.0, 1.0, 1.0), (false, false, false), [], 0.25, 0.25, 0.3, 1.0),
    ("box periodic x", (1.0, 1.5, 1.0), (true, false, false), [], 0.25, 0.25, 0.3, 1.5),
    ("box periodic xyz", (1.0, 1.0, 1.0), (true, true, true), [], 0.25, 0.25, 0.3, 1.0),
    ("box obstacle", (1.0, 1.0, 1.0), (false, false, false), [(2, 0.5, 0.5, 0.5, 0.15, 0.15, 0.15)], 0.2, 0.1, 0.3, 1.0 - 0.027),
    ("sphere", (1.0, 1.0, 1.0), (false, false, false), [(1, 0.5, 0.5, 0.5, 0.2, 0.2, 0.2)], 0.2, 0.1, 0.3, nothing),
    ("sphere grading 0", (1.0, 1.0, 1.0), (false, false, false), [(1, 0.5, 0.5, 0.5, 0.25, 0.25, 0.25)], 0.2, 0.1, 0.0, nothing),
    ("two obstacles periodic", (2.0, 1.0, 1.0), (false, true, true),
        [(1, 0.7, 0.5, 0.5, 0.2, 0.2, 0.2), (2, 1.4, 0.5, 0.5, 0.15, 0.2, 0.15)], 0.2, 0.1, 0.3, nothing),
    ("sphere near a wall", (1.0, 1.0, 1.0), (false, false, false), [(1, 0.5, 0.5, 0.3, 0.2, 0.2, 0.2)], 0.2, 0.1, 0.3, nothing),
]
for (name, L, per, obs, h, hObs, g, exact) in configs3
    mesh = build3(L, per, obs, h, hObs, g)
    check_mesh3(name, L, per, obs, mesh; exactVolume = exact, hObs = hObs)
    meshA = build3(L, per, obs, h, hObs, g)
    check("$name: two runs give the same mesh", meshA[1] == mesh[1] && meshA[2] == mesh[2])
    meshD = build3(L, per, obs, h, hObs, g; sweeps = 0)
    check_mesh3(name * " (no smoothing)", L, per, obs, meshD; delaunay = true, hObs = hObs)
end

# ── Flip pass ────────────────────────────────────────────────────────────

# Return the faces that belong to exactly one cell, as sorted vertex triples.
function boundary_faces(cells)
    return sort!([k for (k, v) in face_table(cells) if length(v) == 1])
end

# Return the total volume of a cell list, and whether every cell has a positive volume.
function volume_info(pts, cells)
    meas = zeros(size(cells, 1))
    WCore.signed_measure_3d_kernel(cells, pts, size(cells, 1), meas)
    return sum(meas), all(>(0), meas)
end

println("Flip pass")
# Unit tests on two tetrahedra with a common face.
let
    # The segment between the two apexes crosses the common triangle, so a 2-3 flip is valid.
    pts = [0.0 0.0 0.0; 1.0 0.0 0.0; 0.5 1.0 0.0; 0.5 0.3 0.5; 0.5 0.3 -0.5]
    cells = [1 2 3 4; 1 3 2 5]
    vol0, pos0 = volume_info(pts, cells)
    check("flip unit: the two starting cells have a positive volume", pos0)
    nbr = WCore._cell_neighbors3(cells)
    check("flip unit: the two cells are neighbors across one face", count(!=(0), nbr[1, :]) == 1 && count(!=(0), nbr[2, :]) == 1 && 2 in nbr[1, :] && 1 in nbr[2, :])
    out, nTried, nDone = WCore._flip_pass3(pts, cells, deg2rad(89.0), 2.0, -Inf)
    vol1, pos1 = volume_info(pts, out)
    check("flip unit: a valid 2-3 flip turns 2 cells into 3", nDone == 1 && size(out, 1) == 3)
    check("flip unit: the 2-3 flip keeps the volume and the positive orientation", isapprox(vol1, vol0; rtol = 1e-12) && pos1)
    check("flip unit: the 2-3 flip uses the same vertices", sort(unique(vec(out))) == [1, 2, 3, 4, 5])
    check("flip unit: the new edge joins the two apexes in all 3 cells", all(4 in out[c, :] && 5 in out[c, :] for c in 1:3))
    # The 3 cells now form a closed ring around the edge 4-5.
    nbr3 = WCore._cell_neighbors3(out)
    cA = 1
    iA = findfirst(==(4), out[cA, :]); jA = findfirst(==(5), out[cA, :])
    check("flip unit: the edge between the apexes has a closed ring of 3 cells", length(WCore._edge_ring3(out, nbr3, cA, iA, jA)) == 3)
    other = [(i, j) for (i, j) in WCore._tet_local_edges() if (out[cA, i], out[cA, j]) != (4, 5) && (out[cA, i], out[cA, j]) != (5, 4)]
    check("flip unit: no other edge of that cell has a ring of 3", all(isempty(WCore._edge_ring3(out, nbr3, cA, i, j)) for (i, j) in other))
    # Flip back. The result has 2 or 3 cells, with the same volume.
    out2, _, _ = WCore._flip_pass3(pts, out, deg2rad(89.0), 2.0, -Inf)
    vol2, pos2 = volume_info(pts, out2)
    check("flip unit: a second pass keeps the volume and the orientation", isapprox(vol2, vol0; rtol = 1e-12) && pos2)
    # The segment between the apexes misses the common triangle. The flip is not valid.
    ptsN = [0.0 0.0 0.0; 1.0 0.0 0.0; 0.0 1.0 0.0; 0.1 0.1 1.0; 2.0 2.0 -1.0]
    cellsN = [1 2 3 4; 1 3 2 5]
    why = zeros(Int, 4)
    outN, nTriedN, nDoneN = WCore._flip_pass3(ptsN, cellsN, deg2rad(89.0), 2.0, -Inf, why)
    check("flip unit: a flip that does not fill the same space is rejected", nDoneN == 0 && outN == cellsN && why[1] >= 1)
    # No cell is below the goal. Nothing happens.
    outZ, _, nDoneZ = WCore._flip_pass3(pts, cells, deg2rad(1.0), 2.0, 0.0)
    check("flip unit: no flip when no cell is below the goal", nDoneZ == 0 && outZ == cells)
end

# Invariants on real meshes. The comparison uses no smoothing, so only the flips change the mesh.
for (name, L, per, obs, h, hObs, g) in (("box", (1.0, 1.0, 1.0), (false, false, false), [], 0.2, 0.2, 0.3),
                                        ("sphere", (1.0, 1.0, 1.0), (false, false, false), [(1, 0.5, 0.5, 0.5, 0.2, 0.2, 0.2)], 0.2, 0.1, 0.3),
                                        ("periodic two obstacles", (2.0, 1.0, 1.0), (false, true, true),
                                         [(1, 0.7, 0.5, 0.5, 0.2, 0.2, 0.2), (2, 1.4, 0.5, 0.5, 0.15, 0.2, 0.15)], 0.2, 0.1, 0.3))
    m0 = build3(L, per, obs, h, hObs, g; sweeps = 0)
    m1 = build3(L, per, obs, h, hObs, g; sweeps = 0, minDihedral = 10.0)
    pts0, cells0 = m0[1], m0[2]
    pts1, cells1 = m1[1], m1[2]
    v0, p0 = volume_info(pts0, cells0)
    v1, p1 = volume_info(pts1, cells1)
    check("flip $name: the vertices are unchanged", pts0 == pts1 && m0[5] == m1[5])
    check("flip $name: the boundary faces are unchanged", boundary_faces(cells0) == boundary_faces(cells1))
    check("flip $name: the total volume is unchanged", isapprox(v0, v1; rtol = 1e-12))
    check("flip $name: every cell has a positive volume", p1)
    check("flip $name: every face belongs to 1 or 2 cells", all(length(v) <= 2 for v in values(face_table(cells1))))
    check("flip $name: the smallest dihedral angle does not shrink", m1[6].minDihedralDeg >= m0[6].minDihedralDeg - 1e-9)
    check("flip $name: the number of cells below the goal does not grow", m1[6].nBelowGoal <= m0[6].nBelowStart + 0)
    check("flip $name: the radius-edge ratio stays within the bound", m1[6].maxRatio <= max(m0[6].maxRatio, 2.0) + 1e-12)
    check("flip $name: the flips reduce the number of cells below 10 degrees", m1[6].nBelowGoal < m1[6].nBelowStart || m1[6].nBelowStart == 0)
    check_mesh3("$name with flips", L, per, obs, build3(L, per, obs, h, hObs, g; minDihedral = 10.0); hObs = hObs)
end
check("min_dihedral above 60 is an error", raises(() -> WCore._resolve_mesh_options3(Dict("mesh_size" => 0.1, "min_dihedral" => 61.0))))
check("min_dihedral below 0 is an error", raises(() -> WCore._resolve_mesh_options3(Dict("mesh_size" => 0.1, "min_dihedral" => -1.0))))
check("min_dihedral defaults to 10", WCore._resolve_mesh_options3(Dict("mesh_size" => 0.1)).minDihedral == 10.0)

# ── Random configurations ────────────────────────────────────────────────

# Build random boxes with random obstacles. Check every mesh.
function random_tests3(nCases)
    println("Random configurations (3-D)")
    rng = MersenneTwister(2025)
    nOk = 0
    for t in 1:nCases
        L = (1.0 + 0.5 * rand(rng), 0.8 + 0.4 * rand(rng), 0.8 + 0.4 * rand(rng))
        h = 0.18 + 0.1 * rand(rng)
        hObs = h * (0.4 + 0.6 * rand(rng))
        g = rand(rng) < 0.2 ? 0.0 : 0.1 + 0.5 * rand(rng)
        per = (rand(rng) < 0.4, rand(rng) < 0.4, rand(rng) < 0.4)
        obs = Tuple[]
        for _ in 1:rand(rng, 0:2)
            kind = rand(rng, 1:2)
            hx = 0.1 + 0.1 * rand(rng)
            hy = kind == 1 ? hx : 0.1 + 0.1 * rand(rng)
            hz = kind == 1 ? hx : 0.1 + 0.1 * rand(rng)
            cx = hx + 0.05 + (L[1] - 2 * hx - 0.1) * rand(rng)
            cy = hy + 0.05 + (L[2] - 2 * hy - 0.1) * rand(rng)
            cz = hz + 0.05 + (L[3] - 2 * hz - 0.1) * rand(rng)
            ok = true
            for o in obs
                sep = maximum([o[2] - o[5] - (cx + hx), cx - hx - (o[2] + o[5]), o[3] - o[6] - (cy + hy),
                               cy - hy - (o[3] + o[6]), o[4] - o[7] - (cz + hz), cz - hz - (o[4] + o[7])])
                sep <= 0.05 && (ok = false)
            end
            ok && push!(obs, (kind, cx, cy, cz, hx, hy, hz))
        end
        name = "random $t"
        local mesh
        try
            mesh = build3(L, per, obs, h, hObs, g)
        catch e
            check("$name: the mesh builds ($(sprint(showerror, e)))", false)
            continue
        end
        before = nFail
        check_mesh3(name, L, per, obs, mesh; hObs = hObs)
        nFail == before && (nOk += 1)
    end
    println("  random cases without a failed check: $nOk of $nCases")
end
random_tests3(parse(Int, get(ENV, "N_RANDOM", "12")))

# ── Option and obstacle checks ───────────────────────────────────────────

println("Option and obstacle checks (3-D)")
check("mesh_size is required", raises(() -> WCore._resolve_mesh_options3(Dict("grading" => 0.1))))
check("unknown mesh_options key is an error", raises(() -> WCore._resolve_mesh_options3(Dict("mesh_size" => 0.1, "foo" => 1))))
check("radius_edge_bound below 2 is an error", raises(() -> WCore._resolve_mesh_options3(Dict("mesh_size" => 0.1, "radius_edge_bound" => 1.5))))
box3 = (2.0, 1.0, 1.0)
check("obstacle that crosses the box is an error",
      raises(() -> WCore._validate_obstacles3(IOBuffer(), [Dict("shape" => "sphere", "center" => [0.1, 0.5, 0.5], "radius" => 0.2)], box3, 0.1)))
check("overlapping bounding boxes are an error",
      raises(() -> WCore._validate_obstacles3(IOBuffer(), [Dict("shape" => "sphere", "center" => [0.8, 0.5, 0.5], "radius" => 0.2),
                                               Dict("shape" => "box", "center" => [1.1, 0.5, 0.5], "half_extents" => [0.2, 0.2, 0.2])], box3, 0.1)))
check("a 2-D shape is an error in 3-D",
      raises(() -> WCore._validate_obstacles3(IOBuffer(), [Dict("shape" => "circle", "center" => [1.0, 0.5, 0.5], "radius" => 0.2)], box3, 0.1)))
check("a valid obstacle passes",
      !raises(() -> WCore._validate_obstacles3(IOBuffer(), [Dict("shape" => "sphere", "center" => [1.0, 0.5, 0.5], "radius" => 0.2)], box3, 0.1)))

# ── Helper tests for all dimensions ──────────────────────────────────────

println("Helper tests (2-D and 3-D)")
# The slow reference version of the edge registry (a Dict).
function reference_registry(cells, localEdges)
    edgeIndex = Dict{Tuple{Int,Int},Int}()
    rows = Tuple{Int,Int}[]
    for ci in 1:size(cells, 1), (a, b) in localEdges
        va, vb = cells[ci, a], cells[ci, b]
        key = va < vb ? (va, vb) : (vb, va)
        if !haskey(edgeIndex, key)
            push!(rows, key)
            edgeIndex[key] = length(rows) - 1
        end
    end
    edges = Matrix{Int}(undef, length(rows), 2)
    for (i, (a, b)) in enumerate(rows)
        edges[i, 1], edges[i, 2] = a, b
    end
    c2e = Matrix{Int}(undef, size(cells, 1), length(localEdges))
    for ci in 1:size(cells, 1), (li, (a, b)) in enumerate(localEdges)
        va, vb = cells[ci, a], cells[ci, b]
        c2e[ci, li] = edgeIndex[va < vb ? (va, vb) : (vb, va)]
    end
    return edges, c2e
end
# The slow reference version of the vertex-to-cell lists.
function reference_v2c(nPts, cells)
    adjacency = [Int[] for _ in 1:nPts]
    for ci in 1:size(cells, 1), k in 1:size(cells, 2)
        push!(adjacency[cells[ci, k]], ci - 1)
    end
    offsets = zeros(Int, nPts + 1)
    for i in 1:nPts
        offsets[i + 1] = offsets[i] + length(adjacency[i])
    end
    return offsets, reduce(vcat, adjacency; init = Int[])
end
# The slow reference version of the periodic search: the nearest vertex
# over all vertices of the high side, as in the old files.
function reference_periodic(pts, box, periodic; bndTol = 1e-9, matchTol = 1e-7)
    pairs = Dict{String,Matrix{Int}}()
    for d in 1:length(periodic)
        periodic[d] || continue
        L, dname = box[d], WCore._dirs()[d]
        other = [k for k in 1:size(pts, 2) if k != d]
        lo = [i for i in 1:size(pts, 1) if abs(pts[i, d]) < bndTol]
        hi = [i for i in 1:size(pts, 1) if abs(pts[i, d] - L) < bndTol]
        matched = Tuple{Int,Int}[]
        for vi in lo
            bestJ, bestD = 0, Inf
            for vj in hi
                dist = norm(pts[vi, other] .- pts[vj, other])
                dist < bestD && (bestD = dist; bestJ = vj)
            end
            bestD < matchTol && push!(matched, (vi - 1, bestJ - 1))
        end
        pairs[dname] = isempty(matched) ? zeros(Int, 0, 2) : reduce(vcat, [reshape([a, b], 1, 2) for (a, b) in matched])
    end
    return pairs
end

# A 2-D mesh: the Delaunay mesh from the 2-D tests, and a grid mesh.
mesh2 = build2(2.0, 1.0, true, false, [(1, 1.0, 0.5, 0.2, 0.2)], 0.12, 0.05, 0.3)
e2, c2 = WCore._build_edge_registry(mesh2[2], 2)
eRef2, cRef2 = reference_registry(mesh2[2], WCore._tri_local_edges())
check("2-D: sorted edge registry equals the Dict registry", e2 == eRef2 && c2 == cRef2)
check("2-D: the registry equals the edges that the mesher built", e2 == mesh2[3])
off2, idx2 = WCore._vertex_to_cell(size(mesh2[1], 1), mesh2[2])
offRef2, idxRef2 = reference_v2c(size(mesh2[1], 1), mesh2[2])
check("2-D: counted vertex-to-cell lists equal the old lists", off2 == offRef2 && idx2 == idxRef2)
pts2 = hcat(mesh2[1], zeros(size(mesh2[1], 1)))
redirect_stdout(devnull) do
    global ppNew2 = WCore._periodic_pairs(IOBuffer(), pts2, (2.0, 1.0), (true, false))
end
ppRef2 = reference_periodic(pts2, (2.0, 1.0), (true, false))
check("2-D: hashed periodic pairs equal the brute-force pairs", ppNew2 == ppRef2 && size(ppNew2["x"], 1) > 5)

meshK = WCore._build_kuhn(IOBuffer(), WCore._check_config(Dict("mesh_type" => "kuhn", "box_length" => [1.0, 1.0, 1.0],
                                                    "n_divisions" => [3, 4, 3], "diagonal" => "alternate_layers")))
eK, cK = WCore._build_edge_registry(meshK.cells, 3)
eRefK, cRefK = reference_registry(meshK.cells, WCore._tet_local_edges())
check("kuhn: sorted edge registry equals the Dict registry", eK == eRefK && cK == cRefK)
offK, idxK = WCore._vertex_to_cell(size(meshK.pts, 1), meshK.cells)
offRefK, idxRefK = reference_v2c(size(meshK.pts, 1), meshK.cells)
check("kuhn: counted vertex-to-cell lists equal the old lists", offK == offRefK && idxK == idxRefK)
redirect_stdout(devnull) do
    global ppNewK = WCore._periodic_pairs(IOBuffer(), meshK.pts, (1.0, 1.0, 1.0), (true, true, true))
end
ppRefK = reference_periodic(meshK.pts, (1.0, 1.0, 1.0), (true, true, true))
check("kuhn: hashed periodic pairs equal the brute-force pairs", ppNewK == ppRefK && size(ppNewK["z"], 1) == 20)
meshE = WCore._build_equilateral(IOBuffer(), WCore._check_config(Dict("mesh_type" => "equilateral", "box_length" => [2.0, 1.0],
                                                          "n_divisions" => [6, 7], "periodic" => [true, true])))
ptsE = hcat(meshE.pts, zeros(size(meshE.pts, 1)))
redirect_stdout(devnull) do
    global ppNewE = WCore._periodic_pairs(IOBuffer(), ptsE, (2.0, 1.0), (true, true))
end
ppRefE = reference_periodic(ptsE, (2.0, 1.0), (true, true))
check("equilateral with odd rows: hashed pairs equal the brute-force pairs", ppNewE == ppRefE)


# ── Equilateral lattice, diagonals, layout (version 2) ───────────────────

# Return the dihedral extremes (degrees) and the radius-edge ratio of each tetrahedron,
# from the kernels that the writer uses.
function quality3(cells, pts)
    n = size(cells, 1)
    minD = zeros(n); maxD = zeros(n)
    WCore.tet_dihedral_kernel(cells, pts, n, minD, maxD)
    ids = collect(1:n)
    vx = pts[:, 1]; vy = pts[:, 2]; vz = pts[:, 3]
    ccx = zeros(n); ccy = zeros(n); ccz = zeros(n); ccR = zeros(n)
    WCore.cell_circumcenter_3d_kernel(cells, vx, vy, vz, ids, n, ccx, ccy, ccz, ccR)
    ratio = zeros(n); sizeRatio = zeros(n)
    WCore.cell_quality_3d_kernel(cells, vx, vy, vz, ids, n, ccR, ones(n), ratio, sizeRatio)
    return rad2deg.(minD), rad2deg.(maxD), ratio
end

println("Cell quality kernels")
# One triangle: right isosceles, then equilateral.
let
    pts = [0.0 0.0; 1.0 0.0; 0.0 1.0; 0.5 sqrt(3) / 2]
    cells = [1 2 3; 1 2 4]
    mn = zeros(2); mx = zeros(2)
    WCore.triangle_angle_kernel(cells, pts, 2, mn, mx)
    check("triangle angles: right isosceles 45 and 90", isapprox(mn[1], 45.0; atol = 1e-12) && isapprox(mx[1], 90.0; atol = 1e-12))
    check("triangle angles: equilateral 60 and 60", isapprox(mn[2], 60.0; atol = 1e-12) && isapprox(mx[2], 60.0; atol = 1e-12))
    cells2 = [1 3 2]
    WCore.triangle_angle_kernel(cells2, pts, 1, mn, mx)
    check("triangle angles: the orientation does not matter", isapprox(mn[1], 45.0; atol = 1e-12) && isapprox(mx[1], 90.0; atol = 1e-12))
end
# One regular tetrahedron: the dihedral angle is acos(1/3) and the ratio is sqrt(6)/4.
let
    pts = [1.0 1.0 1.0; 1.0 -1.0 -1.0; -1.0 1.0 -1.0; -1.0 -1.0 1.0]
    for (label, cells) in (("one orientation", [1 2 3 4]), ("the other orientation", [1 2 4 3]))
        mn, mx, ra = quality3(cells, pts)
        want = acos(1 / 3) * 180 / pi
        check("regular tetrahedron ($label): both dihedral extremes are $(round(want; digits = 4))",
              isapprox(mn[1], want; atol = 1e-9) && isapprox(mx[1], want; atol = 1e-9))
        check("regular tetrahedron ($label): the radius-edge ratio is sqrt(6)/4", isapprox(ra[1], sqrt(6) / 4; atol = 1e-12))
    end
end
# One corner tetrahedron of a cube (the Kuhn shape): 45, 60 and 90 degrees, ratio sqrt(3)/2.
let
    pts = [0.0 0.0 0.0; 1.0 0.0 0.0; 1.0 1.0 0.0; 1.0 1.0 1.0]
    mn, mx, ra = quality3([1 2 3 4], pts)
    check("Kuhn tetrahedron: dihedral 45 to 90 degrees", isapprox(mn[1], 45.0; atol = 1e-9) && isapprox(mx[1], 90.0; atol = 1e-9))
    check("Kuhn tetrahedron: the radius-edge ratio is sqrt(3)/2", isapprox(ra[1], sqrt(3) / 2; atol = 1e-12))
end
# Random tetrahedra against a direct computation: the circumcenter from a linear
# system, the dihedral angle from the projections of two vertices on a plane
# normal to the edge.
let
    rng = MersenneTwister(7)
    n = 200
    pts = rand(rng, n * 4, 3)
    cells = reshape(collect(1:n * 4), 4, n)'
    cells = Matrix{Int}(cells)
    mn, mx, ra = quality3(cells, pts)
    okMin = true; okMax = true; okRatio = true
    for c in 1:n
        X = pts[cells[c, :], :]
        A = 2 .* (X[2:4, :] .- X[1:1, :])
        rhs = vec(sum(X[2:4, :] .^ 2 .- X[1:1, :] .^ 2; dims = 2))
        radius = norm(A \ rhs - X[1, :])
        shortest = minimum(norm(X[i, :] - X[j, :]) for i in 1:4 for j in i + 1:4)
        angles = Float64[]
        for i in 1:4, j in i + 1:4
            k, l = [m for m in 1:4 if m != i && m != j]
            e = (X[j, :] - X[i, :]) / norm(X[j, :] - X[i, :])
            u = X[k, :] - X[i, :]; u -= e * dot(u, e)
            v = X[l, :] - X[i, :]; v -= e * dot(v, e)
            push!(angles, acos(clamp(dot(u, v) / (norm(u) * norm(v)), -1, 1)) * 180 / pi)
        end
        okMin &= isapprox(mn[c], minimum(angles); rtol = 1e-8)
        okMax &= isapprox(mx[c], maximum(angles); rtol = 1e-8)
        okRatio &= isapprox(ra[c], radius / shortest; rtol = 1e-8)
    end
    check("random tetrahedra: the smallest dihedral angle matches a direct computation", okMin)
    check("random tetrahedra: the largest dihedral angle matches a direct computation", okMax)
    check("random tetrahedra: the radius-edge ratio matches a direct computation", okRatio)
end

println("Equilateral lattice and 2-D diagonals")
let
    h = sqrt(3.0) / 2
    # dx = 0.25, so dy = dx*sqrt(3)/2 and the box height is 4*dy = sqrt(3)/2
    cfgE = WCore._check_config(Dict("mesh_type" => "equilateral", "box_length" => [1.0, h], "n_divisions" => [4, 4]))
    mE = WCore._build_equilateral(IOBuffer(), cfgE)
    check("equilateral: an odd row has nx+2 points, so 27 vertices for 4 by 4", size(mE.pts, 1) == 27)
    check("equilateral: ny*(2*nx+1) = 36 cells", size(mE.cells, 1) == 36)
    # Every vertex number is used, and every cell has a positive area.
    check("equilateral: every vertex belongs to a cell", sort(unique(vec(mE.cells))) == collect(1:27))
    angs = Float64[]
    nEq = 0
    for c in 1:size(mE.cells, 1)
        p = [mE.pts[mE.cells[c, k], :] for k in 1:3]
        a = Float64[]
        for k in 1:3
            u = p[mod1(k + 1, 3)] - p[k]; v = p[mod1(k + 2, 3)] - p[k]
            push!(a, acosd(clamp(dot(u, v) / (norm(u) * norm(v)), -1.0, 1.0)))
        end
        append!(angs, a)
        all(x -> abs(x - 60.0) < 1e-9, a) && (nEq += 1)
    end
    check("equilateral: 28 interior triangles are exactly equilateral (8 are right-angled)", nEq == 28)
    check("equilateral: the angles lie between 30 and 90 degrees", minimum(angs) > 30.0 - 1e-9 && maximum(angs) < 90.0 + 1e-9)
    check("equilateral: 16 vertices lie on the box boundary", count(!=(0), mE.ptMk) == 16)
    check("equilateral: 16 boundary edges carry a marker", count(!=(0), mE.edgeMk) == 16)
    # slash and backslash
    for (name, wantSlash) in (("slash", true), ("backslash", false))
        cfg = WCore._check_config(Dict("mesh_type" => "gridsplit", "box_length" => [1.0, 1.0],
                                       "n_divisions" => [1, 1], "diagonal" => name))
        m = WCore._build_gridsplit(IOBuffer(), cfg)
        # the corner (0,0) is vertex 1 and (1,1) is vertex 4; a slash joins them
        pairs = Set(Tuple(sort([m.cells[c, a], m.cells[c, b]])) for c in 1:2 for a in 1:3 for b in a+1:3)
        check("gridsplit: diagonal $name gives the " * (wantSlash ? "/" : "\\") * " cut", ((1, 4) in pairs) == wantSlash)
    end
    cfgA = WCore._check_config(Dict("mesh_type" => "gridsplit", "box_length" => [2.0, 2.0],
                                    "n_divisions" => [2, 2], "diagonal" => "alternate"))
    mA = WCore._build_gridsplit(IOBuffer(), cfgA)
    pairsA = Set(Tuple(sort([mA.cells[c, a], mA.cells[c, b]])) for c in 1:8 for a in 1:3 for b in a+1:3)
    # cell (i,j) = (0,0) has i+j even: backslash joins bottom-right (2) and top-left (4); cell (1,0) odd: slash joins 2 and 6
    check("gridsplit: alternate uses backslash where i+j is even", (2, 4) in pairsA)
    check("gridsplit: alternate uses slash where i+j is odd", (2, 6) in pairsA)
    check("gridsplit: the old diagonal names main and flip are errors",
          raises(() -> WCore._diagonal_mode_2d("main")) && raises(() -> WCore._diagonal_mode_2d("flip")))
    check("gridsplit: the default diagonal is backslash",
          WCore._check_config(Dict("mesh_type" => "gridsplit", "box_length" => [1.0, 1.0], "n_divisions" => [2, 2]))["diagonal"] == "backslash")
    # symmetric stretching with an odd n: the right half has the extra cell
    posS = zeros(4)
    WCore.stretch_positions_kernel(3, 1.5, 1, posS)
    check("stretching: symmetric with odd n puts one cell left of the middle and two right of it", posS[2] == 0.5 && posS[3] > 0.5 && posS[3] < 1.0)
    posS1 = zeros(2)
    WCore.stretch_positions_kernel(1, 1.5, 1, posS1)
    check("stretching: symmetric with n = 1 is a single cell", posS1 == [0.0, 1.0])
end

# ── Config check ─────────────────────────────────────────────────────────

println("Config matrix")
# A valid base config for each mesh type, with only the required keys.
baseConfigs = Dict(
    "gridsplit"   => () -> Dict{Any,Any}("mesh_type" => "gridsplit", "box_length" => [2.0, 1.0], "n_divisions" => [4, 3]),
    "equilateral" => () -> Dict{Any,Any}("mesh_type" => "equilateral", "box_length" => [2.0, 1.0], "n_divisions" => [4, 4]),
    "kuhn"        => () -> Dict{Any,Any}("mesh_type" => "kuhn", "box_length" => [1.0, 1.0, 1.0], "n_divisions" => [2, 2, 2]),
    "triangle"    => () -> Dict{Any,Any}("mesh_type" => "triangle", "box_length" => [2.0, 1.0], "mesh_options" => Dict("mesh_size" => 0.2)),
    "tetra"       => () -> Dict{Any,Any}("mesh_type" => "tetra", "box_length" => [1.0, 1.0, 1.0], "mesh_options" => Dict("mesh_size" => 0.3)),
)
dimOf = Dict("gridsplit" => 2, "equilateral" => 2, "kuhn" => 3, "triangle" => 2, "tetra" => 3)
for mt in WCore._mesh_types()
    check("$mt: the base config passes", errmsg(() -> WCore._check_config(baseConfigs[mt]())) === nothing)
end
# key => (value for dimension 2, value for dimension 3, set of types that accept it)
keyTable = [
    ("n_divisions", [4, 4], [3, 3, 3], ("gridsplit", "equilateral", "kuhn")),
    ("diagonal", "backslash", "main", ("gridsplit", "kuhn")),
    ("stretching_ratio", [1.1, 1.1], [1.1, 1.1, 1.1], ("gridsplit", "kuhn")),
    ("stretching_symmetric", [false, false], [false, false, false], ("gridsplit", "kuhn")),
    ("obstacles", [Dict("shape" => "circle", "center" => [1.0, 0.5], "radius" => 0.1)],
                  [Dict("shape" => "sphere", "center" => [0.5, 0.5, 0.5], "radius" => 0.1)], ("triangle", "tetra")),
    ("mesh_options", Dict("mesh_size" => 0.2), Dict("mesh_size" => 0.2), ("triangle", "tetra")),
    ("gmsh_options", Dict("a" => 1), Dict("a" => 1), ()),
    ("periodic", [false, false], [false, false, false], WCore._mesh_types()),
]
for mt in WCore._mesh_types(), (key, v2, v3, accepted) in keyTable
    cfg = baseConfigs[mt]()
    cfg[key] = dimOf[mt] == 2 ? v2 : v3
    msg = errmsg(() -> WCore._check_config(cfg))
    if mt in accepted
        check("$mt accepts '$key'", msg === nothing)
    else
        check("$mt rejects '$key' and names the key and the type", msg !== nothing && occursin("'$key'", msg) && occursin(mt, msg))
    end
end
# An empty obstacles list is no error for a structured type.
cfgEmpty = baseConfigs["gridsplit"](); cfgEmpty["obstacles"] = []
check("a structured type ignores an empty obstacles list", errmsg(() -> WCore._check_config(cfgEmpty)) === nothing)
# Required keys.
for mt in WCore._mesh_types()
    cfg = baseConfigs[mt](); delete!(cfg, "box_length")
    check("$mt: box_length is required", raises(() -> WCore._check_config(cfg)))
    cfg = baseConfigs[mt](); cfg["box_length"] = [1.0]
    check("$mt: box_length needs the right number of entries", raises(() -> WCore._check_config(cfg)))
end
for mt in ("gridsplit", "equilateral", "kuhn")
    cfg = baseConfigs[mt](); delete!(cfg, "n_divisions")
    check("$mt: n_divisions is required", raises(() -> WCore._check_config(cfg)))
end
for mt in ("triangle", "tetra")
    cfg = baseConfigs[mt](); delete!(cfg, "mesh_options")
    check("$mt: mesh_options is required", raises(() -> WCore._check_config(cfg)))
end
# Options of the other dimension are unknown keys.
cfgT = baseConfigs["triangle"](); cfgT["mesh_options"]["radius_edge_bound"] = 2.0
check("triangle rejects radius_edge_bound", occursin("unknown key", something(errmsg(() -> WCore._check_config(cfgT)), "")))
cfgR = baseConfigs["tetra"](); cfgR["mesh_options"]["min_angle"] = 20.0
check("tetra rejects min_angle", occursin("unknown key", something(errmsg(() -> WCore._check_config(cfgR)), "")))
# Defaults of the Delaunay options differ by dimension, as in the old files.
check("triangle default max_vertices is 2000000", WCore._check_config(baseConfigs["triangle"]())["mesh_options"].maxVertices == 2000000)
check("tetra default max_vertices is 1000000", WCore._check_config(baseConfigs["tetra"]())["mesh_options"].maxVertices == 1000000)
# The default type, and the missing type.
println("Default and missing type")
cfgDef = WCore._check_config(Dict{Any,Any}("box_length" => [1.0, 1.0], "n_divisions" => [3, 3]))
check("no mesh_type and no mesh_options gives gridsplit", cfgDef["mesh_type"] == "gridsplit" && cfgDef["dim"] == 2)
msgMissing = errmsg(() -> WCore._check_config(Dict{Any,Any}("box_length" => [1.0, 1.0], "mesh_options" => Dict("mesh_size" => 0.1))))
check("mesh_options without mesh_type is an error that tells the user to set mesh_type",
      msgMissing !== nothing && occursin("mesh_type", msgMissing))
msgBad = errmsg(() -> WCore._check_config(Dict{Any,Any}("mesh_type" => "hexa", "box_length" => [1.0, 1.0])))
check("an unknown mesh_type lists the five valid values",
      msgBad !== nothing && all(occursin(t, msgBad) for t in WCore._mesh_types()))

# ── Names ────────────────────────────────────────────────────────────────

println("Duplicate names")
# List every top-level function (name, number of arguments) in the text.
function function_signatures(src)
    sigs = Tuple{String,Int}[]
    lines = split(src, '\n')
    i = 1
    while i <= length(lines)
        m = match(r"^function\s+([A-Za-z_][A-Za-z0-9_!]*)\s*\(", lines[i])
        if m !== nothing
            txt = lines[i]
            j = i
            while count(==('('), txt) > count(==(')'), txt) && j < length(lines)
                j += 1
                txt *= "\n" * lines[j]
            end
            a = findfirst('(', txt)
            depth = 0
            k = a
            while k <= lastindex(txt)
                ch = txt[k]
                ch == '(' && (depth += 1)
                if ch == ')'
                    depth -= 1
                    depth == 0 && break
                end
                k = nextind(txt, k)
            end
            inner = txt[nextind(txt, a):prevind(txt, k)]
            d = 0
            n = 0
            seen = false
            for ch in inner
                ch in "([{" && (d += 1)
                ch in ")]}" && (d -= 1)
                ch == ';' && d == 0 && break
                ch == ',' && d == 0 && (n += 1)
                isspace(ch) || (seen = true)
            end
            push!(sigs, (String(m.captures[1]), seen ? n + 1 : 0))
        end
        i += 1
    end
    return sigs
end
sigs = function_signatures(text)
check("the scan finds many functions", length(sigs) > 150)
dupSigs = [s for s in unique(sigs) if count(==(s), sigs) > 1]
check("no function name appears twice with the same number of arguments ($(dupSigs))", isempty(dupSigs))
consts = [m.captures[1] for m in eachmatch(r"(?m)^(?:const\s+)?([A-Za-z_][A-Za-z0-9_]*)\s*=[^=]", text)]
dupConsts = [c for c in unique(consts) if count(==(c), consts) > 1]
check("no top-level constant appears twice ($(dupConsts))", isempty(dupConsts))
allNames = Set(n for (n, _) in sigs)
check("the 2-D mesher functions end with 2",
      all(n in allNames for n in ("_new_mesher2", "_locate2", "_cavity2!", "_commit2!", "_insert2!", "_triangulate_domain2",
                                  "_resolve_mesh_options2", "_validate_obstacles2", "_obstacle_record2", "_max_ratio2")))
check("the 3-D mesher functions end with 3",
      all(n in allNames for n in ("_new_mesher3", "_locate3", "_cavity3!", "_commit3!", "_insert3!", "_triangulate_domain3",
                                  "_resolve_mesh_options3", "_validate_obstacles3", "_obstacle_record3", "_max_ratio3")))
check("the old unsuffixed 2-D names are gone", !any(n in ("_locate", "_cavity!", "_commit!", "_new_mesher", "_triangulate_domain",
                                                         "_resolve_mesh_options", "_validate_obstacles", "_obstacle_record", "_max_ratio",
                                                         "circumcenter_kernel", "cell_quality_kernel", "segment_encroach_kernel",
                                                         "split_point_kernel", "blend_positions_kernel", "point_markers_2d_kernel",
                                                         "point_markers_3d_kernel", "_cell_to_edge", "_build_edge_registry3",
                                                         "_vertex_to_cell3") for (n, _) in sigs))

# ── Timing ───────────────────────────────────────────────────────────────

println("Timing")
for h in (0.02, 0.007)
    t = @elapsed m = build2(1.0, 1.0, false, false, [(1, 0.5, 0.5, 0.15, 0.15)], h, h / 2, 0.3)
    println("  2-D h = $h: $(size(m[2], 1)) cells in $(round(t, digits = 2)) s")
end
for h in (0.1, 0.05)
    t = @elapsed m = build3((1.0, 1.0, 1.0), (false, false, false), [(1, 0.5, 0.5, 0.5, 0.15, 0.15, 0.15)], h, h / 2, 0.3)
    println("  3-D h = $h: $(size(m[2], 1)) cells in $(round(t, digits = 2)) s, exact insphere share $(round(100 * m[6].exactShare, digits = 1)) percent")
end

# ── End-to-end tests ─────────────────────────────────────────────────────

haveIO = try
    @eval using YAML
    @eval using HDF5
    true
catch
    false
end
if haveIO
    println("End-to-end tests")
    W = Module(:W)
    Base.include(W, writerPath)
    dir = mktempdir()

    # One small config for each mesh type. The fields: name, YAML text,
    # dimension, expected HDF5 mesh_type, exact box volume, exact volume
    # of the obstacles that the polyhedra keep exactly (rectangle, box),
    # and an upper bound for the volume that the round obstacles remove.
    cases = [
        ("gridsplit", "mesh_type: gridsplit\nbox_length: [2.0, 1.0]\nn_divisions: [6, 4]\nperiodic: [true, false]\n",
         2, "triangle", 2.0, 0.0, 0.0, false),
        ("equilateral", "mesh_type: equilateral\nbox_length: [2.0, 1.0]\nn_divisions: [6, 6]\nperiodic: [false, true]\n",
         2, "triangle", 2.0, 0.0, 0.0, false),
        ("kuhn", "mesh_type: kuhn\nbox_length: [1.0, 2.0, 1.0]\nn_divisions: [3, 4, 3]\nperiodic: [true, false, true]\n",
         3, "tetra", 2.0, 0.0, 0.0, false),
        ("triangle", "mesh_type: triangle\nbox_length: [3.0, 1.0]\nperiodic: [true, true]\nobstacles:\n" *
         "  - shape: circle\n    center: [1.0, 0.5]\n    radius: 0.2\n  - shape: rectangle\n    center: [2.0, 0.5]\n    half_extents: [0.25, 0.15]\n" *
         "mesh_options:\n  mesh_size: 0.1\n  mesh_size_obstacle: 0.04\n",
         2, "triangle", 3.0, 4 * 0.25 * 0.15, pi * 0.2^2, true),
        ("tetra", "mesh_type: tetra\nbox_length: [2.0, 1.0, 1.0]\nperiodic: [false, true, true]\nobstacles:\n" *
         "  - shape: sphere\n    center: [0.7, 0.5, 0.5]\n    radius: 0.2\n  - shape: box\n    center: [1.4, 0.5, 0.5]\n    half_extents: [0.15, 0.2, 0.15]\n" *
         "mesh_options:\n  mesh_size: 0.2\n  mesh_size_obstacle: 0.1\n",
         3, "tetra", 2.0, 8 * 0.15 * 0.2 * 0.15, 4 / 3 * pi * 0.2^3, true),
        ("triangle_plain", "mesh_type: triangle\nbox_length: [2.0, 1.0]\nmesh_options:\n  mesh_size: 0.25\n",
         2, "triangle", 2.0, 0.0, 0.0, false),
        ("tetra_plain", "mesh_type: tetra\nbox_length: [1.0, 1.0, 1.0]\nmesh_options:\n  mesh_size: 0.34\n",
         3, "tetra", 1.0, 0.0, 0.0, false),
    ]
    # Read every dataset and attribute of an HDF5 file into one Dict.
    function read_all(path)
        out = Dict{String,Any}()
        function walk(g, prefix)
            for a in keys(HDF5.attrs(g))
                out["@$prefix/$a"] = HDF5.read_attribute(g, a)
            end
            for k in keys(g)
                o = g[k]
                if o isa HDF5.Dataset
                    out["$prefix/$k"] = W._read_ds(g, k)
                else
                    walk(o, "$prefix/$k")
                end
            end
        end
        HDF5.h5open(path, "r") do f
            walk(f, "")
        end
        return out
    end
    requiredSets = ("points", "cell_to_vertex", "cell_to_edge", "cell_vol", "cell_vertex_normal",
                    "vertex_to_vertex", "vertex_vol", "vertex_to_cell/offsets", "vertex_to_cell/indices",
                    "boundary/point_markers")
    for (name, yamlText, dim, h5type, boxVol, exactObs, roundObs, hasObs) in cases
        yaml = joinpath(dir, "$name.yaml")
        write(yaml, yamlText)
        redirect_stdout(devnull) do
            W.main(yaml)
        end
        d1 = read_all(joinpath(dir, "$name.h5"))
        # A run twice gives the same file content.
        yaml2 = joinpath(dir, "$(name)_again.yaml")
        write(yaml2, yamlText)
        redirect_stdout(devnull) do
            W.main(yaml2)
        end
        d2 = read_all(joinpath(dir, "$(name)_again.h5"))
        check("$name: a run twice gives the same HDF5 content", d1 == d2)

        # Schema (section 5 of the plan).
        for s in requiredSets
            check("$name: dataset $s exists", haskey(d1, "/$s"))
        end
        periodicAll = any(W.load_config(yaml)["periodic"])
        check("$name: all_pairs exists exactly when an axis is periodic", haskey(d1, "/periodicity/all_pairs") == periodicAll)
        check("$name: mesh_type attribute is $h5type", d1["@/mesh_type"] == h5type)
        check("$name: schema_version attribute is 1", d1["@/schema_version"] == 1)
        qNames = dim == 2 ? ["min_angle_deg", "max_angle_deg"] : ["min_dihedral_deg", "max_dihedral_deg", "max_radius_edge_ratio"]
        qOther = dim == 2 ? ["min_dihedral_deg", "max_dihedral_deg", "max_radius_edge_ratio"] : ["min_angle_deg", "max_angle_deg"]
        check("$name: the quality attributes $(join(qNames, ", ")) exist", all(q -> haskey(d1, "@/$q"), qNames))
        check("$name: the quality attributes of the other dimension are absent", !any(q -> haskey(d1, "@/$q"), qOther))
        if dim == 2
            check("$name: the angles are in (0, 180) and ordered",
                  0 < d1["@/min_angle_deg"] <= 60 + 1e-9 && 60 - 1e-9 <= d1["@/max_angle_deg"] < 180)
        else
            check("$name: the dihedral angles and the ratio are sensible",
                  0 < d1["@/min_dihedral_deg"] <= acos(1 / 3) * 180 / pi + 1e-9 &&
                  acos(1 / 3) * 180 / pi - 1e-9 <= d1["@/max_dihedral_deg"] < 180 &&
                  d1["@/max_radius_edge_ratio"] >= sqrt(6) / 4 - 1e-9)
        end
        check("$name: dimension attribute is $dim", d1["@/dimension"] == dim)
        check("$name: box_length attribute has $dim entries", length(d1["@/box_length"]) == dim)
        check("$name: edge markers exist only for a 2-D mesh", haskey(d1, "/boundary/edge_markers") == (dim == 2))
        check("$name: the obstacles group exists only with obstacles", any(startswith(k, "/obstacles") || startswith(k, "@/obstacles") for k in keys(d1)) == hasObs)
        cells = d1["/cell_to_vertex"]
        nPts = d1["@/n_vertices"]
        check("$name: connectivity is 0-based with $(dim + 1) columns",
              minimum(cells) == 0 && maximum(cells) == nPts - 1 && size(cells, 2) == dim + 1)
        check("$name: cell-to-edge has $(dim == 2 ? 3 : 6) columns", size(d1["/cell_to_edge"], 2) == (dim == 2 ? 3 : 6))
        check("$name: points have 3 columns", size(d1["/points"], 2) == 3)
        check("$name: edge array is 0-based", minimum(d1["/vertex_to_vertex"]) == 0 && maximum(d1["/vertex_to_vertex"]) == nPts - 1)
        for (d, dn) in enumerate(("x", "y", "z")[1:dim])
            check("$name: periodic_$dn attribute exists", haskey(d1, "@/periodicity/periodic_$dn"))
        end

        # Pipeline invariants (section 12.3.3 of the plan).
        vol = d1["/cell_vol"]
        check("$name: all cell volumes are positive", all(>(0), vol))
        total = sum(vol)
        if roundObs == 0.0
            check("$name: cell volumes sum to the box volume minus the obstacles",
                  abs(total - (boxVol - exactObs)) < 1e-9 * boxVol)
        else
            # Round obstacles are polyhedra inside the true shape, so the
            # removed volume is below the exact one.
            check("$name: cell volumes sum to a value between the exact bounds",
                  boxVol - exactObs - roundObs - 1e-9 < total < boxVol - exactObs + 1e-9)
        end
        normals = d1["/cell_vertex_normal"]
        check("$name: the vertex normals of each cell sum to zero", maximum(abs.(sum(normals, dims = 2))) < 1e-10)
        vv = d1["/vertex_vol"]
        periodicAny = any(d1["@/periodicity/periodic_$dn"] for dn in ("x", "y", "z")[1:dim])
        if !periodicAny
            check("$name: vertex volumes sum to the cell volumes", abs(sum(vv) - total) < 1e-9 * boxVol)
        end
        pts = d1["/points"]
        Ls = d1["@/box_length"]
        for (d, dn) in enumerate(("x", "y", "z")[1:dim])
            d1["@/periodicity/periodic_$dn"] || continue
            pr = d1["/periodicity/$(dn)_pairs"]
            other = [k for k in 1:3 if k != d]
            lo = count(v -> abs(pts[v, d]) < 1e-9, 1:nPts)
            hi = count(v -> abs(pts[v, d] - Ls[d]) < 1e-6, 1:nPts)
            okPairs = size(pr, 1) == lo && lo == hi &&
                      all(abs(pts[pr[r, 1] + 1, k] - pts[pr[r, 2] + 1, k]) < 1e-7 for r in 1:size(pr, 1), k in other)
            check("$name: periodic axis $dn has matching faces", okPairs)
            check("$name: periodic vertex volumes agree on axis $dn", all(abs(vv[pr[r, 1] + 1] - vv[pr[r, 2] + 1]) < 1e-12 for r in 1:size(pr, 1)))
        end
        if hasObs
            first = dim == 2 ? 10 : 100
            check("$name: the obstacle markers are $first and $(first + 1)",
                  d1["@/obstacles/obstacle_0/marker"] == first && d1["@/obstacles/obstacle_1/marker"] == first + 1)
            check("$name: the markers are present in the vertex markers", first in d1["/boundary/point_markers"] && first + 1 in d1["/boundary/point_markers"])
        end

        # XDMF (section 12.3.5 of the plan).
        xdmf = read(joinpath(dir, "$name.xdmf"), String)
        check("$name: the XDMF topology is $(dim == 2 ? "Triangle" : "Tetrahedron")",
              occursin("TopologyType=\"$(dim == 2 ? "Triangle" : "Tetrahedron")\"", xdmf))
        check("$name: the top XDMF grid is named MeshCollection", occursin("<Grid Name=\"MeshCollection\" GridType=\"Collection\"", xdmf))
        check("$name: the provenance label WRITER.JL SOURCE is present", occursin("WRITER.JL SOURCE", xdmf))
        check("$name: the provenance holds the YAML text", occursin("YAML CONFIG: $name.yaml", xdmf))
    end

    # Markers of a plain grid mesh keep their values.
    dg = read_all(joinpath(dir, "gridsplit.h5"))
    mk = dg["/boundary/point_markers"]
    check("gridsplit: the corner markers are the lowest markers", mk[1] == 1 && mk[end] == 2)
    check("gridsplit: interior vertices have marker 0", count(==(0), mk) == (6 - 1) * (4 - 1))

    # Layout of the HDF5 file (version 2): the stored matrices are row-major,
    # with the shape (N, k) that XDMF, h5py and ParaView expect.
    let
        nV = size(dg["/points"], 1); nC = size(dg["/cell_to_vertex"], 1); nE = size(dg["/vertex_to_vertex"], 1)
        # HDF5.jl shows a row-major (N, k) dataset as a Julia (k, N) array
        check("layout: points are stored as N by 3", HDF5.h5open(joinpath(dir, "gridsplit.h5"), "r") do f
            HDF5.size(f["points"]) == (3, nV) end)
        check("layout: cell_to_vertex is stored as N by 3", HDF5.h5open(joinpath(dir, "gridsplit.h5"), "r") do f
            HDF5.size(f["cell_to_vertex"]) == (3, nC) end)
        xd = read(joinpath(dir, "gridsplit.xdmf"), String)
        check("layout: the XDMF dimensions match the dataset shapes",
              occursin("Dimensions=\"$nC 3\">gridsplit.h5:/cell_to_vertex", xd) &&
              occursin("Dimensions=\"$nV 3\">gridsplit.h5:/points", xd))
        check("layout: _read_ds gives the Julia array (N, k) back", size(dg["/points"]) == (nV, 3) && size(dg["/vertex_to_vertex"]) == (nE, 2))
    end


    # Directory mode (section 12.3.6): all five types and one bad config.
    dirMode = joinpath(dir, "many")
    mkpath(dirMode)
    for (name, yamlText, _...) in cases[1:5]
        write(joinpath(dirMode, "$name.yaml"), yamlText)
    end
    write(joinpath(dirMode, "zz_bad.yaml"), "mesh_type: kuhn\nbox_length: [1.0, 1.0, 1.0]\nn_divisions: [2, 2, 2]\ndiagonal: flip_w\n")
    outFile = joinpath(dir, "many.out")
    errFile = joinpath(dir, "many.err")
    proc = run(pipeline(ignorestatus(`$(Base.julia_cmd()) $writerPath $dirMode`); stdout = outFile, stderr = errFile))
    outText = read(outFile, String)
    errText = read(errFile, String)
    check("directory mode: the bad config makes the exit code 1", proc.exitcode == 1)
    check("directory mode: the bad config is reported", occursin("[error] zz_bad.yaml", outText) && occursin("diagonal", outText))
    check("directory mode: the good configs are built",
          all(isfile(joinpath(dirMode, "$n.h5")) && isfile(joinpath(dirMode, "$n.xdmf")) for (n, _...) in cases[1:5]))
    check("directory mode: the bad config writes no file", !isfile(joinpath(dirMode, "zz_bad.h5")))
    check("directory mode: the failure line is on the screen", occursin("Completed with 1 failure(s): zz_bad.yaml", errText))

    # Config errors from files.
    bad1 = joinpath(dir, "bad1.yaml")
    write(bad1, "box_length: [1.0, 1.0]\ngmsh_options:\n  Mesh.MeshSizeMax: 0.1\nn_divisions: [2, 2]\n")
    check("gmsh_options is an error in a file", raises(() -> W.load_config(bad1)))
    bad2 = joinpath(dir, "bad2.yaml")
    write(bad2, "box_length: [1.0, 1.0]\nmesh_options:\n  mesh_size: 0.1\n")
    check("an old Delaunay config without mesh_type is an error", raises(() -> W.load_config(bad2)))


    # ── Gallery (version 2) ───────────────────────────────────────────
    # The 17 configs of the gallery. Each one is written, then read back.
    # The tests check the schema, the layout, the numbering and the
    # orientation. These are the properties that ParaView and writer.py
    # rely on. The fields: name, YAML text, dimension, vertex count, cell
    # count, whether the counts are exact (the Delaunay counts only need
    # to lie within 10 percent), and the sum of the cell volumes.
    println("Gallery")
    galleryCases = [
        ("equilateral_plain", "mesh_type: equilateral\nbox_length: [1.0, 0.866]\nn_divisions: [4, 4]\nperiodic: [false, false]", 2, 27, 36, true, 0.866),
        ("gridsplit_alternate", "mesh_type: gridsplit\nbox_length: [1.0, 1.0]\nn_divisions: [4, 4]\nperiodic: [false, false]\ndiagonal: alternate", 2, 25, 32, true, 1.0),
        ("gridsplit_backslash", "mesh_type: gridsplit\nbox_length: [1.0, 1.0]\nn_divisions: [4, 4]\nperiodic: [false, false]\ndiagonal: backslash", 2, 25, 32, true, 1.0),
        ("gridsplit_slash", "mesh_type: gridsplit\nbox_length: [1.0, 1.0]\nn_divisions: [4, 4]\nperiodic: [false, false]\ndiagonal: slash", 2, 25, 32, true, 1.0),
        ("gridsplit_stretched_onesided", "mesh_type: gridsplit\nbox_length: [1.0, 1.0]\nn_divisions: [16, 16]\nperiodic: [false, false]\nstretching_ratio: [1.0, 1.5]\nstretching_symmetric: [false, false]", 2, 289, 512, true, 1.0),
        ("gridsplit_stretched_symmetric", "mesh_type: gridsplit\nbox_length: [1.0, 1.0]\nn_divisions: [16, 16]\nperiodic: [false, false]\nstretching_ratio: [1.0, 1.5]\nstretching_symmetric: [false, true]", 2, 289, 512, true, 1.0),
        ("kuhn_alternate_layers", "mesh_type: kuhn\nbox_length: [1.0, 1.0, 1.0]\nn_divisions: [3, 3, 3]\nperiodic: [false, false, false]\ndiagonal: alternate_layers", 3, 64, 162, true, 1.0),
        ("kuhn_flip_x", "mesh_type: kuhn\nbox_length: [1.0, 1.0, 1.0]\nn_divisions: [3, 3, 3]\nperiodic: [false, false, false]\ndiagonal: flip_x", 3, 64, 162, true, 1.0),
        ("kuhn_flip_y", "mesh_type: kuhn\nbox_length: [1.0, 1.0, 1.0]\nn_divisions: [3, 3, 3]\nperiodic: [false, false, false]\ndiagonal: flip_y", 3, 64, 162, true, 1.0),
        ("kuhn_flip_z", "mesh_type: kuhn\nbox_length: [1.0, 1.0, 1.0]\nn_divisions: [3, 3, 3]\nperiodic: [false, false, false]\ndiagonal: flip_z", 3, 64, 162, true, 1.0),
        ("kuhn_main", "mesh_type: kuhn\nbox_length: [1.0, 1.0, 1.0]\nn_divisions: [3, 3, 3]\nperiodic: [false, false, false]\ndiagonal: main", 3, 64, 162, true, 1.0),
        ("kuhn_stretched_onesided", "mesh_type: kuhn\nbox_length: [1.0, 1.0, 1.0]\nn_divisions: [12, 12, 12]\nperiodic: [false, false, false]\nstretching_ratio: [1.0, 1.5, 1.0]\nstretching_symmetric: [false, false, false]", 3, 2197, 10368, true, 1.0),
        ("kuhn_stretched_symmetric", "mesh_type: kuhn\nbox_length: [1.0, 1.0, 1.0]\nn_divisions: [12, 12, 12]\nperiodic: [false, false, false]\nstretching_ratio: [1.0, 1.5, 1.0]\nstretching_symmetric: [false, true, false]", 3, 2197, 10368, true, 1.0),
        ("tetra_plain", "mesh_type: tetra\nbox_length: [1.0, 1.0, 1.0]\nperiodic: [false, false, false]\nmesh_options:\n  mesh_size: 0.16666667\n", 3, 987, 4623, false, 1.0000000000218279),
        ("triangle_circle", "mesh_type: triangle\nbox_length: [1.0, 1.0]\nperiodic: [false, false]\nmesh_options:\n  mesh_size: 0.16666667\n  min_angle: 20.0\nobstacles:\n  - shape: circle\n    center: [0.5, 0.5]\n    radius: 0.15\n", 2, 94, 144, false, 0.9363603896932107),
        ("triangle_plain", "mesh_type: triangle\nbox_length: [1.0, 1.0]\nperiodic: [false, false]\nmesh_options:\n  mesh_size: 0.16666667\n  min_angle: 20.0\n", 2, 92, 146, false, 1.0),
        ("triangle_rectangle", "mesh_type: triangle\nbox_length: [1.0, 1.0]\nperiodic: [false, false]\nmesh_options:\n  mesh_size: 0.16666667\n  min_angle: 20.0\nobstacles:\n  - shape: rectangle\n    center: [0.5, 0.5]\n    half_extents: [0.15, 0.1]\n", 2, 92, 139, false, 0.94),
    ]
    galDir = joinpath(dir, "gallery")
    mkpath(galDir)
    for (name, yamlText, dim, nV, nC, exactCounts, volSum) in galleryCases
        yaml = joinpath(galDir, "$name.yaml")
        write(yaml, yamlText)
        redirect_stdout(devnull) do
            W.main(yaml)
        end
        g = read_all(joinpath(galDir, "$name.h5"))
        pts = g["/points"]; cells = g["/cell_to_vertex"] .+ 1
        nv = size(pts, 1); nc = size(cells, 1)
        if exactCounts
            check("gallery $name: $nV vertices and $nC cells", nv == nV && nc == nC)
        else
            check("gallery $name: the counts are within 10 percent of $nV and $nC", abs(nv - nV) <= 0.1nV && abs(nc - nC) <= 0.1nC)
        end
        volTol = exactCounts ? 1e-9 : 0.05
        check("gallery $name: the cell volumes sum to about $(round(volSum, digits = 4))", abs(sum(g["/cell_vol"]) - volSum) < volTol)
        # Layout: every dataset the XDMF names has the declared shape.
        # HDF5.size shows the reversed shape of a row-major dataset.
        xd = read(joinpath(galDir, "$name.xdmf"), String)
        nItems = Ref(0); nOk = Ref(0)
        HDF5.h5open(joinpath(galDir, "$name.h5"), "r") do f
            for m in eachmatch(r"Dimensions=\"([0-9 ]+)\"[^>]*>[^<:]*:/([^<]+)</DataItem>", xd)
                nItems[] += 1
                haskey(f, m.captures[2]) || continue
                declared = parse.(Int, split(m.captures[1]))
                nOk[] += (collect(reverse(HDF5.size(f[m.captures[2]]))) == declared)
            end
        end
        check("gallery $name: all $(nItems[]) XDMF datasets have the declared shape", nItems[] > 0 && nItems[] == nOk[])
        # Numbering: edges are sorted pairs, and the cell-to-edge rows name the cell's own edges.
        e = g["/vertex_to_vertex"] .+ 1; c2e = g["/cell_to_edge"] .+ 1
        check("gallery $name: every edge lists its smaller vertex first", all(e[:, 1] .< e[:, 2]))
        localEdges = dim == 2 ? WCore._tri_local_edges() : WCore._tet_local_edges()
        okC2E = true
        for c in 1:nc, (li, (la, lb)) in enumerate(localEdges)
            a, b = minmax(cells[c, la], cells[c, lb])
            (e[c2e[c, li], 1] == a && e[c2e[c, li], 2] == b) || (okC2E = false)
        end
        check("gallery $name: the cell-to-edge rows name the edges of the cell", okC2E)
        off = g["/vertex_to_cell/offsets"]; idx = g["/vertex_to_cell/indices"]
        okV2C = length(off) == nv + 1 && off[1] == 0 && off[end] == length(idx)
        for v in 1:nv
            okV2C || break
            for k in off[v] + 1:off[v + 1]
                (v in cells[idx[k] + 1, :]) || (okV2C = false)
            end
        end
        check("gallery $name: the vertex-to-cell lists agree with the cells", okV2C)
        # Orientation: every cell has a positive signed measure.
        p3 = size(pts, 2) == 3 ? pts : hcat(pts, zeros(nv))
        meas = zeros(nc)
        dim == 2 ? WCore.signed_measure_2d_kernel(cells, p3, nc, meas) : WCore.signed_measure_3d_kernel(cells, p3, nc, meas)
        check("gallery $name: every cell has a positive orientation", all(>(0), meas))
        check("gallery $name: no vertex is unused", length(unique(vec(cells))) == nv)
    end

    # The stored folder gallery_julia/ holds one YAML config and one PNG per
    # gallery case. Every config must load, and the names must match the cases above.
    storedDir = joinpath(dirname(writerPath), "gallery_julia")
    if isdir(storedDir)
        storedNames = sort([splitext(f)[1] for f in readdir(storedDir) if endswith(f, ".yaml")])
        check("gallery_julia: the YAML names are the gallery case names", storedNames == sort([c[1] for c in galleryCases]))
        check("gallery_julia: every YAML has a PNG", all(n -> isfile(joinpath(storedDir, "$n.png")), storedNames))
        loads(n) = try
            W.load_config(joinpath(storedDir, "$n.yaml"))
            true
        catch
            false
        end
        okLoad = all(loads, storedNames)
        check("gallery_julia: every YAML loads as a config", okLoad)
    end

    # Regression against the reference outputs of the old files
    # (section 12.2.1 of the plan), when the folder is given.
    refDir = get(ENV, "WRITER_REF_DIR", "")
    if !isempty(refDir) && isdir(refDir)
        println("Regression against $refDir")
        nCmp = Ref(0)
        for group in readdir(refDir)
            gdir = joinpath(refDir, group)
            isdir(gdir) || continue
            work = joinpath(dir, "reg_$group")
            mkpath(work)
            for f in readdir(gdir)
                endswith(f, ".yaml") && cp(joinpath(gdir, f), joinpath(work, f))
            end
            for f in sort(readdir(work))
                redirect_stdout(devnull) do
                    W.main(joinpath(work, f))
                end
                refH5 = joinpath(gdir, replace(f, ".yaml" => ".h5"))
                isfile(refH5) || continue
                a = read_all(refH5)
                b = read_all(joinpath(work, replace(f, ".yaml" => ".h5")))
                same = keys(a) == keys(b) && all(a[k] isa AbstractArray{Float64} ? reinterpret(UInt64, a[k]) == reinterpret(UInt64, b[k]) : a[k] == b[k] for k in keys(a))
                check("regression: $group/$f equals the reference output", same)
                nCmp[] += 1
            end
        end
        println("  $(nCmp[]) outputs compared")
    else
        println("Regression skipped: set WRITER_REF_DIR to a folder of reference outputs")
    end
else
    println("End-to-end tests skipped: the YAML or HDF5 package is missing")
end

println()
println("passed: $nPass   failed: $nFail")
exit(nFail == 0 ? 0 : 1)
