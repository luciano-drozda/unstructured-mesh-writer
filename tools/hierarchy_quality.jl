#!/usr/bin/env julia
# Compare the diagonal rules of the 3-D refinement on a few meshes. For
# each mesh and each rule, the script refines the mesh again and again. For
# each level it prints the smallest and largest dihedral angle, the
# largest radius-edge ratio and the smallest mean-ratio shape measure.
# The rules are maxmin (the writer default) and fixed1, fixed2, fixed3.
# Usage: julia --project=. tools/hierarchy_quality.jl
using Random, LinearAlgebra
text = read(joinpath(@__DIR__, "..", "writer.jl"), String)
i0 = first(findfirst("# package-free part: begins", text)); i1 = last(findfirst("# package-free part: ends", text))
W = Module(:WG); Base.include_string(W, "using LinearAlgebra\nusing Random\n" * text[i0:i1])

function level0(cfgDict)
    cfg = W._check_config(cfgDict)
    res = redirect_stdout(devnull) do
        W._build_mesh(IOBuffer(), cfg)
    end
    edges, c2e = res.edges === nothing ? W._build_edge_registry(res.cells, 3) : (res.edges, res.c2e)
    return merge(res, (edges = edges, c2e = c2e))
end
function metrics(res)
    pts = res.pts; cells = res.cells
    minD, maxD, vol, ratio = W._cell_numbers3(pts, cells)
    eta = zeros(size(cells, 1)); W.tet_shape_kernel(cells, pts, size(cells, 1), eta)
    return (minDihedral = rad2deg(minimum(minD)), maxDihedral = rad2deg(maximum(maxD)), maxRatio = maximum(ratio), minEta = cbrt(minimum(eta)), nCells = size(cells, 1))
end
function run(name, cfgDict, nLev)
    println("== $name")
    for rule in (0, 1, 2, 3)
        res = level0(deepcopy(cfgDict))
        line = rule == 0 ? "maxmin " : "fixed$rule "
        for k in 0:nLev - 1
            m = metrics(res)
            line *= " L$k: dih $(round(m.minDihedral, digits = 2))-$(round(m.maxDihedral, digits = 1)) rho $(round(m.maxRatio, digits = 2)) eta $(round(m.minEta, digits = 3)) |"
            k < nLev - 1 && (res = W._refine_mesh3(res, rule)[1])
        end
        println(line)
    end
end
for mode in ("main", "flip_x", "flip_y", "flip_z", "alternate_layers")
    run("kuhn $mode", Dict("mesh_type" => "kuhn", "box_length" => [1.0, 1.0, 1.0], "n_divisions" => [3, 3, 3], "diagonal" => mode), 4)
end
run("kuhn stretched", Dict("mesh_type" => "kuhn", "box_length" => [1.0, 1.0, 1.0], "n_divisions" => [3, 3, 3], "stretching_ratio" => [1.3, 1.5, 1.0], "stretching_symmetric" => [false, true, false]), 4)
run("tetra box", Dict("mesh_type" => "tetra", "box_length" => [1.0, 1.0, 1.0], "mesh_options" => Dict("mesh_size" => 0.4)), 4)
run("tetra sphere", Dict("mesh_type" => "tetra", "box_length" => [1.0, 1.0, 1.0], "mesh_options" => Dict("mesh_size" => 0.4, "mesh_size_obstacle" => 0.2),
    "obstacles" => [Dict("shape" => "sphere", "center" => [0.5, 0.5, 0.5], "radius" => 0.2)]), 3)
run("tetra sphere+box", Dict("mesh_type" => "tetra", "box_length" => [2.0, 1.0, 1.0], "mesh_options" => Dict("mesh_size" => 0.5, "mesh_size_obstacle" => 0.3),
    "obstacles" => [Dict("shape" => "sphere", "center" => [0.7, 0.5, 0.5], "radius" => 0.2), Dict("shape" => "box", "center" => [1.4, 0.5, 0.5], "half_extents" => [0.15, 0.2, 0.15])]), 3)
