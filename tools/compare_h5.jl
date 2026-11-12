# Compare every dataset and attribute of two HDF5 trees, bit by bit.
# Usage: julia compare_h5.jl REFDIR NEWDIR [--ref-column-major]
# It compares all *.h5 files of REFDIR with the files of the same name in
# NEWDIR. The flag is for a reference written before schema version 1, whose
# matrices are column-major: the script then permutes the axes of the
# reference arrays first.
using HDF5
function walk(g, prefix, out)
    for a in keys(attrs(g)); out["@$prefix/$a"] = read_attribute(g, a); end
    for k in keys(g)
        o = g[k]
        if o isa HDF5.Dataset; out["$prefix/$k"] = read(o)
        else; walk(o, "$prefix/$k", out); end
    end
end
function load(p, flip)
    d = Dict{String,Any}(); h5open(p) do f; walk(f, "", d); end
    flip && for (k, v) in d
        v isa AbstractArray && ndims(v) >= 2 && (d[k] = permutedims(v, ndims(v):-1:1))
    end
    return d
end
function same(a, b)
    typeof(a) == typeof(b) || return false
    if a isa AbstractArray
        size(a) == size(b) || return false
        return a isa AbstractArray{Float64} ? reinterpret(UInt64, a) == reinterpret(UInt64, b) : a == b
    end
    return a isa Float64 ? reinterpret(UInt64, a) == reinterpret(UInt64, b) : a == b
end
ref, new = ARGS[1], ARGS[2]
flipRef = "--ref-column-major" in ARGS
nbad = 0; n = 0
for f in sort(filter(x -> endswith(x, ".h5"), readdir(ref)))
    global n += 1
    if !isfile(joinpath(new, f)); println("MISSING $f"); global nbad += 1; continue; end
    a, b = load(joinpath(ref, f), flipRef), load(joinpath(new, f), false)
    bad = String[]
    for k in union(keys(a), keys(b))
        (haskey(a, k) && haskey(b, k)) || (push!(bad, "$k only in " * (haskey(a, k) ? "ref" : "new")); continue)
        same(a[k], b[k]) || push!(bad, k)
    end
    if !isempty(bad); println("DIFF $f: ", join(sort(bad), ", ")); global nbad += 1; end
end
println("$n files compared, $nbad differ")
exit(nbad == 0 ? 0 : 1)
