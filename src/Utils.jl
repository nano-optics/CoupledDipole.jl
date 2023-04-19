## utils

"""
expand_grid(; kws...)

Expand all combinations of named arguments

- `kws...`: named parameters

# Examples

```jldoctest
julia> expand_grid(a=[0.1,0.2,0.3], b = [1,2])
6×2 DataFrame
 Row │ a        b     
     │ Float64  Int64 
─────┼────────────────
   1 │     0.1      1
   2 │     0.2      1
   3 │     0.3      1
   4 │     0.1      2
   5 │     0.2      2
   6 │     0.3      2
```

"""
function expand_grid(; kws...)
    names, vals = keys(kws), values(kws)
    return DataFrame(NamedTuple{names}(t) for t in Iterators.product(vals...))
end

"""
pmap_df(p, f, kws...; join = true)

Expand all combinations of named arguments

- `p`: DataFrame-like array of parameters
- `f`: function to apply to each row vector
- `kws...`: extra arguments passed to f
- `join`: logical, join parameters and results

"""
function pmap_df(p, f, kws...; join=true, showprogress=false)
    if showprogress
        tmp = @showprogress map(f, eachrow(p), kws...)
    else
        tmp = map(f, eachrow(p), kws...)
    end
    all = reduce(vcat, tmp, source="id")
    if !join
        return all
    end
    p[!, :id] = 1:nrow(p)
    return DataFrames.leftjoin(p, all, on=:id)
end

## orientation-averaging

@doc raw"""
    quadrature_lgwt(N::Int, a::Real, b::Real)

N-point Gauss-Legendre quadrature over [a,b] interval
- `N`: number of nodes
- `a,b`: bounds

```math
\int_a^b f(x)\,dx=\frac{b-a}{2}\int_{-1}^{1}f\left(\frac{b-a}{2} x + \frac{a+b}{2}\right)\,dx.
```

# Examples

```
julia> quadrature_lgwt(6,0,3)
(nodes = [0.10129572869527204, 0.5081859203006033, 1.1420712208752046, 1.8579287791247954, 2.4918140796993966, 2.898704271304728], weights = [0.25698673856875537, 0.5411423595722078, 0.7018709018590369, 0.7018709018590369, 0.5411423595722078, 0.25698673856875537])
```

"""
function quadrature_lgwt(N, a, b)
    n, w = gausslegendre(N)
    (nodes=(b - a) / 2 * n .+ (a + b) / 2,
        weights=(b - a) / 2 * w)
end

"""
    cubature_sphere(N::Int, method::String)

N-point cubature on the sphere
- `N`: number of nodes
- `method`: cubature method (only 'gl' currently implemented)

Returns a Cubature object containing 2 arrays (N'x3 nodes and N'x1 weights), N'≈N
Note: using array instead of tuple for weights because
we'll use them in a scalar product in orientation-averaging
For nodes there is less of a reason, but it can be convenient to visualise the nodes.

The cubature is normalised by 4π such that a unit integrand approximates 1.

# Examples

```
julia> cubature_sphere(6)
(nodes = SVector{3, Float64}[[0.43625314334650644, 2.1862760354652844, 0.0], [0.43625314334650644, 0.9553166181245093, 0.0], [2.0735107047038173, 2.1862760354652844, 0.0], [2.0735107047038173, 0.9553166181245093, 0.0], [4.2096746024757685, 2.1862760354652844, 0.0], [4.2096746024757685, 0.9553166181245093, 0.0], [5.846932163833079, 2.1862760354652844, 0.0], [5.846932163833079, 0.9553166181245093, 0.0]], weights = [0.08696371128436346, 0.08696371128436346, 0.16303628871563655, 0.16303628871563655, 0.16303628871563655, 0.16303628871563655, 0.08696371128436346, 0.08696371128436346])
```

"""
function cubature_sphere(N, method="gl")
    #might have slightly more than N total points
    rndN = Integer(ceil(sqrt(N / 2.0)))

    φ = quadrature_lgwt(2rndN, 0, 2π) # longitude
    cθ = quadrature_lgwt(rndN, 0, 1)   # cos(colatitude)

    nodes = [
        SVector{3}(a, acos(2b - 1), 0.0) for
        b in cθ.nodes, a in φ.nodes
    ]
    weights = hcat([a * b for b in cθ.weights, a in φ.weights]...)

    (nodes=nodes[:], weights=1 / (2π) * weights[:])
end

function spheroid_ar(a₀, χ)

    a = a₀ * χ^(-1 / 3)
    c = χ * a

    return a, c
end


ellipsoid(origin, size) = (origin[1] / size[1])^2 + (origin[2] / size[2])^2 + (origin[3] / size[3])^2

function is_inside(probe, positions, sizes, ParticleRotations)
    tests = map((p, s, r) -> ellipsoid(r' * (probe - p), s) <= 1, positions, sizes, ParticleRotations)
    overall = reduce(|, tests)
    id = findall(tests)
    return overall, id
end