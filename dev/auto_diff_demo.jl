using LinearAlgebra, ForwardDiff
# using BenchmarkTools

struct MyParams{T}
    # general enough to handle dual numbers
    x::T
    y::T
end

struct MyParams2
    # restricts parameter values to a specific concrete type; not good for autodiff
    x::Float64
    y::Float64
end

A = rand(3,3)
B = rand(3,3)
C = rand(3,3)
D = rand(3,3)

"""
    physics(p)
Do "physics" with parameters `p`.
"""
function physics(p)
    X = A*exp(im*p.x)\(B .+ p.y)
    a = X[:,1]
    b = X[:,2]

    r = imag(dot(a,b))
    return r
end

 physics(MyParams(1.3,2.5))


"objective function suitable for autodiff"
function objective(v::AbstractVector)
    p = MyParams(v[1],v[2])
    r = physics(p)
    return r
end


objective([1.3,2.5])



gradFD = ForwardDiff.gradient(objective, [1.3,2.5])




using ReverseDiff

gradRD = ReverseDiff.gradient(objective, [1.3,2.5])

gradFD ≈ gradRD
@benchmark  ReverseDiff.gradient(objective, $([1.3,2.5]))
# 79 μs

ReverseDiff.gradient(objective_2, [1.3,2.5]) # fails
# fails as well

using Zygote

gradZY = Zygote.gradient(objective, [1.3,2.5])
# gives complex number results - looks nonsense
# funnily, the real parts are the same as above

gradFD ≈ gradZY[1]

@benchmark  Zygote.gradient(objective, $([1.3,2.5]))
# 13 μs

Zygote.gradient(objective_2, [1.3,2.5])
# this does work but still gives the same wrong result

# Let's try to use block matrices constructed with `mortar()` from `BlockArrays.jl`
# see (https://juliaarrays.github.io/BlockArrays.jl/stable/lib/public/#BlockArrays.mortar)
mortar((A,B),(C,D))

using BlockArrays

"""
    physics_blocks(p)
Do "physics" with parameters `p` and block matrices.
"""
function physics_blocks(p)
    # construct block matrices using  complex numbers and nonlinear dependence on parameters
    left_matrix = mortar((A,exp(im*p.x+p.y) .+ B),(B, A .* p.y))
    right_matrix = mortar((C,D),(D',D .- p.y^2))
    # call a linear solver
    X = left_matrix\right_matrix
    return X
end

"objective function suitable for autodiff"
function objective_blocks(v::AbstractVector)
    p = MyParams(v[1],v[2])
    X = physics_blocks(p)
    return norm(X)
end

objective_blocks([1.3,2.5])
f2 = ForwardDiff.gradient(objective_blocks, [1.3,2.5])
r2 = ReverseDiff.gradient(objective_blocks, [1.3,2.5])
f2 ≈ r2 # OK
# so both ForwardDiff and ReverseDiff seem to run fine

z2 = Zygote.gradient(objective_blocks, [1.3,2.5]) # fails
# error message is unhelpful ...

# differentiating a matrix eponential
ForwardDiff.gradient(v -> norm(exp(A*v[1])), [1.3,2.5]) #fails
# apparently differentiating through a matrix exponential is too much for ForwardDiff
ReverseDiff.gradient(v -> norm(exp(A*v[1])), [1.3,2.5]) #fails
# same here
Zygote.gradient(v -> norm(exp(A*v[1])), [1.3,2.5])
# gives a result that looks reasonable (but I have not checked it)
