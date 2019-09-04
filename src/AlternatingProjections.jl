#=
 Alternating Projections methods, first try
- Julia version: 1.1.0
- Author: Oleg Soloviev
- Date: 2019-09-01
=#

"""
    AlternatingProjections

Module contains types describing different feasible sets and corresponding projection methods on them.
"""
module AlternatingProjections

using LinearAlgebra
using FFTW

"""
# abstract type APMethod

# Examples

```jldoctest
julia>
```
"""
abstract type  APMethod end

"""
    FeasibleSet

Abstract type representing a set for feasibility problem.
"""
abstract type FeasibleSet end


import Base.size # to add size method for FeasibleSet's subtypes


"""
    project(x, A)

Project `x` on set `A`.
"""
function project(x, feasset::FeasibleSet)
    error("Don't know how to project on ", typeof(feasset))
end

"""
    ConvexSet

General type, no projection method is specified.

- Julia version:
- Author: Oleg Soloviev
- Date: 2019-09-01

# Examples

```jldoctest
```
"""
abstract type ConvexSet <: FeasibleSet end

"""
    apstep(xᵏ, A::FeasibleSet, B::FeasibleSet, f, b)

Iterate `xᵏ` by ``xᵏ⁺¹  = P_A( b_(P_B(f(xᵏ))))``

"""
function apstep(xᵏ, A::FeasibleSet, B::FeasibleSet, forward, backward)
    ỹᵏ = forward(xᵏ)
    yᵏ = project(ỹᵏ, B)
    x̃ᵏ⁺¹ = backward(yᵏ)
    xᵏ⁺¹ = project(x̃ᵏ⁺¹, A)
end

function apsolve(A, B, ::Type{T}; x⁰=zeros(size(A)), maxit = 20, maxϵ =0.01) where {T<:APMethod}
    alg = T(A,B)
    xprev = x⁰
    x = xprev
    i = 0
    ϵ = Inf

    while i < maxit && ϵ > maxϵ 
        x = apstep(xprev, alg.a, alg.A, alg.forward, alg.backward)
        ϵ = LinearAlgebra.norm(x - xprev)
        xprev = x
#         println(ϵ)
        i += 1
    end

    println("To converge with $ϵ accuracy, it took me $i iterations")
    return x
end

export APMethod, FeasibleSet, project, ConvexSet, apsolve

# Constraints
include("SupportConstraint.jl")
include("AmplitudeConstraint.jl")

# algortihms
include("GerschbergSaxton.jl")





end # module
