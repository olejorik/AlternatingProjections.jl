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
    AP <: APMethod
Classical Alternating projection method

Project on one set, than on another, stop after fixed number of iterations or if the desired accuracy is achieved
"""
struct AP <: APMethod
    maxit
    maxϵ
end

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

abstract type Problem end

struct FeasibilityProblem <: Problem
    A::FeasibleSet
    B::FeasibleSet
    forward
    backward
end

function solve(p::Problem,x⁰,alg::APMethod)
    error("Don't know how to solve ", typeof(p), " with method ", typeof(alg))
end


#TODO This is better than below, finish this approach
function solve(p::FeasibilityProblem, x⁰, alg::AP)
    A = p.A
    B = p.B
    forward = p.forward
    backward = p.backward
    maxit = alg.maxit
    maxϵ =alg.maxϵ

    k = 0
    xᵏ = x⁰
    ϵ = Inf

    while k < maxit && ϵ > maxϵ
        ỹᵏ = forward(xᵏ)
        yᵏ = project(ỹᵏ, B)
        x̃ᵏ⁺¹ = backward(yᵏ)
        xᵏ⁺¹ = project(x̃ᵏ⁺¹, A)
        ϵ = LinearAlgebra.norm(xᵏ⁺¹ - xᵏ)
        xᵏ = xᵏ⁺¹
    #         println(ϵ)
        k += 1
    end

    println("To converge with $ϵ accuracy, it took me $k iterations")
    return xᵏ

end

#TODO this should replace the apsolve function
function findfeasible(A::FeasibleSet,B::FeasibleSet,forward,backward, alg::APMethod)
end



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

export APMethod, FeasibleSet, project, ConvexSet, apsolve,solve,FeasibilityProblem,AP

# Constraints
include("SupportConstraint.jl")
include("AmplitudeConstraint.jl")

# algortihms
include("GerschbergSaxton.jl")





end # module
