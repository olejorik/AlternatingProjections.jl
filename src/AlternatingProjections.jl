#=
 Alternating Projections methods, first try
- Julia version: 1.1.0
- Author: Oleg Soloviev
- Date: 2019-09-01
=#


"""
    AlternatingProjections

Module contains types describing different feasible sets and corresponding projection methods on them.
Being an auto-pedagogical package, it explorers the possibilities, for instance, adding more functionality
by describing a new type of FeasibleSets in a new separate file and tries not to be perfect.
"""
module AlternatingProjections

using LinearAlgebra
using FFTW



"""
# abstract type APMethod

# Examples

There should be some examples, for sure. At least of the old good Gerchberg-Saxton.


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


# Now the sets

"""
    FeasibleSet

Abstract type representing a set for feasibility problem.
"""
abstract type FeasibleSet end



"""
    ConvexSet

General type, no projection method is specified.

# Examples

```jldoctest
```
"""
abstract type ConvexSet <: FeasibleSet end

import Base.size # to add size method for FeasibleSet's subtypes


# Can you project on any feasible set?

"""
    project!(xp, x, A)

Project `x` on set `A` using preallocated `xp`.
"""
function project!(xp, x, feasset::FeasibleSet)
    error("Don't know how to project on ", typeof(feasset))
end


"""
    project(x, A)

Project `x` on set `A`.
"""
project(x, feasset::FeasibleSet) = project!(copy(x), x, feasset)




reflect(x, feasset::FeasibleSet) = 2 * project(x, feasset) - x


function _actualise_projection(set)
    function actualised(xp,x)
        project!(xp,x,set)
    end
    return actualised
end

# Problems :-)
abstract type Problem end

struct FeasibilityProblem <: Problem
    A::FeasibleSet
    B::FeasibleSet
    forward
    backward
end



"""
    solve(p::Problem,x⁰,alg::APMethod)

    solve(p::Problem,x⁰,alg::APMethod, keephistory::Bool)

    solve(p::Problem,x⁰,alg::APMethod, keephistory::Bool, snapshots)

Solve problem `p`, using method `alg`. Optionally keep the error history and the iteration snapshots



# Examples

```jldoctest
julia>
```
"""
function solve(p::Problem,x⁰::T,alg::APMethod,keephistory::Bool, snapshots::Vector{Int64}) where T
    error("Don't know how to solve ", typeof(p), " with method ", typeof(alg))
end


solve(p::Problem,x⁰,alg::APMethod)  = solve(p, x⁰, alg, false, Int64[]) 

#TODO This is better than below, finish this approach
function solve(p::FeasibilityProblem, x⁰::T, alg::AP, keephistory::Bool, snapshots::Vector{Int64}) where {T}
    A = p.A
    B = p.B
    forward! = p.forward # TODO not sure if forward! = p.forward! would be better
    backward! = p.backward

    # # this is to avoid materialisation of A and B
    # projectB! = _actualise_projection(B)
    # projectA! = _actualise_projection(A)

    maxit = alg.maxit
    maxϵ =alg.maxϵ

    k = 0
    ϵ = Inf
    
    xᵏ = copy(x⁰)
    xᵏ⁺¹ = similar(xᵏ)
    x̃ᵏ⁺¹ = similar(xᵏ)
    ỹᵏ = similar(xᵏ)
    yᵏ = similar(xᵏ)

    #  preallocated error vector
    err = similar(xᵏ)

    if keephistory
        errhist = Vector{Float64}(undef, maxit)
    else
        errhist = Float64[]
    end

    if length(snapshots) != 0
        xhist =fill(x⁰,length(snapshots))
        j=1
    else
        xhist = T[]
    end

    while k < maxit && ϵ > maxϵ
        forward!(ỹᵏ, xᵏ)
        # projectB!(yᵏ, ỹᵏ)
        project!(yᵏ, ỹᵏ, B)
        backward!(x̃ᵏ⁺¹, yᵏ)
        # projectA!(xᵏ⁺¹ , x̃ᵏ⁺¹)
        project!(xᵏ⁺¹ , x̃ᵏ⁺¹, A)

        err .= xᵏ⁺¹ .- xᵏ
        ϵ = LinearAlgebra.norm(err)
        xᵏ .= xᵏ⁺¹
        k += 1

    #         println(ϵ)
        if keephistory
            errhist[k] = ϵ
        end

        if k ∈ snapshots
            xhist[j] .= xᵏ
             j += 1
        end

    end

    # println("To converge with $ϵ accuracy, it took me $k iterations")
    # if keephistory
    #     if length(snapshots) != 0
    #         return xᵏ, errhist, xhist
    #     else
    #         return xᵏ, errhist, Vector{T, (0,)}
    #     end
    # else
    #     return xᵏ, Vector{Float64, (0,)}, Vector{T, (0,)}
    # end
    return xᵏ, errhist, xhist

end


### Below is the old approach
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

export APMethod, FeasibleSet, project, project!, ConvexSet, apsolve,solve, FeasibilityProblem, AP

# Constraints
include("SupportConstraint.jl")
include("AmplitudeConstraint.jl")

# algortihms
include("GerchbergSaxton.jl")





end # module
