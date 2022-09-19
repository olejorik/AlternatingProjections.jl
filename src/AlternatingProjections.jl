#= 
Alternating Projections methods, first try
- Julia version: 1.1.0
- Author: Oleg Soloviev
- Date: 2019-09-01 =#

"""
AlternatingProjections

Module contains types and functions for formulation and solving feasibility problem
    using projections-based methods.
"""
module AlternatingProjections

using LinearAlgebra
using FFTW
import Base.size # to add size method for FeasibleSet's subtypes

include("AbstractProblems.jl")
export Problem,
    Algorithm,
    IterativeAlgorithm,
    solve,
    solution,
    history,
    lasty,
    errhist,
    disthist,
    distgthist,
    xhist,
    itersteps,
    combinehist

export ProjectionsMethod,
    FeasibleSet,
    FeasibilityProblem,
    project,
    project!,
    reflect,
    reflect!,
    ConvexSet,
    apsolve,
    TwoSetsFP,
    TransformedSet,
    LinearTransformedSet,
    FourierTransformedSet,
    forward!,
    backward!,
    generatingset

# Common type of a problem is feasibility problem, for which we need to introduce concept of
# a feasuible set
"""
    FeasibleSet

Abstract type representing a set for feasibility problem. 
    For any set we should be able to take it's representative element.

For any concrete subtype, we define a project operator (maybe set-valued).

"""
abstract type FeasibleSet end

getelement(s::FeasibleSet) = error("Don't know how to take an element of $(typeof(s))")

"""
Big class of feasibility problems./
"""
abstract type FeasibilityProblem <: Problem end

"""
    TwoSetsFP(A,B)

Two sets feasibility problem: for given sets A and B find x∈A∩B, if A∩B≠∅, or x∈A closest to B in some sense.

"""
struct TwoSetsFP <: FeasibilityProblem
    A::FeasibleSet
    B::FeasibleSet
end

"""
    ConvexSet

General type, no projection method is specified. For a convex set a projection is unique.

"""
abstract type ConvexSet <: FeasibleSet end

# Projections

"""
    project!(xp, x, A)

Project `x` on set `A` using preallocated `xp`.
"""
function project!(xp, x, feasset::FeasibleSet)
    return error(
        "Don't know how to project  $(typeof(x)) on ",
        typeof(feasset),
        "(please implement the two-element project)",
    )
end

"""
    project!(x, A)

Project `x` on set `A` and storing result in `x`.
"""
function project!(x, feasset::FeasibleSet)
    return error(
        "Don't know how to project $(typeof(x)) on ",
        typeof(feasset),
        "(please implement the one-element project)",
    )
end

"""
    project(x, A)

Project `x` on set `A` and return the result of projection.
"""
project(x, feasset::FeasibleSet) = project!(copy(x), x, feasset)

# Reflections
"""
    reflect!(xr, x, A)

Reflect `x` on set `A` using preallocated `xp`, R(x) = 2P(x) -x, where P is projection on the set.
"""
function reflect!(xr, x, feasset::FeasibleSet)
    project!(xr, x, feasset)
    #  xr .= 2 .* xr .- x
    @. xr = 2 * xr - x
    return xr
end

reflect(x, feasset::FeasibleSet) = 2 * project(x, feasset) - x

"""
ProjectionsMethod is class of iterative algorithms for solving feasibility problems based on projections on the feasible sets.
    Examples are Alternating Projections (AP), also known in the literature as Projections on the Convex Sets (POCS),
    Douglas-Rachford algorithm (DR), their combination DRAP and many others.

"""
abstract type ProjectionsMethod <: IterativeAlgorithm end

include("AP.jl")
include("DR.jl")
include("DRAP.jl")

include("TransformedSet.jl")

# Constraints
include("SupportConstraint.jl")
include("AmplitudeConstraint.jl")

# Iterators
include("iterators.jl")

end # module
