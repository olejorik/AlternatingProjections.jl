#= 
Alternating Projections methods, first try
- Julia version: 1.1.0
- Author: Oleg Soloviev
- Date: 2019-09-01 =#


"""
AlternatingProjections

Module contains types describing different feasible sets and corresponding projection methods on them.
Being an auto-pedagogical package, it explorers the possibilities, for instance, adding more functionality
by describing a new type of FeasibleSets in a new separate file and tries not to be perfect.
"""
module AlternatingProjections

using LinearAlgebra
using FFTW
import Base.size # to add size method for FeasibleSet's subtypes

"""
Abstract type containing any algorithm
"""
abstract type Algorithm end

# Problems :-)
abstract type Problem end


# Common type of a problem is feasibility problem, for which we need to introduce concept of
# a feasuible set

"""
    FeasibleSet

Abstract type representing a set for feasibility problem. For any set we should be able to take it's representative element.
"""
abstract type FeasibleSet end

getelement(s::FeasibleSet) = error("Don't know how to take an element of $typeof(s)")

abstract type FeasibilityProblem <: Problem end
struct TwoSetsFP <: FeasibilityProblem
    A::FeasibleSet
    B::FeasibleSet
end


"""
    solve(p::Problem,x⁰,alg::Algorithm)

    solve(p::Problem,alg::IterativeAlgorithm; x⁰, ϵ, maxit, keephistory, snapshshots)


Solve problem `p`, using method `alg`. For iterative algorithms the arguments maybe specified 
Optionally keep the error history and the iteration snapshots.



# Examples

```jldoctest
julia>
```
"""
function solve(p::Problem, alg::Algorithm)
    error("Don't know how to solve ", typeof(p), " with method ", typeof(alg))
end


"""
Iterative methods form an important class of algorithms with iteratively adjust solution.

Concrete types of the `IterativeAlgorithm` should contain the initial value, tolerance and maximum number of iterations and 
are used more for convenience. These can be obtained by fucntions `initial`, `tolerance`, `maxit` fuctions. If applied the
abstract types, these fucntions return `missing` and can trigger using of the default values.

In addition, the concrete types contain instructions on whether or not to keeep the history of the convergence and 
some snapshots of the inner state of an iterator. These values are given by functions `keephistory` and `snapshots`.

"""
abstract type IterativeAlgorithm <: Algorithm end

initial(alg::IterativeAlgorithm) = missing
tolerance(alg::IterativeAlgorithm) = missing
maxit(alg::IterativeAlgorithm) = missing
keephistory(alg::IterativeAlgorithm) = false
snapshots(alg::IterativeAlgorithm) = Int64[]

"""
ProjectionsMethod is class of iterative algorithms for solving feasibility problems based on projections on the feasible sets.
    Examples are Alternating Projections (AP), also known in the literature as Projections on the Convex Sets (POCS),
    Douglas-Rachford algorithm (DR), their combination DRAP and many others.

"""
abstract type  ProjectionsMethod <: IterativeAlgorithm end

include("AP.jl")



# Projections

"""
    project!(xp, x, A)

Project `x` on set `A` using preallocated `xp`.
"""
function project!(xp, x, feasset::FeasibleSet)
    error("Don't know how to project on ", typeof(feasset))
end

"""
    project!(x, A)

Project `x` on set `A` and storing result in `x`.
"""
function project!(x, feasset::FeasibleSet)
    error("Don't know how to project on ", typeof(feasset))
end

"""
    project(x, A)

Project `x` on set `A` and return the result of projection.
"""
project(x, feasset::FeasibleSet) = project!(copy(x), x, feasset)

reflect(x, feasset::FeasibleSet) = 2 * project(x, feasset) - x



"""
    ConvexSet

General type, no projection method is specified.

# Examples

```jldoctest
```
"""
abstract type ConvexSet <: FeasibleSet end

include("TransformedSet.jl")




#unpack the values of iterative algorithm
solve(p::Problem, alg::IterativeAlgorithm; x⁰ = initial(alg), ϵ =tolerance(alg), maxit = maxit(alg), keephistory = keephistory(alg), snapshots = snapshots(alg)) = solve(p, alg, x⁰, ϵ, maxit, keephistory, snapshots)
# solve(p::Problem, alg::IterativeAlgorithm) = solve(p, alg, initial(alg), tolerance(alg), maxit(alg), keephistory(alg), snapshots(alg))

# Now proceed with the unpacked default or specified by the user values
function solve(p::Problem, alg::IterativeAlgorithm,  args...) 
    error("Don't know how to solve ", typeof(p), " with method ", typeof(alg))
end


export ProjectionsMethod, FeasibleSet, project, project!, ConvexSet, apsolve,solve, TwoSetsFP, TransformedSet, LinearTransformedSet,
FourierTransformedSet, forward!, backward!, generatingset

# Constraints
include("SupportConstraint.jl")
include("AmplitudeConstraint.jl")

# algortihms
# include("GerchbergSaxton.jl")

# Iterators
include("iterators.jl")



end # module
