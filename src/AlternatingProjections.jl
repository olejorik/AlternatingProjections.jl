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

function apstep(xᵏ, A::FeasibleSet, B::FeasibleSet, forward, backward)
    ỹᵏ = forward(xᵏ)
    yᵏ = project(ỹᵏ, B)
    x̃ᵏ⁺¹ = backward(yᵏ)
    xᵏ⁺¹ = project(x̃ᵏ⁺¹, A)
end

export APMethod, FeasibleSet, project, ConvexSet

# Constraints
include("SupportConstraint.jl")
include("AmplitudeConstraint.jl")

# algortihms




end # module
