"""
TransformedSet

Set obtained by some transformation from a feasible set (`generatingset`).
Should support `forward!` and `bacward!` methods. 
"""
abstract type TransformedSet <: FeasibleSet end

function generatingset(ts::TransformedSet)
error("Cannot find the generating set for $ts.")
end


"""
forward!(ts::TransformedSet)

In-place forward transform assosiated with the set `ts`: `x ∈ ts ⇔ f(x) ∈ ts.set`
Use `forward!(ts)(q,p)` to obtain `q = f(p)`.
"""
function forward!(ts::TransformedSet)
error("The forward transorm is not defined for $ts.")
end

function backward!(ts::TransformedSet)
error("The backward transorm is not defined for $ts.")
end

# it's easy to get an element of a transformed set
function getelement(ts::TransformedSet) 
    el = getelement(generatingset(ts))
    backward!(ts)(getelement(generatingset(ts)))
end


"""
LinearTransformedSet

Subtype of TransformedSet where `forward` and `backward` transformations are given by multiplication 
by forward and backward "plans" (i.e. -- precomputed matrices). 
"""
abstract type LinearTransformedSet <: TransformedSet end

generatingset(s::LinearTransformedSet) = s.set
forward!(s::LinearTransformedSet) = ((xf, x) -> mul!(xf, s.fplan, x))
backward!(s::LinearTransformedSet) = (x -> mul!(x, s.bplan, x))

project!(xp, x, s::LinearTransformedSet) = (forward!(s)(xp, x); project!(xp, s.set); backward!(s)(xp))

struct FourierTransformedSet{TS,PF,PB} <: LinearTransformedSet where {TS <: FeasibleSet,PF <: AbstractFFTs.Plan,PB <: AbstractFFTs.Plan}
set::TS
fplan::PF
bplan::PB
end

