"""
    TransformedSet

Set obtained by some transformation from a feasible set (`generatingset`).
Should support `forward!` and `backward!` methods. 
"""
abstract type TransformedSet <: FeasibleSet end

function generatingset(ts::TransformedSet)
error("Cannot find the generating set for $ts.")
end


"""
    forward!(ts::TransformedSet)

Two-argument forward transform assosiated with the set `ts`, `ts = F(s)`: `p ∈ s ⇔ q = F(p) ∈ ts`,
it updates the first argument:
use `forward!(ts)(q,p)` to obtain `q = F(p)`.
"""
function forward!(ts::TransformedSet)
error("The forward transform is not defined for $ts.")
end

"""
    backward!(ts::TransformedSet)

Two-argument forward transform assosiated with the set ``ts`, `ts = F(s)`: `p ∈ s ⇔ q = F(p) ∈ ts`,
it updates the first argument:
use `backward!(ts)(p,q)` to obtain `p: q = F(p)`.
"""
function backward!(ts::TransformedSet)
error("The backward transform is not defined for $ts.")
end

# it's easy to get an element of a transformed set
function getelement(ts::TransformedSet) 
    genel = getelement(generatingset(ts))
    el=copy(genel)
    forward!(ts)(el,genel)
    return el
end


"""
    AbstractLinearTransformedSet

Subtype of TransformedSet where `forward` and `backward` transformations are given by multiplication 
by forward and backward "plans" (i.e. --- precomputed matrices). 
"""
abstract type AbstractLinearTransformedSet <: TransformedSet end

generatingset(s::AbstractLinearTransformedSet) = s.set
forward!(s::AbstractLinearTransformedSet) = ((q, p) -> mul!(q, s.fplan, p))
backward!(s::AbstractLinearTransformedSet) = ((p,q) -> mul!(p, s.bplan, q))
bufer(s::AbstractLinearTransformedSet) = s.bufer

function project!(xp, x, s::AbstractLinearTransformedSet) # we don't want to destroy x
    buf = bufer(s)
    backward!(s)(buf, x)
    project!(buf, s.set)
    forward!(s)(xp, buf)
    return xp
end

function project!(x, s::AbstractLinearTransformedSet) # we can discard the value of x
    buf = bufer(s)
    backward!(s)(buf, x)
    project!(buf, s.set)
    forward!(s)(x, buf)
    return x
end

# concrete types
struct LinearTransformedSet{TS,T,N} <: AbstractLinearTransformedSet 
    set::TS
    fplan::Array{T,N}
    bplan::Array{T,N}
    bufer
end

transform(set::FeasibleSet, fplan, bplan)  = LinearTransformedSet(set, fplan, bplan, getelement(set))

struct FourierTransformedSet{TS,PF,PB} <: AbstractLinearTransformedSet where {TS <: FeasibleSet,PF <: AbstractFFTs.Plan,PB <: AbstractFFTs.Plan}
    set::TS
    fplan::PF
    bplan::PB
    bufer
end

# TODO make a flag for direction
FourierTransformedSet(s::FeasibleSet) = 
    FourierTransformedSet(s, FFTW.plan_ifft(getelement(s)), FFTW.plan_fft(getelement(s)), getelement(s))

# struct UnitaryTransformedSet{TS,T,N} <: AbstractLinearTransformedSet where {TS <: FeasibleSet,T,N}
#     set::TS
#     fplan::Array{T,N}
#     bplan::Array{T,N}
#     bufer
# end

