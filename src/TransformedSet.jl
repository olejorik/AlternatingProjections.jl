"""
    TransformedSet

Set obtained by some transformation from a feasible set (`generatingset`).
Should support `forward!` and `backward!` methods. 
"""
abstract type TransformedSet <: FeasibleSet end

function generatingset(ts::TransformedSet)
    return error("Cannot find the generating set for $ts.")
end

"""
    forward!(ts::TransformedSet)

Two-argument forward transform assosiated with the set `ts`, `ts = F(s)`: `p ∈ s ⇔ q = F(p) ∈ ts`,
it updates the first argument:
use `forward!(ts)(q,p)` to obtain `q = F(p)`.
"""
function forward!(ts::TransformedSet)
    return error("The forward transform is not defined for $ts.")
end

"""
    backward!(ts::TransformedSet)

Two-argument forward transform assosiated with the set ``ts`, `ts = F(s)`: `p ∈ s ⇔ q = F(p) ∈ ts`,
it updates the first argument:
use `backward!(ts)(p,q)` to obtain `p: q = F(p)`.
"""
function backward!(ts::TransformedSet)
    return error("The backward transform is not defined for $ts.")
end

# it's easy to get an element of a transformed set
function getelement(ts::TransformedSet)
    genel = getelement(generatingset(ts))
    el = copy(genel)
    forward!(ts)(el, genel)
    return el
end

"""
    AbstractLinearTransformedSet

Subtype of TransformedSet where `forward` and `backward` transformations are given by multiplication 
by forward and backward "plans" (i.e. --- precomputed matrices). 
"""
abstract type AbstractLinearTransformedSet <: TransformedSet end
abstract type AbstractScaledCopiesSet <: AbstractLinearTransformedSet end
abstract type AbstractUnitairyTransformedSet <:AbstractLinearTransformedSet end

# and the abstract plans
abstract type AbstractLTPlan{T,N,M} end
abstract type AbstractSCPlan{T,N,M} <: AbstractLTPlan{T,N,M} end

size(p::AbstractLTPlan) = error("size of plan $(typeof(p)) is not defined")
eltype(p::AbstractLTPlan{T,N,M})  where {T,N,M} = T

generatingset(s::AbstractLinearTransformedSet) = s.set

# the forward and backward transforms of the LinearTransformedSet are given by the linear operations
# and thus can be cacluated by `mul!` from LinearAlgebra.
# In some cases this can be calucalte faster, and `plan` is a structure that contains
# all the required information.
# Base.mul! should be defined for each plan
import LinearAlgebra: mul!
fplan(s::AbstractLinearTransformedSet) = s.fplan
bplan(s::AbstractLinearTransformedSet) = s.bplan
forward!(s::AbstractLinearTransformedSet) = ((q, p) -> mul!(q, s.fplan, p))
backward!(s::AbstractLinearTransformedSet) = ((p, q) -> mul!(p, s.bplan, q))

getplanelement(p::AbstractLTPlan) = zeros(eltype(p),size(p))
getdomainelement(s::AbstractLinearTransformedSet) = getplanelement(fplan(s))
getimageelement(s::AbstractLinearTransformedSet) = getplanelement(bplan(s))

bufer(s::AbstractLinearTransformedSet) = s.bufer

#sizes of the elements in the domain and the image spaces are given by the plan szes
fsize(s::AbstractLinearTransformedSet) = size(fplan(s))
bsize(s::AbstractLinearTransformedSet) = size(bplan(s))

"""
    Calcuate image of the "projection in the orignal space" by transforming back, projecting, and transforming forward.

"""
function backproject!(xp, x, s::AbstractLinearTransformedSet) # we don't want to destroy x
    buf = bufer(s)
    backward!(s)(buf, x)
    project!(buf, s.set)
    forward!(s)(xp, buf)
    return xp
end

function backproject!(x, s::AbstractLinearTransformedSet) # we can discard the value of x
    buf = bufer(s)
    backward!(s)(buf, x)
    project!(buf, s.set)
    forward!(s)(x, buf)
    return x
end

# For unitair trnasforms, the back projections is the way to calculate the projection
project!(x, s::AbstractUnitairyTransformedSet) = backproject!(x,s)
project!(xp, x, s::AbstractUnitairyTransformedSet) = backproject!(xp, x,s)


# concrete types

struct plan_LT{T,N,M} <: AbstractLTPlan{T,N,M}
    matrix::Array{T} # TODO {T, N+M}
    dims_domain::Array{Int,1}
    dims_image::Array{Int,1}
end

eltype(p::plan_LT{T,N}) where {T,N} = T

Base.size(p::plan_LT) = size(p.matrix[p.dims_image])
struct LinearTransformedSet{TS,T,N} <: AbstractLinearTransformedSet
    set::TS
    fplan::Array{T,N}
    bplan::Array{T,N}
    bufer
end

function transform(set::FeasibleSet, fplan, bplan)
    return LinearTransformedSet(set, fplan, bplan, getelement(set))
end

# ScaledCopies
struct plan_SC{T,N,M} <: AbstractSCPlan{T,N,M}
    scales::Array{Array{T,N},M}
end

struct plan_iSC{T,N,M} <: AbstractSCPlan{T,N,M}
    scales::Array{Array{T,N},M}
    norm::Array{T,N}
end

function invert(p::plan_SC)
    normarray = similar(p.scales[1])
    for i in eachindex(normarray)
        norm = 0
        for j in eachindex(p.scales)
            norm += abs2(p.scales[j][i])
        end
        normarray[i]=norm
    end
    return plan_iSC(conj.(p.scales),normarray)
end


eltype(p::AbstractSCPlan{T,N, M}) where {T,N, M} = T
size(p::AbstractSCPlan) = (size(p.scales[1])...,size(p.scales)...)

function LinearAlgebra.mul!(y, p::plan_SC, x)
    for i in CartesianIndices(p.scales)
        for indx in CartesianIndices(x)
            y[indx,i] = p.scales[i][indx] * x[indx]
        end
    end
    return y
end

struct ScaledCopies{TS, PF, PB} <: AbstractScaledCopiesSet where
    {TS<:FeasibleSet, PF<:AbstractSCPlan, PB <: AbstractSCPlan}
    set::TS
    fplan::PF
    bplan::PB
    bufer
end

# Fourier-transformed
struct FourierTransformedSet{TS,PF,PB} <: AbstractUnitairyTransformedSet where {
    TS<:FeasibleSet,PF<:AbstractFFTs.Plan,PB<:AbstractFFTs.Plan
}
    set::TS
    fplan::PF
    bplan::PB
    bufer
end

# TODO make a flag for direction
function FourierTransformedSet(s::FeasibleSet)
    return FourierTransformedSet(
        s, FFTW.plan_ifft(getelement(s)), FFTW.plan_fft(getelement(s)), getelement(s)
    )
end

# struct UnitaryTransformedSet{TS,T,N} <: AbstractLinearTransformedSet where {TS <: FeasibleSet,T,N}
#     set::TS
#     fplan::Array{T,N}
#     bplan::Array{T,N}
#     bufer
# end
