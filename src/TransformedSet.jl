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

function forward(ts::TransformedSet)
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

function backward(ts::TransformedSet)
    return error("The backward transform is not defined for $ts.")
end

# it's easy to get an element of a transformed set
function getelement(ts::TransformedSet)
    genel = getelement(generatingset(ts))
    return forward(ts)(genel)
end

"""
    AbstractLinearTransformedSet

Subtype of TransformedSet where `forward` and `backward` transformations are given
by multiplication
by forward and backward "plans" (i.e. --- precomputed matrices, or fast algorithms
implementing this multiplication).
"""
abstract type AbstractLinearTransformedSet <: TransformedSet end
abstract type AbstractScaledCopiesSet <: AbstractLinearTransformedSet end
abstract type AbstractUnitairyTransformedSet <: AbstractLinearTransformedSet end

# and the abstract plans
abstract type AbstractLTPlan{T,N,M} end
abstract type AbstractSCPlan{T,N,M} <: AbstractLTPlan{T,N,M} end

size(p::AbstractLTPlan) = error("size of plan $(typeof(p)) is not defined")
eltype(p::AbstractLTPlan{T,N,M}) where {T,N,M} = T
# forward and backward sizes of the plans are input and output dimensions of the linear transform
fsize(p::AbstractLTPlan) = error("forward size of plan $(typeof(p)) is not defined")
bsize(p::AbstractLTPlan) = error("backward size of plan $(typeof(p)) is not defined")

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

function forward(s::AbstractLinearTransformedSet)
    return (p -> mul!(similar(p, size(s.fplan)...), s.fplan, p))
end
function backward(s::AbstractLinearTransformedSet)
    return (q -> mul!(similar(bplan.scales[1]), s.bplan, q))
end

getplanelement(p::AbstractLTPlan) = zeros(eltype(p), size(p))
getdomainelement(s::AbstractLinearTransformedSet) = zeros(eltype(fplan(s)), bsize(s))
getimageelement(s::AbstractLinearTransformedSet) = zeros(eltype(fplan(s)), fsize(s))

bufer(s::AbstractLinearTransformedSet) = s.bufer

#sizes of the elements in the domain and the image spaces are given by the plan sizes
fsize(s::AbstractLinearTransformedSet) = fsize(fplan(s))
bsize(s::AbstractLinearTransformedSet) = bsize(fplan(s))

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
project!(x, s::AbstractUnitairyTransformedSet) = backproject!(x, s)
project!(xp, x, s::AbstractUnitairyTransformedSet) = backproject!(xp, x, s)

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
    # norm::Array{T,N}
end

fsize(p::plan_SC) = size(p)
bsize(p::plan_SC{T,N,M}) where {T,N,M} = size(p)[1:N]

function invert(p::plan_SC)
    normarray = similar(p.scales[1])
    for i in eachindex(normarray)
        norm = 0
        for j in eachindex(p.scales)
            norm += abs2(p.scales[j][i])
        end
        normarray[i] = norm != 0 ? norm : 1
    end
    return plan_iSC([conj.(p.scales[i]) ./ normarray for i in CartesianIndices(p.scales)])
end

eltype(p::AbstractSCPlan{T,N,M}) where {T,N,M} = T
size(p::AbstractSCPlan) = (size(p.scales[1])..., size(p.scales)...)
# norming(p::plan_iSC) = p.norm

function LinearAlgebra.mul!(y, p::plan_SC, x)
    for i in CartesianIndices(p.scales)
        for indx in CartesianIndices(x)
            y[indx, i] = p.scales[i][indx] * x[indx]
        end
    end
    return y
end

function LinearAlgebra.mul!(x, p::plan_iSC, y)
    for indx in CartesianIndices(x)
        x[indx] = 0
        for i in CartesianIndices(p.scales)
            x[indx] += p.scales[i][indx] * y[indx, i]
        end
        # x[indx] /= p.norm[indx]
    end
    return x
end

struct ScaledCopies{TS,PF,PB} <:
       AbstractScaledCopiesSet where {TS<:FeasibleSet,PF<:AbstractSCPlan,PB<:AbstractSCPlan}
    set::TS
    fplan::PF
    bplan::PB
    bufer
end

function ScaledCopies(A::FeasibleSet, scales)
    return ScaledCopies(A, plan_SC(scales), invert(plan_SC(scales)), getelement(A))
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
