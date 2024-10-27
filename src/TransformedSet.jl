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
# struct plan_SC{T,N,M} <: AbstractSCPlan{T,N,M}
#     scales::Array{Array{T,N},M}
#     stackedscales #tmp solution
#     dims::Tuple
#     size::Tuple
# end
struct plan_SC{T,NSS, N,M,P,R,S} <: AbstractSCPlan{T,N,M} where {NSS, P,R,S}
    stackedscales::Array{T,NSS} #tmp solution
    dims_p::NTuple{P,Int}
    dims_r::NTuple{R,Int}
    dims_s::NTuple{S,Int}
    size::NTuple{M,Int}
end

function plan_SC(
    scales::Array{Array{T,L},N}, origin, dims=collect(1:ndims(origin))
) where {T,L,N}
    size(scales[1]) == size(origin)[dims] || error("Sizes mismatch when creating `plan_SC`")
    O = ndims(origin)
    M = O + ndims(scales)
    dims_p = tuple(dims...)
    dims_r = Tuple((length(dims) + 1):ndims(origin))
    dims_s = Tuple(length(dims) .+ (1:ndims(scales)))
    s = (size(origin)..., size(scales)...)
    return plan_SC{T,L+N, O,M,length(dims_p),length(dims_r),length(dims_s)}(
        stack(scales), dims_p, dims_r, dims_s, s
    )
end

scales(p::AbstractSCPlan) = eachslice(p.stackedscales; dims=p.dims_s)
stackedscales(p::AbstractSCPlan) = p.stackedscales

struct plan_iSC{T,NSS, N,M,P,R,S} <: AbstractSCPlan{T,N,M} where {NSS, P,R,S}
    stackedscales::Array{T,NSS}
    dims_p::NTuple{P,Int}
    dims_r::NTuple{R,Int}
    dims_s::NTuple{S,Int}
    size::NTuple{M,Int}
    # norm::Array{T,N}
end

fsize(p::plan_SC) = size(p)
bsize(p::plan_SC{T,NSS, N,M,P,R,S}) where {T,NSS, N,M,P,R,S} = size(p)[1:N]

function invert(p::plan_SC{T,NSS, N,M,P,R,S}) where {T,NSS,N,M,P,R,S}
    normarray = reshape(sum(abs2, p.stackedscales; dims=p.dims_s), p.size[[p.dims_p...]])
    for j in eachindex(normarray)
        if normarray[j] == 0
            normarray[j] = 1
        end
    end

    return plan_iSC{T,P+S, M,N,P,R,S}(
        stack([s ./ normarray for s in scales(p)]), p.dims_p, p.dims_r, p.dims_s, bsize(p)
    )
end

# scales(p::plan_iSC) = p.scales

eltype(p::AbstractSCPlan{T,N,M}) where {T,N,M} = T
size(p::AbstractSCPlan) = p.size
# norming(p::plan_iSC) = p.norm

# function LinearAlgebra.mul!(y, p::plan_SC, x)
#     dimsx = (size(x)..., ones(Int64, length(size(scales(p))))...)
#     dimsscales = ones(Int64, length(dimsx))
#     dimsscales[[p.dims...]] .= size(x)[[p.dims...]]
#     dimsscales[(ndims(x) + 1):end] .= size(scales(p))
#     y .= reshape(x, dimsx) .* reshape(p.stackedscales, (dimsscales...)) # plans have now the stacked scales
#     return y
# end

function LinearAlgebra.mul!(
    y::Array, p::plan_SC{T,NSS,N,M,NP,NR,NS}, x::Array
) where {T,NSS,N,M,NP,NR,NS}
    # size(y) == size(p) || error("size of the destination array does not mathc to the plan size")

    # S = scales(p)
    # P = eachslice(x; dims=p.dims_r)
    # @views for ip in eachindex(IndexCartesian(), P[1]),
    #     ir in eachindex(IndexCartesian(), P),
    #     is in eachindex(IndexCartesian(), S)

    # for  ip in  CartesianIndices(size(p)[1:NP]), ir in CartesianIndices(size(p)[(NP+1):N]), is in CartesianIndices(size(p)[(N+1):M]) #attemt to avoid allocations (but works opposite)

    #     # y[ip,ir,is] = (P[ir][ip]) * (S[is][ip])
    #     # @show (ip, ir, is)
    #    @inbounds y[ip, ir, is] = x[ip, ir] * stackedscales(p)[ip, is] # This is less allocations
    # end

    # outer_multiply!(y, stackedscales(p),x, NP, N, M, size(p))

    # outer_multiply2!(y, stackedscales(p),x, NP, N, M, size(p))

    # size1 = ntuple(i -> i <N+1 ? size(p)[i] : 1 ,M )
    # size2 = ntuple(i -> NP < i <N+1 ? 1 : size(p)[i] ,M)
    # outer_multiply3!(y, stackedscales(p),x, NP, N, M,size1, size2)


    sizep = slicedims(size(p), 1, NP)
    sizer = slicedims(size(p), NP + 1, N)
    sizes = slicedims(size(p), N + 1, M)
    # @show size(p), sizep, sizer, sizes
    outer_multiply5!(y, p.stackedscales, x, sizep, sizer, sizes)

    return y
end

function slicedims(dims::NTuple{n,Int}, dfirst::Int, dlast::Int) where {n}
    return ntuple(i -> dims[i + dfirst - 1], dlast - dfirst + 1)
end

function outer_multiply5!(y, s, x, sizep, sizer, sizes)
    for ip in CartesianIndices(sizep),
        ir in CartesianIndices(sizer),
        is in CartesianIndices(sizes)

        @inline y[ip, ir, is] = x[ip, ir] * s[ip, is]
    end
    return y
end



function outer_multiply!(y, s, x, NP, N, M, sizep)
    for ip in CartesianIndices(sizep[1:NP]),
        ir in CartesianIndices(sizep[(NP + 1):N]),
        is in CartesianIndices(sizep[(N + 1):M])

        @inbounds y[ip, ir, is] = x[ip, ir] * s[ip, is]
    end
    return y
end

function outer_multiply2!(y, s, x, NP, N, M, sizep)
    # size1 = ones(Int,M)
    # size2 = ones(Int,M)
    # size1[1:N] .= sizep[1:N]
    # size2[1:NP] .= sizep[1:NP]
    # size2[(N+1):M] .= sizep[(N+1):M]

    size1 = (i < N + 1 ? sizep[i] : 1 for i in 1:M)
    size2 = (NP < i < N + 1 ? 1 : sizep[i] for i in 1:M)

    # @show size1
    # @show size2
    y .= reshape(x, Tuple(size1)) .* reshape(s, Tuple(size2))
    return y
end

function outer_multiply3!(y, s, x, NP, N, M, size1, size2)
    y .= reshape(x, size1) .* reshape(s, size2)
    return y
end




# function LinearAlgebra.mul!(x, p::plan_iSC, y)
#     indp = zeros(Int64, length(p.dims))
#     for indx in CartesianIndices(x)
#         # indp = CartesianIndex((indx.I[[p.dims...]]...))
#         indp .= indx.I[[p.dims...]]
#         x[indx] = 0
#         for i in CartesianIndices(scales(p))
#             x[indx] += scales(p)[i][indp...] * y[indx, i]
#         end
#         # x[indx] /= p.norm[indx]
#     end
#     return x
# end

function LinearAlgebra.mul!(x, p::plan_iSC{T,NSS,M,N,NP,NR,NS}, y) where {T, NSS, N,M,NP,NR,NS}
    # S = scales(p)
    # P = eachslice(y; dims= N .+ ntuple(i->i, M-N))
    # for ip in eachindex(IndexCartesian(), S[1]),
    #     ir in CartesianIndices(p.size[[p.dims_r...]])

    #     x[ip, ir] = 0
    #     for is in eachindex(IndexCartesian(), S)
    #         x[ip, ir] += y[ip, ir, is] * p.stackedscales[ip, is]
    #     end
    # end

    sizep = slicedims(size(p), 1, NP)
    sizer = slicedims(size(p), NP + 1, N)
    sizes = slicedims(size(p.stackedscales), NP + 1, NP + NS)
    # @show size(p), sizep, sizer, sizes
    inner_multiply2!(x, p.stackedscales, y, sizep, sizer, sizes)

    return x
end

function inner_multiply!(x, s, y, sizep, sizer, sizes)
    for ip in CartesianIndices(sizep), ir in CartesianIndices(sizer)
        x[ip, ir] = 0
        for is in CartesianIndices(sizes)
            # @inline y[ip, ir, is] = x[ip, ir] *s[ip, is]
            @inline x[ip, ir] += y[ip, ir, is] * s[ip, is]
        end
    end
    return x
end



function inner_multiply2!(x, s, y, sizep, sizer, sizes)
    x .= 0
    for ip in CartesianIndices(sizep),
        ir in CartesianIndices(sizer),
        is in CartesianIndices(sizes)
        # @inline y[ip, ir, is] = x[ip, ir] *s[ip, is]
        @inline x[ip, ir] += y[ip, ir, is] * s[ip, is]

    end
    return x
end

struct ScaledCopies{TS,PF,PB,TBF} <:
       AbstractScaledCopiesSet where {TS<:FeasibleSet,PF<:AbstractSCPlan,PB<:AbstractSCPlan, TBF}
    set::TS
    fplan::PF
    bplan::PB
    bufer::TBF
end

function ScaledCopies(A::FeasibleSet, scales)
    a = getelement(A)
    fplan = plan_SC(scales, getelement(A))
    return ScaledCopies(A, fplan, invert(fplan), a)
end

project!(xp, x, feasset::ScaledCopies) = backproject!(xp, x, feasset)
project!(x, feasset::ScaledCopies) = backproject!(x, feasset)

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
function FourierTransformedSet(s::FeasibleSet, dims)
    return FourierTransformedSet(
        s,
        FFTW.plan_ifft(getelement(s), dims),
        FFTW.plan_fft(getelement(s), dims),
        getelement(s),
    )
end

# struct UnitaryTransformedSet{TS,T,N} <: AbstractLinearTransformedSet where {TS <: FeasibleSet,T,N}
#     set::TS
#     fplan::Array{T,N}
#     bplan::Array{T,N}
#     bufer
# end
