# When scaling of the amplitude is allowed, we call it's shape. 
# TODO Technically, it should be supertype or subtype of amplitude constrain? Just now I do it for easyness as subtype. 

"""
ConstrainedByShape{T,N}

Set of abstract arryas with elemet typ `T` and dimensions `N defined by the shape constraint `|x| = s amp` for some `s`. 
Field n contains the sum of square of all elements.
"""
struct ConstrainedByShape{T,N} <: AmplitudeConstrainedSet where {T <: Real, N}
    amp::Array{T,N}  #
    n::T
end

"""
ConstrainedByShape(a::AbstractArray{T,N})

Construct set defined by the shape constraint `|x| = s a` for some scaling `s`. Type and dimension of the set are inhereited from the array.
"""
function ConstrainedByShape(a::AbstractArray{T,N}) where {T,N} 
    ConstrainedByShape{T,N}(a, sum(a .^2))
end

# function ConstrainedByShape(a::AbstractArray{Union{T,Nothing},N}) where {T,N} 
#     ConstrainedByShape{Union{T,Nothing},N}(a)
# end

# function ConstrainedByAmplitude(a::AbstractArray{T,N}) 
#     ConstrainedByAmplitude{T,N}(a)    
# end
export ConstrainedByShape

"""
ConstrainedByAmplitudeMasked(a, mask::Vector)
ConstrainedByAmplitudeMasked(a, AbstractArray{Bool})

The amplitude constrained set only in the indexes given by mask:  `|xᵢ| = aᵢ` for `i ∈ mask`.
"""
struct ConstrainedByShapeMasked{T,N} <: AmplitudeConstrainedSet where {T <: Real, N}
    amp::Array{T,N}  #
    mask::Union{Vector{CartesianIndex{N}}, Vector{Int}} #because 1D arrays are indexed as Vector
    n::T
end


ConstrainedByShapeMasked(amp, mask::Vector) = ConstrainedByShapeMasked(amp, mask, sum(abs2,amp[mask]))

ConstrainedByShapeMasked(a, m::AbstractArray{Bool}) = ConstrainedByShapeMasked(a, findall(m))

export ConstrainedByShapeMasked


"""
ConstrainedByShapeSaturated(a, mask::Vector)
ConstrainedByShapeSaturated(a, AbstractArray{Bool})

The amplitude constrained set only in the indexes given by mask:  `|xᵢ| = s aᵢ` for `i ∈ mask` for some `s`.
Outside the mask function should be larger than the satruation level `|xᵢ| > s ` for `i ∉ mask`.

For this set, the amplitude should be provided in the range from 0 to 1 (1 corresponding to the saturated values).
"""
struct ConstrainedByShapeSaturated{T,N} <: AmplitudeConstrainedSet where {T <: Real, N}
    amp::Array{T,N}  #
    mask::Union{Vector{CartesianIndex{N}}, Vector{Int}} #because 1D arrays are indexed as Vector
    sat::Union{Vector{CartesianIndex{N}}, Vector{Int}} # complementary set to mask
    n::T
end


ConstrainedByShapeSaturated(amp, mask::Vector, sat::Vector) = ConstrainedByShapeSaturated(amp, mask,sat, sum(abs2,amp[mask]))

ConstrainedByShapeSaturated(a, m::AbstractArray{Bool}) = ConstrainedByShapeSaturated(a, findall(m), findall(.!(m)))

export ConstrainedByShapeSaturated

struct ConstrainedByShapeClipped{T,N} <: AmplitudeConstrainedSet where {T <: Real, N}
    amp::Array{T,N}  #
    mid::Union{Vector{CartesianIndex{N}}, Vector{Int}} #because 1D arrays are indexed as Vector
    high::Union{Vector{CartesianIndex{N}}, Vector{Int}} # complementary set to mask
    low::Union{Vector{CartesianIndex{N}}, Vector{Int}} # complementary set to mask
    vhigh::T
    vlow::T
    n::T
end


ConstrainedByShapeClipped(amp, mid, high, low, vhigh, vlow) = ConstrainedByShapeClipped(amp, mid, high, low, vhigh, vlow, sum(abs2,amp[mid]))

ConstrainedByShapeClipped(amp, vhigh, vlow) = ConstrainedByShapeClipped(amp, findall( amp .>= vlow .&& amp .<= vhigh),findall( amp .>= vhigh),findall( amp .<= vlow ), vhigh, vlow)



export ConstrainedByShapeClipped



function project!(xp, x, feasset::ConstrainedByShape)
    s = abs.(x)[:]' * feasset.amp[:] / feasset.n
    @inbounds for i in eachindex(xp)
        xp[i] = update_amplitude(s * feasset.amp[i], x[i])
        end
    return xp
end


function project!(x, feasset::ConstrainedByShape)
    s = abs.(x)[:]' * feasset.amp[:] / feasset.n
    for i in eachindex(x)
        x[i] = update_amplitude(s * feasset.amp[i], x[i])
    end

    return x
end

function project!(xp, x, feasset::ConstrainedByShapeMasked)    
    s = abs.(x)[feasset.mask]' * feasset.amp[feasset.mask] / feasset.n
    # s = x[feasset.mask]' * feasset.amp[feasset.mask] / feasset.n
    @inbounds for i in feasset.mask
        xp[i] = update_amplitude(s * feasset.amp[i], x[i])
    end
    return xp
end

function project!(x, feasset::ConstrainedByShapeMasked)    
    s = abs.(x)[feasset.mask]' * feasset.amp[feasset.mask] / feasset.n
    @inbounds for i in feasset.mask
        x[i] = update_amplitude(s * feasset.amp[i], x[i])
    end
    return x
end

using Roots

function project!(xp, x, feasset::ConstrainedByShapeSaturated)    
    s0 = abs.(x)[feasset.mask]' * feasset.amp[feasset.mask] / feasset.n

    xs = sort(abs.(x[feasset.sat]))
    function rp(s)
        j = searchsortedlast(xs,s)
        return (s - s0) * feasset.n + j*s - sum(xs[1:j])
    end
    # println(s0, extrema(xs)) # debug
    smin, smax = extrema(xs)
    sopt = find_zero(rp,(min(smin,s0), max(smax,s0)))

    @inbounds for i in feasset.mask
        xp[i] = update_amplitude(sopt * feasset.amp[i], x[i])
    end
    @inbounds for i in feasset.sat
        xp[i] = update_amplitude(max.(abs(x[i]),sopt), x[i])
    end

    return xp
end

project!(x, feasset::ConstrainedByShapeSaturated) = project!(x, x, feasset)

function _project!(xp, x, feasset::ConstrainedByShapeClipped)    #introduced helper function to get access to sopt
    s0 = abs.(x)[feasset.mid]' * feasset.amp[feasset.mid] / feasset.n
    b=feasset.vhigh
    a=feasset.vlow

    xhigh = sort(abs.(x[feasset.high]))
    xlow = sort(abs.(x[feasset.low]), rev = true)
    function rp(s)
        jhigh = searchsortedlast(xhigh, s*b)
        jlow = searchsortedlast(xlow,s*a)
        return (s - s0) * feasset.n + jhigh*s*b - sum(xhigh[1:jhigh]) - jlow*s*a + sum(xlow[1:jlow])  
    end
    # println(s0, extrema(xs)) # debug
    sminhigh,smaxhigh = extrema(xhigh)
    sminlow,smaxlow = extrema(xlow)
    smax = maximum((s0, smaxhigh/b, smaxlow/a))
    smin = minimum((s0, sminhigh/b, sminlow/a))
    # smax = maximum((s0, smaxhigh/b))
    # smin = minimum((s0, sminlow/a))
    sopt = find_zero(rp, (smin, smax))

    #debug
    # sopt = try
    #     find_zero(rp, (smin, smax))
    #     # sopt = find_zero(rp, (0, smax))
    #     # sopt = find_zero(rp, s0)
    # catch e
    #     @info "s0 = $s0, smin = $smin smax = $smax"
    #     @info "rp(s0) = $(rp(s0)), rp(smin) = $(rp(smin)) rp(smax) = $(rp(smax))"
    #     s0
    # end


    @inbounds for i in feasset.mid
        xp[i] = update_amplitude(sopt * feasset.amp[i], x[i])
    end
    @inbounds for i in feasset.low
        xp[i] = update_amplitude(min.(abs(x[i]),sopt * a), x[i])
    end
    @inbounds for i in feasset.high
        xp[i] = update_amplitude(max.(abs(x[i]),sopt * b), x[i])
    end

    return xp, sopt
end

project!(xp, x, feasset::ConstrainedByShapeClipped) = _project!(xp, x, feasset)[1]

project!(x, feasset::ConstrainedByShapeClipped) = project!(x, x, feasset)
