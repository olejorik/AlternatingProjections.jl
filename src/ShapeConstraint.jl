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
Outside the mask function should be larger than the satruation level `|xᵢ| > s ` for `i ∉ mask`
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

upperthreshold(a, th) = a < th ? th : a

function project!(xp, x, feasset::ConstrainedByShapeSaturated)    
    s = abs.(x)[feasset.mask]' * feasset.amp[feasset.mask] / feasset.n
    # s = x[feasset.mask]' * feasset.amp[feasset.mask] / feasset.n
    @inbounds for i in feasset.mask
        xp[i] = update_amplitude(s * feasset.amp[i], x[i])
    end
    @inbounds for i in feasset.sat
        xp[i] = update_amplitude(upperthreshold.(abs(x[i]),s), x[i])
    end
    return xp
end

function project!(x, feasset::ConstrainedByShapeSaturated)    
    s = abs.(x)[feasset.mask]' * feasset.amp[feasset.mask] / feasset.n
    @inbounds for i in feasset.mask
        x[i] = update_amplitude(s * feasset.amp[i], x[i])
    end
    @inbounds for i in feasset.sat
        x[i] = update_amplitude(upperthreshold.(abs(x[i]),s), x[i])
    end
    return x
end