# When scaling of the amplitude is allowed, we call it's shape. 
# TODO Techincally, it should be supertype or subtype of amplitude constrain? Just now I do it for easyness as subtype. 

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
ConstrainedByAmplitude(a::AbstractArray{T,N})

Construct set defined by the amplitude constraint `|x| = a`. Type and dimension of the set are inhereited from the array.
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
ConstrainedByAmplitudeMasked(a, mask::Vector{CartesianIndex{N}})
ConstrainedByAmplitudeMasked(a, AbstractArray{Bool})

The amplitude constrained set only in the indexes given by mask:  `|xᵢ| = aᵢ` for `i ∈ mask`.
"""
struct ConstrainedByShapeMasked{T,N} <: AmplitudeConstrainedSet where {T <: Real, N}
    amp::Array{T,N}  #
    mask::Vector{CartesianIndex{N}}  #
    n::T
end


ConstrainedByShapeMasked(amp, mask) = ConstrainedByShapeMasked(amp, mask, sum(amp[mask] .^2))

ConstrainedByShapeMasked(a, m::AbstractArray{Bool}) = ConstrainedByShapeMasked(a, findall(m))

export ConstrainedByShapeMasked



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