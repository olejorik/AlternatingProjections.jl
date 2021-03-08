"""
# type AmplitudeConstraint


# Examples

```jldoctest
julia>
```
"""
abstract type  AmplitudeConstrainedSet <: FeasibleSet end
export AmplitudeConstrainedSet

"""
    ConstrainedByAmplitude(a)

Set defined by the amplitude constraint `|x| = a`.

"""
struct ConstrainedByAmplitude{T,N} <: AmplitudeConstrainedSet where {T <: Real, N}
    amp::Array{T,N}  #
end

function ConstrainedByAmplitude(a::AbstractArray{Union{T,Nothing},N}) where {T,N} 
    ConstrainedByAmplitude{Union{T,Nothing},N}(a)
end

function ConstrainedByAmplitude(a::AbstractArray{T,N}) where {T,N} 
    ConstrainedByAmplitude{T,N}(a)
end

# function ConstrainedByAmplitude(a::AbstractArray{T,N}) 
#     ConstrainedByAmplitude{T,N}(a)    
# end
export ConstrainedByAmplitude

"""
ConstrainedByAmplitudeMasked(a, mask::Vector{CartesianIndex{N}})

The amplitude constrained set only in the indexes given by mask:  `|xᵢ| = aᵢ` for `i ∈ mask`.
"""
struct ConstrainedByAmplitudeMasked{T,N} <: AmplitudeConstrainedSet where {T <: Real, N}
    amp::Array{T,N}  #
    mask::Vector{CartesianIndex{N}}  #
end

export ConstrainedByAmplitudeMasked

function project!(xp, x, feasset::ConstrainedByAmplitude)
    # @inbounds for i in eachindex(xp)
    #     xp[i] = feasset.amp[i] * exp( im * angle(x[i]))
    # end
    # xp .= feasset.amp .* exp.(1.0im .* angle.(x)) # this is 1000 times slower
    # xp .= _replace_amp.(feasset.amp, x)

    # xp .= feasset.amp .* _unit_amp.(x)
    # xp .= update_amplitude.(feasset.amp,x)

    @inbounds for i in eachindex(xp)
        xp[i] = update_amplitude(feasset.amp[i], x[i])
    end

    return xp
end

function project!(xp, x, feasset::ConstrainedByAmplitudeMasked)
    @inbounds for i in feasset.mask
        xp[i] = update_amplitude(feasset.amp[i], x[i])
    end
    return xp
end

update_amplitude(amp, x) = isnothing(amp) ? x : amp * _unit_amp(x)

update_amplitude!(xp, amp, x) = isnothing(amp) ? x : amp * _unit_amp(x)

@inline function _unit_amp(z)
    abs(z)≈0 ? one(z) : z/abs(z)
end

@inline function _replace_amp(amp,z)
    abs(z)≈0 ? zero(z) : z* (amp/abs(z))
end


size(feasset::ConstrainedByAmplitude) = size(feasset.amp)
