"""
    AmplitudeConstrainedSet 

Abstract set of complex-values `x` with a given absolute value `|x| = amp`. The amplitude can be extracted by function `amp`.

"""
abstract type  AmplitudeConstrainedSet <: FeasibleSet end
export AmplitudeConstrainedSet

amp(s::AmplitudeConstrainedSet) = s.amp #|| error("amplitude is not defined for $(typeof(s))")
getelement(s::AmplitudeConstrainedSet) = complex(float(amp(s))) # element of the set is a complex array

"""
    ConstrainedByAmplitude{T,N}

Set of abstract arryas with element type `T` and dimensions `N` defined by the amplitude constraint `|x| = amp`. 

"""
struct ConstrainedByAmplitude{T,N} <: AmplitudeConstrainedSet where {T <: Real, N}
    amp::Array{T,N}  #
end

amp(s::ConstrainedByAmplitude) = s.amp

"""
ConstrainedByAmplitude(a::AbstractArray{T,N})

Construct set defined by the amplitude constraint `|x| = a`. Type and dimension of the set are inhereited from the array.
"""
function ConstrainedByAmplitude(a::AbstractArray{T,N}) where {T,N} 
    ConstrainedByAmplitude{T,N}(a)
end

function ConstrainedByAmplitude(a::AbstractArray{Union{T,Nothing},N}) where {T,N} 
    ConstrainedByAmplitude{Union{T,Nothing},N}(a)
end

# function ConstrainedByAmplitude(a::AbstractArray{T,N}) 
#     ConstrainedByAmplitude{T,N}(a)    
# end
export ConstrainedByAmplitude


"""
    ACset{T,N,M,K}(amp,projdims) 

Constructs an extended version of the AmplitudeCosntrained set. Here, `amp` can be a tuple
of amplitudes, and `projdims` are the dimeshoins, along with the length of the vector is measured.

TODO extend description
"""
struct ACset{T,N,M,K} <: AmplitudeConstrainedSet where {T <: Real, N,M,K}
    amp::NTuple{M,Array{T,N}} # tuple of arrays; TODO expand to more dimesnions via array? Array{Array{T,N},M} #
    projdims::NTuple{K, Int} # dimensions along which the projection is made
end

amp(s::ACset) = s.amp
projdims(s::ACset) = s.projdims

function ACset(amp::NTuple{M, AbstractArray{T,N}})  where {M,T,N}
    ACset{T,N,M,0}(amp, ())
end

ACset(amp::AbstractArray{T,N}, dims...)  where {T,N,K} =   ACset((amp,), dims...)

export ACset

"""
    ConstrainedByAmplitudeMasked(a, mask::Vector{CartesianIndex{N}})
    ConstrainedByAmplitudeMasked(a, AbstractArray{Bool})

The amplitude constrained set only in the indexes given by mask:  `|xᵢ| = aᵢ` for `i ∈ mask`.
"""
struct ConstrainedByAmplitudeMasked{T,N} <: AmplitudeConstrainedSet where {T <: Real, N}
    amp::Array{T,N}  #
    mask::Vector{CartesianIndex{N}}  #
end

ConstrainedByAmplitudeMasked(a, m::AbstractArray{Bool}) = ConstrainedByAmplitudeMasked(a, findall(m))

export ConstrainedByAmplitudeMasked
amp(s::ConstrainedByAmplitudeMasked) = s.amp


# here the backward plan is in place
# FourierTransformedSet(s::AmplitudeConstrainedSet) = 
#     FourierTransformedSet(s, FFTW.plan_fft(getelement(s)), FFTW.plan_ifft!(getelement(s)))

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

function project!(x, feasset::ConstrainedByAmplitude)
    for i in eachindex(x)
        x[i] = update_amplitude(feasset.amp[i], x[i])
    end

    return x
end

function project!(xp, x, feasset::ConstrainedByAmplitudeMasked)
    @inbounds for i in feasset.mask
        xp[i] = update_amplitude(feasset.amp[i], x[i])
    end
    return xp
end

function project!(x, feasset::ConstrainedByAmplitudeMasked)
    @inbounds for i in feasset.mask
        x[i] = update_amplitude(feasset.amp[i], x[i])
    end
    return x
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


include("ShapeConstraint.jl")