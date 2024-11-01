"""
    LCSet{T, N, K}(amp, projdims, addsize)

Construct length-constrained set such that length of the slices along the axes given by `projdims` is equal to `amp`. The size of the set element is given by `addsize` along the projected dimensions `projdims` and by the size of `amp` along the other dimensions.

Parameters:
    amp::Array{T,N}
    projdims::NTuple{K,Int}
    addsize::NTuple{K,Int}

Example
=======
```jldoctes
julia> amp = [1 2; 3 4]; projdims = (2,3);addsize= (5,4); A = LCSet(amp, projdims, addsize)
LCSet{Int64, 2, 2}([1 2; 3 4], (2, 3), (5, 4), (1, 4))

julia> z = getelement(A)
2×5×4×2 Array{ComplexF64, 4}:
[:, :, 1, 1] =
 1.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 3.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im

[:, :, 2, 1] =
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im

[:, :, 3, 1] =
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im

[:, :, 4, 1] =
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im

[:, :, 1, 2] =
 2.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 4.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im

[:, :, 2, 2] =
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im

[:, :, 3, 2] =
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im

[:, :, 4, 2] =
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im

julia> vec(Int.(sqrt.(sum(abs2, z, dims = projdims)))) == vec(amp)
true
```


"""
struct LCSet{T,N,K,M} <: AmplitudeConstrainedSet where {T<:Real,N,K,M}
    amp::Array{T,N}
    projdims::NTuple{K,Int}
    addsize::NTuple{K,Int}
    adddims::NTuple{N,Int}
    buffer::Array{Float64,M} # Todo use the same function as in getelement
end

LCSet(
    amp::Array{T,N}, projdims::NTuple{K,Int}, addsize::NTuple{K,Int}
) where {T<:Real,N,K} = LCSet(
    amp,
    projdims,
    addsize,
    sorted_setdiff(ntuple(i -> i, N + K), projdims),
    similar(amp, Float64, insert_ones(size(amp), projdims)),
)

# This was slow because of type instability of eachslice
# struct LCSet{T,N,K} <: AmplitudeConstrainedSet where {T<:Real,N,K}
#     amp::Array{T,N}
#     projdims::NTuple{K,Int}
#     addsize::NTuple{K,Int}
#     adddims::NTuple{N, Int}
# end

# LCSet(amp::Array{T,N}, projdims::NTuple{K,Int}, addsize::NTuple{K,Int}) where {T<:Real,N,K} = LCSet(amp, projdims, addsize, sorted_setdiff(ntuple(i->i,N+K), projdims))
# 


amp(s::LCSet) = s.amp
projdims(s::LCSet) = s.projdims
addsize(s::LCSet) = s.addsize
# adddims(S::LCSet{T,N,K}) where {T,N,K} = NTuple{N, Int}(setdiff(1:(N + K), projdims(S)))
# adddims(S::LCSet{T,N,K}) where {T,N,K} = NTuple{N, Int}(sorted_setdiff(ntuple(i->i,N+K), projdims(S))) # moved to the set structure
adddims(s::LCSet) = s.adddims


function getelement(S::LCSet)
    eldims = collect(size(amp(S)))
    for i in eachindex(projdims(S))
        insert!(eldims, projdims(S)[i], addsize(S)[i])
    end
    ret = zeros(ComplexF64, eldims...)
    for (i, s) in enumerate(eachslice(ret; dims=adddims(S)))
        s[1] = amp(S)[i]
    end
    return ret
end


# This was for the first impementation of LCSet
# function project!(x::Array{T,N}, S::LCSet) where{T,N}
#     update_amp_slice!(x, amp(S), adddims(S))
# end
# 


function project!(x::Array{T,N}, S::LCSet) where {T,N}
    # S.buffer .= sum!(abs2, x, dims = S.projdims)
    sum!(abs2, S.buffer, x)
    return x .= x ./ S.buffer
end


function update_amp_slice!(x::Array{T,M}, amp, dims::NTuple{N,Int}) where {T,M,N}
    for (i, s) in enumerate(eachslice(x; dims=dims))
        sss = sum(abs2, s)
        #     @show s
        #     @show i
        #     @show sum(abs2, s)
        if sss != 0
            s .*= (amp[i] / sqrt(sss))
        end
    end
    return x
end

function update_amp_slice2!(x, amp, dims)
    for (i, s) in enumerate(eachslice(x; dims=dims))
        # sss = 0
        # for j in eachindex(s)
        #     sss += s[j]
        # end
        # if sss != 0
        #     s .*= (amp[i] / sqrt(sss))
        # end
    end
    return x
end

viewdims(M, dims) = ntuple(i -> i ∈ dims ? i : Colon(), M)

function update_amp_slice3!(x, amp, viewdims)
    for (i, s) in enumerate(eachslice(x; dims=dims))
        sss = sum(abs2, s)
        #     @show s
        #     @show i
        #     @show sum(abs2, s)
        if sss != 0
            s .*= (amp[i] / sqrt(sss))
        end
    end
    return x
end

export LCSet

