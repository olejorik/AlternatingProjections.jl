struct LCSet{T,N,K} <: AmplitudeConstrainedSet where {T<:Real,N,K}
    amp::Array{T,N}
    projdims::NTuple{K,Int}
    addsize::NTuple{K,Int}
end

amp(s::LCSet) = s.amp
projdims(s::LCSet) = s.projdims
addsize(s::LCSet) = s.addsize
addims(S::LCSet{T,N,K}) where {T,N,K} = Tuple(setdiff(1:(N + K), projdims(S)))

function getelement(S::LCSet)
    eldims = collect(size(amp(S)))
    for i in eachindex(projdims(S))
        insert!(eldims, projdims(S)[i], addsize(S)[i])
    end
    ret = zeros(ComplexF64, eldims...)
    for (i, s) in enumerate(eachslice(ret; dims=addims(S)))
        s[1] = amp(S)[i]
    end
    return ret
end

function project!(x, S::LCSet)
    for (i, s) in enumerate(eachslice(x; dims=addims(S)))
        #     @show s
        #     @show i
        #     @show sum(abs2, s)
        if sum(abs2, s) != 0
            s .*= (amp(S)[i] / sqrt(sum(abs2, s)))
        end
    end
    return x
end


export LCSet
