
import LinearAlgebra: mul!

struct PhaseDiversityPlan{T,N,M} <: AbstractSCPlan{T,N,M}
    diversity::Array{Array{T,N},M}
end
struct InversePhaseDiversityPlan{T,N,M} <: AbstractSCPlan{T,N,M}
    diversity::Array{Array{T,N},M}
end

diversities(p::PhaseDiversityPlan) = p.diversity

function mul!(y, p::PhaseDiversityPlan, x)
    for i in CartesianIndices(y)
        y[i] .= p.diversity[i] .* x
    end
    return y
end

struct PhaseDiversedSet{TS,PF,PB} <: AbstractScaledCopiesSet
    set::TS
    fplan::PF
    bplan::PB
    bufer
end

function PhaseDiversedSet(
    s::FeasibleSet, phases::Array{T}
) where {T<:Array{N} where {N<:Real}}
    scales = map(x -> exp.(im .* x), phases)
    return PhaseDiversedSet(s, plan_SC(scales), invert(plan_SC(scales)), getelement(s))
end

project!(xp, x, feasset::PhaseDiversedSet) = backproject!(xp, x, feasset)
project!(x, feasset::PhaseDiversedSet) = backproject!(x, feasset)

# """
#     Set created by several point-wise scaled copies of each element of the set.
#     `{x} → {(a₁x,…,aₙx)}`.

# """
# struct ScaledCopies{TS,T,N} <: LinearTransformedSet
#     set::TS
#     diversities::Array{T,N}

# end
