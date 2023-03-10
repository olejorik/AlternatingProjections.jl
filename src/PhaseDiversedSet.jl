
import LinearAlgebra: mul!

struct PhaseDiversityPlan{T, N} <: AbstractFFTs.Plan{T}
    diversity::Array{T,N}
end

diversities(p::PhaseDiversityPlan) = p.diversity

function mul!(y,p::PhaseDiversityPlan, x)
    for i in CartesianIndices(y)
        y[i] = p.diversity[i] .* x
    end
    return y
end

struct PhaseDiversedSet{TS, TP} <: AbstractLinearTransformedSet where {
    TS<:FeasibleSet , TP <: PhaseDiversityPlan
}
    set:: TS
    fplan::TP
    bplan::TP
    bufer
end

function PhaseDiversedSet(s::FeasibleSet, phases::Array{Real})
end


# """
#     Set created by several point-wise scaled copies of each element of the set.
#     `{x} → {(a₁x,…,aₙx)}`.


# """
# struct ScaledCopies{TS,T,N} <: LinearTransformedSet
#     set::TS
#     diversities::Array{T,N}

# end


