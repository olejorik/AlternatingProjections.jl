import Base: iterate

struct APiterate{TA, TB, TX} # where {TA <:FeasibleSet, TB <:FeasibleSet, TX}
    setA::TA
    setB::TB
    x⁰::AbstractArray{TX}
end

mutable struct APstate{TX}
    xᵏ⁻¹::AbstractArray{TX}
    yᵏ::AbstractArray{TX}
    xᵏ::AbstractArray{TX}
end

function iterate(iter::APiterate) 
    # state = APstate(copy(iter.x⁰), project(iter.x⁰,iter.setB), project(project(iter.x⁰,iter.setB),iter.setA))
    state = APstate(copy(iter.x⁰), copy(iter.x⁰), copy(iter.x⁰))
    return state, state 
end

function iterate(iter::APiterate, state::APstate)
    state.xᵏ⁻¹ .= state.xᵏ
    project!(state.yᵏ, state.xᵏ⁻¹, iter.setB)
    project!(state.xᵏ, state.yᵏ, iter.setB)
    return state, state
end