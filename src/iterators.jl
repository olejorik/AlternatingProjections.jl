import Base: iterate

struct APiterator{TA, TB, TX} # where {TA <:FeasibleSet, TB <:FeasibleSet, TX}
    setA::TA
    setB::TB
    x⁰::TX
end



mutable struct APstate{TX}
    xᵏ⁻¹::TX
    yᵏ::TX
    xᵏ::TX
end

function iterate(iter::APiterator) 
    # state = APstate(copy(iter.x⁰), project(iter.x⁰,iter.setB), project(project(iter.x⁰,iter.setB),iter.setA))
    if iter.x⁰ === nothing x0 = complex(one(iter.setA.amp)) 
    else x0 = copy(iter.x⁰)
    end
    y1 = project(x0, iter.setB)
    x1 = project(y1, iter.setA)
    state = APstate(x0, y1, x1)
    return state, state 
end

function iterate(iter::APiterator, state::APstate)
    state.xᵏ⁻¹ .= state.xᵏ
    project!(state.yᵏ, state.xᵏ⁻¹, iter.setB)
    project!(state.xᵏ, state.yᵏ, iter.setA)
    return state, state
end

Iterators.IteratorSize(::APiterator) = Iterators.IsInfinite()

# wrappers (see https://lostella.github.io/2018/07/25/iterative-methods-done-right.html)

struct HaltingIterable{I, F}
    iter::I
    fun::F
end

function iterate(iter::HaltingIterable)
    next = iterate(iter.iter)
    return dispatch(iter, next)
end

function iterate(iter::HaltingIterable, (instruction, state))
    if instruction == :halt return nothing end
    next = iterate(iter.iter, state)
    return dispatch(iter, next)
end

function dispatch(iter::HaltingIterable, next)
    if next === nothing return nothing end
    return next[1], (iter.fun(next[1]) ? :halt : :continue, next[2])
end

halt(iter::I, fun::F) where {I, F} = HaltingIterable{I, F}(iter, fun)

struct TeeIterable{I, F}
    iter::I
    fun::F
end

function iterate(iter::TeeIterable, args...)
    next = iterate(iter.iter, args...)
    if next !== nothing iter.fun(next[1]) end
    return next
end

tee(iter::I, fun::F) where {I, F} = TeeIterable{I, F}(iter, fun)

struct SamplingIterable{I}
    iter::I
    period::UInt
end

function iterate(iter::SamplingIterable, state=iter.iter)
    current = iterate(state)
    if current === nothing return nothing end
    for i = 1:iter.period-1
        next = iterate(state, current[2])
        if next === nothing return current[1], Iterators.rest(state, current[2]) end
        current = next
    end
    return current[1], Iterators.rest(state, current[2])
end

sample(iter::I, period) where I = SamplingIterable{I}(iter, period)

struct StopwatchIterable{I}
    iter::I
end

function iterate(iter::StopwatchIterable)
    t0 = time_ns()
    next = iterate(iter.iter)
    return dispatch(iter, t0, next)
end

function iterate(iter::StopwatchIterable, (t0, state))
    next = iterate(iter.iter, state)
    return dispatch(iter, t0, next)
end

function dispatch(iter::StopwatchIterable, t0, next)
    if next === nothing return nothing end
    return (time_ns()-t0, next[1]), (t0, next[2])
end

stopwatch(iter::I) where I = StopwatchIterable{I}(iter)

function loop(iter)
    x = nothing
    for y in iter x = y end
    return x
end

