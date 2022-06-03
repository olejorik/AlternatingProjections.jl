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

#= 
A simple (and in the same time the most general) iterator based on repetion of some function or operator
=#

"""
    Opiterator1(x⁰, f)  

Create a simple iteration with operator `f` and inital value `x⁰`.

# Examples

```jldoctest
julia> ttt = Opiterator([3], x -> x .+ 2)
AlternatingProjections.Opiterator1{Vector{Int64}}([3], var"#11#12"())

julia> for state in Base.Iterators.take(ttt,5)
       println(state.xᵏ[1])
       end
5
7
9
11
13

julia> ttt20 = Iterators.take(ttt,20);

julia> for state in ttt20
       print(state.xᵏ[1])
       end
5791113151719212325272931333537394143
```

## Example with Images.jl
```
julia> using Images

julia> img = Gray.(rand(Float64, (256,256)))

julia> blurop = Opiterator(img, x -> imfilter(x, Kernel.gaussian(2)));

julia> blurop20 = Iterators.take(blurop, 20);

julia> mosaic([copy(b.xᵏ⁻¹) for b in blurop20]; nrow=4, rowmajor = true)
```


"""
struct Opiterator1{TX}  # where TX
    x⁰::TX
    f::Function 
end

Opiterator = Opiterator1


mutable struct Opstate1{TX}  # where TX
    xᵏ⁻¹::TX
    xᵏ::TX
end

Opstate = Opstate1

function iterate(iter::Opiterator) 
    # state = Opstate(copy(iter.x⁰), project(iter.x⁰,iter.setB), project(project(iter.x⁰,iter.setB),iter.setA))
    if iter.x⁰ === nothing error("intial value is required for $iter")
    else 
        x0 = copy(iter.x⁰)
        x1 = copy(iter.x⁰)
    end
    x1 .= iter.f(x0)
    state = Opstate(x0, x1)
    return state, state 
end

function iterate(iter::Opiterator, state::Opstate)
    state.xᵏ⁻¹ .= state.xᵏ
    state.xᵏ .= iter.f(state.xᵏ⁻¹)
    return state, state
end

Iterators.IteratorSize(::Opiterator) = Iterators.IsInfinite()

export Opstate, Opiterator

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

