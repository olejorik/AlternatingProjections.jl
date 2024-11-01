# Type stable set differerence
# 
# from https://discourse.julialang.org/t/type-stable-difference-of-tuples/3933/4
import Base: tail
@inline function sorted_setdiff(t1::Tuple, t2::Tuple)
    if t1[1] == t2[1]
        sorted_setdiff(tail(t1), tail(t2))
    else
        (t1[1], sorted_setdiff(tail(t1), t2)...)
    end
end
@noinline sorted_setdiff(t1::Tuple{}, t2::Tuple) = error("did not find $(t2[1])")
sorted_setdiff(t1::Tuple, ::Tuple{}) = t1
sorted_setdiff(::Tuple{}, ::Tuple{}) = ()


# type stable insertion of ones but with allocations
function insert_ones(ordims::NTuple{K,Int}, indims::NTuple{M,Int}) where {K,M}
    arr = collect(ordims)
    for i in eachindex(indims)
        insert!(arr, indims[i], 1)
    end
    return NTuple{K + M,Int}(arr)
end

# function to update amplitude of a complex number
# x-> amp * x/|x|
update_amplitude(amp, x) = isnothing(amp) ? x : amp * _unit_amp(x)

update_amplitude!(xp, amp, x) = isnothing(amp) ? x : amp * _unit_amp(x)

@inline function _unit_amp(z)
    return abs(z) ≈ 0 ? one(z) : z / abs(z)
end

@inline function _replace_amp(amp, z)
    return abs(z) ≈ 0 ? zero(z) : z * (amp / abs(z))
end

