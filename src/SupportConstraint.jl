abstract type SupportConstrained <: ConvexSet end

amp(s::SupportConstrained) = s.support #|| error("amplitude is not defined for $(typeof(s))")
getelement(s::SupportConstrained) = complex(float(amp(s))) # element of the set is a complex array

"""
    ConstrainedBySupport(support)

Special type of convex set.

For continuous case: consists of all functions that equals zero outside some fixed area called support:
`` 𝓐_S = \\{f : f(x) = 0, x ∉ S \\} .``

For discrete case: all arrays that equals zero for indexes outside some index set:
`` 𝓐_S = \\{x : x[i] = 0, i ∉ S \\} .``

Currently supports only discrete case, with the support defined as a boolean array.

# Examples

```jldoctest
julia> S =  ConstrainedBySupport([true, false,true])
ConstrainedBySupport(Bool[1, 0, 1])

julia> x = [1, 2, 3]; project(x, S)
3-element Vector{Int64}:
 1
 0
 3

julia> S = ConstrainedBySupport([x^2 + y^2 <=1  for x in -2:.5:2, y in -2:.5:2]);

julia> x = ones(size(S.support));

julia> project(x,S)
9×9 Matrix{Float64}:
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  1.0  1.0  1.0  0.0  0.0  0.0
 0.0  0.0  1.0  1.0  1.0  1.0  1.0  0.0  0.0
 0.0  0.0  0.0  1.0  1.0  1.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
```
"""
struct ConstrainedBySupport <: SupportConstrained
    support::Array{Bool}
end
export ConstrainedBySupport

# function project(x, feasset::ConstrainedBySupport)
#     return feasset.support .* x
# end

function project!(xp, x, feasset::ConstrainedBySupport)
    @inbounds for i in eachindex(xp)
        xp[i] = feasset.support[i] * x[i]
    end

    return xp
end

function project!(x, feasset::ConstrainedBySupport)
    @inbounds for i in eachindex(x)
        x[i] = feasset.support[i] * x[i]
    end

    return x
end

"""
    ConstrainedBySupportNormed(support, norm)

`ConstrainedBySupportNormed(support, norm)` is a set represented by intersection of the 
`ConstrainedBySupport(support)` set and a hypersphere of radius `norm`.

# Examples

```jldoctest
julia>  a = rand((1,10),(5,5));

julia> mask = a .> 5;

julia> aset= ConstrainedBySupportNormed(mask, 10);

julia> b = project(Float64.(a), aset);

julia> sum(abs2,b) ≈ 10^2
true
```

"""
struct ConstrainedBySupportNormed <: SupportConstrained
    support::Array{Bool}
    n::Float64
end

function ConstrainedBySupportNormed(support::Array{T} where {T<:AbstractFloat})
    return ConstrainedBySupportNormed(support .> 0, sqrt(sum(abs2, support)))
end
export ConstrainedBySupportNormed

# function project(x, feasset::ConstrainedBySupport)
#     return feasset.support .* x
# end

function project!(xp, x, feasset::ConstrainedBySupportNormed)
    @inbounds for i in eachindex(xp)
        xp[i] = feasset.support[i] * x[i]
    end

    normratio = feasset.n / sqrt(sum(abs2, xp))
    xp .*= normratio
    return xp
end