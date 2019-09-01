"""
# type AmplitudeConstraint

- Julia version: 1.1.0
- Author: Oleg Soloviev
- Date: 2019-09-01

# Examples

```jldoctest
julia>
```
"""
struct AmplitudeConstraint <: FeasibleSet
    amp::Array{T} where T <: Real #todo nongegative
end
export AmplitudeConstraint

function project(x, feasset::AmplitudeConstraint)
    return feasset.amp .* exp.( im * angle.(x))
end

size(feasset::AmplitudeConstraint) = size(feasset.amp)
