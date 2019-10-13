"""
# type AmplitudeConstraint


# Examples

```jldoctest
julia>
```
"""
abstract type  AmplitudeConstrainedSet <: FeasibleSet end
export AmplitudeConstrainedSet

struct ConstrainedByAmplitude <: AmplitudeConstrainedSet
    amp::Array{T} where T <: Real #todo nongegative
end
export ConstrainedByAmplitude

function project(x, feasset::ConstrainedByAmplitude)
    return feasset.amp .* exp.( im * angle.(x))
end

size(feasset::ConstrainedByAmplitude) = size(feasset.amp)
