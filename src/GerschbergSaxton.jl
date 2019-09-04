"""
    GS

Type representing Gerchberg-Saxton phase retrieval algorithm.

Requires (in the current realisation) two amplitude constraints or two arrays of equal sizes.
"""
struct GS <: APMethod #todo should be sets part of this or added to the step! only?
    a::AmplitudeConstraint
    A::AmplitudeConstraint
    forward
    backward
end
export GS

function GS(a::Array,A::Array)
    size(a) == size(A) ? (P, PB) = (plan_fft(A), plan_ifft(A)) : error("Array sizes do not match")
    GS(AmplitudeConstraint(a), AmplitudeConstraint(A), x -> P * x, x -> PB * x)
end

function init!(alg::GS,  x⁰)
    a = alg.a
    A = alg.A
    if !( size(a) == size(A) == size(x⁰) )
        println("cannot intialise Gerschberg-Saxton, dimensions do not match")
        # raise error or pad a smaller array with zeroes
        # also can include realingment of the sets to remove any linear pahse  term. But then need to make GS mutable
        return
    end
end

function step!(alg::GS, xᵏ)
    a = alg.a
    A = alg.A

#     # step by step explanation
#     yk = fft(xᵏ)
#     println(yk)
#     zk = project(A, yk)
#     println('z',zk)
#     yk = ifft(zk)
#     println(yk)
#     zk = project(a, yk)
#     println('z',zk)
#

    return project(a, ifft(project(A, fft(xᵏ))))
    #todo replace with plan
end
