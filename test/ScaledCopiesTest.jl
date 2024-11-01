# Test of allocations in scaled copies ditect and inverse sets
#

using AlternatingProjections, FFTW
using Profile, BenchmarkTools
using LinearAlgebra: mul!
using Test
#
# We try to impement the general formule that  if q is obtained from p using scales sc and (first) dimensions of p
# q[i_p, i_r, i_sc] = p[i_p, i_r] ⋅ sc[i_sc][i_p]
#
p2 = rand(3, 4, 2) # the dims for multiplication is the first 1
sc2 = reshape([collect(i:(i + 2)) for i in 1.0:3:10], (2, 2)) # 2 x 2 array of 3-vectors
plan = AlternatingProjections.plan_SC(sc2, p2, [1])
q = AlternatingProjections.getplanelement(plan)
mul!(q, plan, p2)

iplan = AlternatingProjections.invert(plan)
x = zero(p2)
mul!(x, iplan, q) ≈ p2


function fplan_test(n, q=q, plan=plan, p2=p2)
    for i in 1:n
        for ind in eachindex(p2)
            p2[ind] = rand()
        end
        mul!(q, plan, p2)
    end
end

@profview fplan_test(1)

@profview fplan_test(100000)

@btime fplan_test(10, $q, $plan, $p2)



function iplan_test(n, x=x, iplan=iplan, q=q)
    for i in 1:n
        for ind in eachindex(q)
            q[ind] = rand()
        end
        mul!(x, iplan, q)
    end
end

@profview iplan_test(10000)
@btime iplan_test(10, $x, $iplan, $q)

#it works, now the sets
A = ConstrainedByAmplitude([1, 2.0, 0])
B = AlternatingProjections.ScaledCopies(A, sc2)
z = AlternatingProjections.getelement(B)
@test z[:, :, 1] == [
    1.0+0.0im 4.0+0.0im
    4.0+0.0im 10.0+0.0im
    0.0+0.0im 0.0+0.0im
]
@test z[:, :, 2] == [
    7.0+0.0im 10.0+0.0im
    16.0+0.0im 22.0+0.0im
    0.0+0.0im 0.0+0.0im
]

x = rand(ComplexF64, size(z)...)
AlternatingProjections.backproject!(x, B)

@btime project!($x, $B)
@profview_allocs for i in 1:10000
    project!(x, B)
end
@bprofile for i in 1:1000
    project!($x, $B)
end

@btime project!($x, $B)

# check projections
A = ConstrainedByAmplitude([1 2.0; 0 3])
sc = [[1 0; 0 1], [1 1; 0 0.0]]
B = AlternatingProjections.ScaledCopies(A, sc)
z = AlternatingProjections.getelement(B)
y = project(-ones(size(z)), B)

x = rand(ComplexF64, size(z))
@bprofile for i in 1:10
    for i in eachindex($x)
        $x[i] = rand()
    end
    project!($x, $B)
end


diversities = [cispi.([0.25 0.25; 0.25 0.25]), cispi.([0.5 0.5; -0.5 0.5])]
C = AlternatingProjections.ScaledCopies(B, diversities, [1, 2])
z2 = AlternatingProjections.getelement(C)
y2 = project(-ones(ComplexF64, size(z2)), C)

abs.(y2)
angle.(y2)

x = rand(ComplexF64, size(z2))
@bprofile for i in 1:10
    for i in eachindex($x)
        $x[i] = rand()
    end
    project!($x, $C)
end
# good, zero allocations!
# 

