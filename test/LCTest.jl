# Test of allocations in scaled copies ditect and inverse sets
#

using AlternatingProjections, FFTW
using Profile, BenchmarkTools
using LinearAlgebra: mul!
using Test
#

amp = [1 2; 3 4];
projdims = (2, 3);
addsize = (5, 4);

A = LCSet(amp, projdims, addsize)

z = getelement(A)
@test vec(Int.(sqrt.(sum(abs2, z; dims=projdims)))) == vec(amp)


x = rand(ComplexF64, size(z))
@btime project!($x, $A)


@bprofile for i in 1:10
    for i in eachindex($x)
        $x[i] = rand()
    end
    project!($x, $A)
end

# Relative to the type stability
arr = reshape(collect(1:27), 3, 3, 3);

function f3(arr, d)
    for (i, x) in enumerate(eachslice(arr; dims=d))
        @show sum(abs2, x)
    end
end

function f4(arr)
    for (i, x) in enumerate(eachslice(arr; dims=(1, 3)))
        @show sum(abs2, x)
    end
end

function f5(arr, d)
    for (i, x) in enumerate(@inline eachslice(arr; dims=d))
        @show sum(abs2, x)
    end
end

@code_warntype f3(arr, (1, 3))

@code_warntype f4(arr)

@code_warntype f5(arr, (1, 3))

## use buffer?
buf .= sum(abs2, x; dims=projdims)
x .= x ./ buf

@bprofile for i in 1:10
    for i in eachindex($x)
        $x[i] = rand()
    end
    # project!($x, $A)
    $buf .= sum(abs2, $x; dims=$projdims)
    $x .= $x ./ $buf

end


