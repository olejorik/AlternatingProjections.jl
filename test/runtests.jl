using AlternatingProjections, FFTW, LinearAlgebra
using Test

@testset "AlternatingProjections.jl" begin
    # Write your own tests here.

    S = ConstrainedBySupport([true, false,true])
    A = ConstrainedByAmplitude([1, sqrt(2), 5])
    x = [1, 2, 3]
    y = [2im, -2. + 2im, 6 - 8im]
    @test project(x, S) == [1, 0, 3]
    @test project(y, A) ≈  [im, -1 + im, 3 - 4im]

    y = zeros(ComplexF64,10,10)  # for ComplexF64 increas the number of iterations
    y[1:5,1:5] = randn(ComplexF32, 5,5)
    Y = fft(y)
    # z = apsolve(abs.(y),abs.(Y), GS, maxit =3000,maxϵ = 1e-18)
    # @test abs.(z) ≈ abs.(y)
    # @test abs.(fft(z)) ≈  abs.(Y)


    # p= PR(abs.(y),abs.(Y))
    # gs=AP(3000,1e-6)
    # z= solve(p, one(y),gs)
    # @test abs.(z) ≈ abs.(y)
    # z, h, snap = solve(p, one(y),gs, true, [1, 5, 10])
    # @test length(h) == 3000
    # @test size(snap) == (10,10,3)

end


@testset "FourierTransformedSet" begin
    a = ConstrainedByAmplitude(ones(5,5))
    A = FourierTransformedSet(a);
    x = rand(Complex{Float64}, 5,5);
    x = project!(x,a)
    @test abs.(x) ≈ ones(5,5)
    
    y = similar(x);
    project!(y,x,A)
    @test abs.(fft(y)) ≈ ones(5,5)

end