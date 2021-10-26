using AlternatingProjections, FFTW
# using LinearAlgebra
using Test

@testset "AlternatingProjections.jl" begin
    # Write your own tests here.

    S = ConstrainedBySupport([true, false,true])
    A = ConstrainedByAmplitude([1, sqrt(2), 5])
    x = [1, 2, 3]
    y = [2im, -2. + 2im, 6 - 8im]
    @test project(x, S) == [1, 0, 3]
    @test project(y, A) ≈  [im, -1 + im, 3 - 4im]

    xr=copy(x)
    reflect!(xr,x,S)
    @test xr == [1, -2, 3]

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

    @test AlternatingProjections.getelement(a) == ones(5,5)
    B = zeros(Complex,(5,5)); B[1] +=1;
    @test AlternatingProjections.getelement(A) ≈ B
     
    
end

@testset "Amplitude and Shape Constraints" begin
    a = [ 8 - i^2 - j^2 for i in -2:2, j in -2:2]
    m = [ i >= 0 && j >= 0 for i in -2:2, j in -2:2]

    A = ConstrainedByAmplitude(a)
    B = ConstrainedByAmplitudeMasked(a,m)
    C = ConstrainedByShape(a)
    D = ConstrainedByShapeMasked(a, m)

    x = randn(ComplexF64, 5,5)
    y = project(x,A)
    z = project(x, B)
    v = project(x,C)
    w = project(x, D)

    @test abs.(y) ≈ a
    @test abs.(z)[m] ≈ a[m]
    @test abs.(v)/sum(abs.(v)) ≈ a /sum(a)

    s = abs.(x)[D.mask]' * D.amp[D.mask] / D.n

    @test abs.(w[m])/s ≈ a[m]


end

@testset "ConstrainedByShapeSaturated" begin
    a = map(x -> 1-x^2, -1:0.1:1)
    meas = 2 * a
    measset = ConstrainedByShapeSaturated(meas, meas .<=1)
    meassat = copy(meas)
    meassat[meassat .> 1] .= 1

    p = project(a, measset)
    @test p == a

    p2 = project(meas, measset)
    @test p2 == meas

end


@testset "solve" begin

    # Quick Gerchberg-Saxton
    n = 50
    m = 10

    amp = 0.35
    x = zeros(ComplexF64, n, n); x[1:m,1:m] .= exp.(amp * 2π * im * randn(Float64, m, m))
    a = abs.(x)
    X = fft(x)
    A = abs.(X)
    aset = ConstrainedByAmplitude(a)
    Aset = FourierTransformedSet(ConstrainedByAmplitude(A))
    problem = TwoSetsFP(aset,Aset)
    # start with a good initial guess
    sol = solve(problem, APparam(), x⁰ = x + 1.5 * randn(ComplexF64, n, n), maxit = 1000)
    
    function testsolution()
        phasediff = -pi .+ mod2pi.(pi .+ angle.(sol[1][1:m,1:m]) .- angle.(x[1:m,1:m]))
        phasediff .-= sum(phasediff)/m^2
        phaserms = sqrt(sum(abs2,phasediff))/m
        #phase difference can be visualised 
        # heatmap(phasediff)
        println("phase rms = $phaserms")
        @test phaserms < 1e-6
    end

    testsolution()

    # and the same problem solved with DR method (it requires more iterations and other starting point)
    sol = solve(problem, DRparam(), x⁰ = x - 1.5 * randn(ComplexF64, n, n), maxit = 3000)
    
    testsolution()

    # and with DRAP algoritm
    drap = DRAPparam(missing,500,missing,true,[1],0.1)
    # sol = solve(problem, drap, maxit=1000);
    sol = solve(problem, drap, x⁰ = x + 1.5 * randn(ComplexF64, n, n), maxit=1000);
    testsolution()
end

