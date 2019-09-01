using AlternatingProjections
using Test

@testset "AlternatingProjections.jl" begin
    # Write your own tests here.

    S = SupportConstraint([true, false,true])
    A = AmplitudeConstraint([1, sqrt(2), 5])
    x = [1, 2, 3]
    y = [2im, -2 + 2im, 6 - 8im]
    @test project(x, S) == [1, 0, 3]
    @test project(y, A) ≈  [im, -1 + im, 3 - 4im]
end
