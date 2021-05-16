using CoupledDipole
using Test

@testset "CoupledDipole" begin
    q = cubature_sphere(8, "gl")
    @test length(q.nodes) == 8
end


@testset "approx equal" begin

    @test π ≈ 3.14 atol = 0.01

end

@testset "trigonometric identities" begin
    θ = 2 / 3 * π
    @test sin(-θ) ≈ -sin(θ)
    @test cos(-θ) ≈ cos(θ)
    @test sin(2θ) ≈ 2 * sin(θ) * cos(θ)
    @test cos(2θ) ≈ cos(θ)^2 - sin(θ)^2
end;


@testset verbose = true "Foo Tests" begin
    @testset "Animals" begin
        @test 9 == 9
        @test "dog" == "dog"
    end
    @testset "Animals" begin
        @test 9 == 9
        @test "dog" == "dog"
    end
end;
