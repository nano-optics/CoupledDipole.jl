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
        @test foo("cat") == 9
        @test foo("dog") == foo("cat")
    end
    @testset "Arrays $i" for i = 1:3
        @test foo(zeros(i)) == i^2
        @test foo(fill(1.0, i)) == i^2
    end
end;
