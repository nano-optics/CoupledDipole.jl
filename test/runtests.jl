using CoupledDipole
using Test

@testset "CoupledDipole.jl" begin
    q = CoupledDipole.cubature_sphere(3, "gl")
    @test length(q.nodes) == 3
end

