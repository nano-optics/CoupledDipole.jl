using CoupledDipole
using Test

@testset "CoupledDipole.jl" begin
    q = CoupledDipole.cubature_sphere(8, "gl")
    @test length(q.nodes) == 8
end

