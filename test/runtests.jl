using Test
using DIVAnd_hfradar

@testset "adjoint" begin
    include("test_adjoint.jl")
end

@testset "adjoint misfit" begin
    include("test_adjoint_misfit.jl")
end

@testset "DIVAndrun analysis" begin
    include("test_DIVAndrun.jl")
end
