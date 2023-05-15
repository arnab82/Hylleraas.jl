#using hylleraas
include("./../src/hylleraas.jl")
using Test

@testset "hylleraas.jl" begin
    @test hylleraas.main() == -2.90372420
end
