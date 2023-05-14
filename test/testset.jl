using Hyllerus
using Test

@testset "YourPackageName.jl" begin
    @test Hyllerus.greet_your_package_name() == "Hello hyllerus!"
    @test Hyllerus.greet_your_package_name() != "Hello world!"
end
