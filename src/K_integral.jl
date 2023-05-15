using LinearAlgebra
using DoubleFloats
using SpecialFunctions
include("./basis.jl")
pi=Double64(π)
function K(n::Int128, l::Int128, m::Int128, alpha::Double64, beta::Double64, gamma::Double64)::Double64
    one = Int128(1)
    two = Int128(2)
    sixteen = Int128(16)
    prefactor = sixteen * pi^two * factorial(n+one) * factorial(l+one) * factorial(m+one)
    k = Double64(0.0)
    for a = 0:n+1, b = 0:l+1, c = 0:m+1
        numerator = binomial(l+one-b+a, a) * binomial(m+one-c+b, b) * binomial(n+one-a+c, c)
        denominator = (alpha+beta)^(l-b+a+two) * (alpha+gamma)^(n-a+c+two) * (beta+gamma)^(m-c+b+two)
        k += numerator / denominator
    end
    return prefactor * k
end