using LinearAlgebra
using SpecialFunctions
include("./basis.jl")
pi=BigFloat(π)
function K(n::BigInt, l::BigInt, m::BigInt, alpha::BigFloat, beta::BigFloat, gamma::BigFloat)::BigFloat
    one = BigInt(1)
    two = BigInt(2)
    sixteen = BigInt(16)
    prefactor = sixteen * pi^two * factorial(n+one) * factorial(l+one) * factorial(m+one)
    k = BigFloat(0.0)
    for a = 0:n+1, b = 0:l+1, c = 0:m+1
        numerator = binomial(l+one-b+a, a) * binomial(m+one-c+b, b) * binomial(n+one-a+c, c)
        denominator = (alpha+beta)^(l-b+a+two) * (alpha+gamma)^(n-a+c+two) * (beta+gamma)^(m-c+b+two)
        k += numerator / denominator
    end
    return prefactor * k
end
binomial(BigInt(13),BigInt(4))