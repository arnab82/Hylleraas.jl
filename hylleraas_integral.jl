using LinearAlgebra
include("./K_integral.jl")

@inline S_ij(bf1::basis, bf2::basis)::BigFloat= K(bf1.n + bf2.n, bf1.l + bf2.l, bf1.m + bf2.m,
                                                 bf1.alpha, bf1.beta, bf1.gamma)

function Vne_ij(bf1::basis, bf2::basis)::BigFloat
    one = BigInt(1)
    k1 = K(bf1.n + bf2.n - one, bf1.l + bf2.l, bf1.m + bf2.m, bf1.alpha, bf1.beta, bf1.gamma)
    k2 = K(bf1.n + bf2.n, bf1.l + bf2.l - one, bf1.m + bf2.m, bf1.alpha, bf1.beta, bf1.gamma)
    return -one * bf1.Z * (k1 + k2)
end

function Vee_ij(bf1::basis, bf2::basis)::BigFloat 
    one=BigInt(1)
    return K(bf1.n + bf2.n, bf1.l + bf2.l, bf1.m + bf2.m - one,
                                                 bf1.alpha, bf1.beta, bf1.gamma)
end

function T_ij(bf1::basis, bf2::basis)::BigFloat

    alpha = bf1.alpha
    beta = bf1.beta
    gamma = bf1.gamma
    nj = bf2.n
    lj = bf2.l
    mj = bf2.m

    screen_prefactor(prefactor::BigFloat)::Bool = prefactor != 0.0

    function Knlm(prefactor::BigFloat, n::BigInt, l::BigInt, m::BigInt)::BigFloat
        if screen_prefactor(prefactor)
            return prefactor * K(bf1.n + bf2.n + n, bf1.l + bf2.l + l, bf1.m + bf2.m + m, alpha, beta, gamma)
        else
            return BigFloat(0.0)
        end
    end

    zero = BigInt(0)
    one = BigInt(1)
    two = BigInt(2)
    half = BigFloat(0.5)
    one_fourth = BigFloat(0.25)
    one_eighth = BigFloat(0.125)

    Tij = -one_eighth* ((alpha^two) +(beta^two) + (gamma^two)) * S_ij(bf1, bf2)

    Tij += Knlm((half * nj * alpha + half * alpha), -one, zero, zero)
    Tij += Knlm((half * lj * beta + half * alpha) ,zero, -one, zero)
    Tij += Knlm((mj * gamma + gamma ),zero, zero, -one)

    Tij -= Knlm((half * nj * (nj - one) + nj),-two, zero, zero)
    Tij -= Knlm((half * lj * (lj - one) + lj),zero, -two, zero)
    Tij -= Knlm(BigFloat(mj * (mj - one) + two * mj),zero, zero, -two)

    Tij -= Knlm((one_eighth * alpha * gamma),-one, zero, one)
            + Knlm((one_eighth * alpha * gamma),one, zero, -one)
            - Knlm((one_eighth * alpha * gamma),-one, two, -one)

    Tij -= Knlm((one_eighth * beta * gamma),zero, -one, one)
            + Knlm((one_eighth * beta * gamma),zero, one, -one)
            - Knlm((one_eighth * beta * gamma),two, -one, -one)

    Tij += Knlm((one_fourth * nj * gamma),-two, zero, one)
            + Knlm((one_fourth * nj * gamma),zero, zero, -one)
            - Knlm((one_fourth * nj * gamma),-two, two, -one)

    Tij += Knlm((one_fourth * mj * alpha),-one, zero, zero)
            + Knlm((one_fourth * mj * alpha),one, zero, -two)
            - Knlm((one_fourth * mj * alpha),-one, two, -two)

    Tij -= Knlm((half * nj * mj),-two, zero, zero)
            + Knlm((half * nj * mj),zero, zero, -two)
            - Knlm((half * nj * mj),-two, two, -two)

    Tij += Knlm((one_fourth * lj * gamma),zero, -two, one)
            + Knlm((one_fourth * lj * gamma),zero, zero, -one)
            - Knlm((one_fourth * lj * gamma),two, -two, -one)

    Tij += Knlm((one_fourth * mj * beta),zero, -one , zero)
            + Knlm((one_fourth * mj * beta),zero, one, -two)
            - Knlm((one_fourth * mj * beta),two,-one, -two)

    Tij -= Knlm((half * lj * mj),zero, -two, zero)
            + Knlm((half * lj * mj),zero, zero, -two)
            - Knlm((half * lj * mj),two, -two, -two)

    return Tij
end
"""
Computes integral matrix given a vector of basis functions

Args:
- bfn: vector of basis functions
- integral_ij: function that computes the elements of the matrix

Returns:
- matrix: computed integral matrix
"""
function compute_integral(bfn::Vector{basis}, integral_ij)
    n = length(bfn)
    matrix = zeros(n, n)
    for i = 1:n
        for j = 1:n
            matrix[i, j] = integral_ij(bfn[i], bfn[j])
        end
    end
    return matrix
end
