using DoubleFloats
struct basis
    Z::Int128# we will be only studying He atom
    n::Int128
    l::Int128
    m::Int128
    alpha::Double64
    beta::Double64
    gamma::Double64
end

function test_basis()
    alpha = 2.0 * 1.8
    beta = 2.0 * 1.8
    gamma = 0.0
    test_basis = basis[]
    push!(test_basis, basis(2,0, 0, 0, alpha, beta, gamma))
    push!(test_basis, basis(2,1, 1, 0, alpha, beta, gamma))
    push!(test_basis, basis(2,0, 0, 1, alpha, beta, gamma))
    return test_basis
end

function lambda_N(N::Int128, alpha::Double64, gamma::Double64)
    lambda_n = basis[]
    for n in 0:N
        for l in 0:N-n
            for m in 0:N-l-n
                bs = basis(2,n, l, m, alpha, alpha, gamma)
                push!(lambda_n, bs)
            end
        end
    end
    return lambda_n
end
println(test_basis())