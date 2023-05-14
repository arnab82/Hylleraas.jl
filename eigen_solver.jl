
using LinearAlgebra

function solve_HC_SCE(H::Matrix, S::Matrix, verbose::Bool=false)
    eigenvalue, eigenvec= eigen(H, S)
    gs_energy = eigenvalue
    if verbose
        println("ground state energy = ", gs_energy)
    end
    return gs_energy
end



