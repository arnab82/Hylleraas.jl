module hylleraas
include("./hylleraas_integral.jl")
include("./eigen_solver.jl")
export main
function roundDecimal!(number, digit)
    return round(number; digits=digit)
end

function main()
	setprecision(15)

	alpha = Double64(3.6)  # actually 2 * 1.8
	gamma = Double64(0.0)
	println("N\tbasis size\tenergy (Eh)")
	e=Double64(0.0)
	for N=Int128(0):Int128(10)
    	b_set = lambda_N(N,alpha,gamma)
    	S = compute_integral(b_set, S_ij)
    	Vne = compute_integral(b_set, Vne_ij)
    	Vee = compute_integral(b_set, Vee_ij)
    	Te = compute_integral(b_set, T_ij)
    	H = Vne + Vee + Te
    	e = solve_HC_SCE(H, S, false)
    	println("$N\t$(length(b_set))\t$e")
	end
	return roundDecimal!(e,8)
end
end
