using ArkRoyal

using ReverseDiff

afl = Airfoil(NACA4(3612))

CD =  α -> begin
	res = solve(afl; α = first(α))

	res.CD
end

@show rad2deg.(ReverseDiff.gradient(CD, [0.0]))
@show CD(0.0)
