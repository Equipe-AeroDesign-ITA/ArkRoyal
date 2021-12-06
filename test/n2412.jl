@assert begin
	afl = Airfoil("n2412.dat")

	αs = collect(-10.0:1.0:15.0)

	polar = vcat(
		[
			begin
				res = solve(afl; α = α)

				[res.CL res.CD res.Cm]
			end for α in αs
		]...
	)

	writedlm(
		"n2412.plr",
		[αs polar]
	)

	res = solve(afl; α = 0.0)

	writedlm(
		"n2412.cp",
		[afl.pts[:, 1] res.Cₚ res.θ res.H res.N res.Cf afl.pts[:, 2]]
	)

	@time res = solve(afl; α = 0.0)
	@time res = solve(afl; α = 0.0)
	@time res = solve(afl; α = 0.0)
	@time res = solve(afl; α = 0.0)
	@time res = solve(afl; α = 0.0)

	true
end
