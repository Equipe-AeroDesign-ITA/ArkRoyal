using ArkRoyal
using DelimitedFiles

afl = Airfoil("GAW1.dat")

αs = collect(0.0:1.0:15.0)

converged = Bool[]
polar = vcat(
	[
		begin
			res = solve(afl; α = α)

			push!(converged, res.converged)

			[res.CL res.CD res.Cm]
		end for α in αs
	]...
)
αs = αs[converged]
polar = polar[converged, :]

writedlm(
	"GAW1.plr",
	[αs polar]
)

res = solve(afl; α = 0.0, n_iter = 15)

writedlm(
	"GAW1.cp",
	[afl.pts[:, 1] res.Cₚ res.θ res.H res.N res.Cf afl.pts[:, 2]]
)
