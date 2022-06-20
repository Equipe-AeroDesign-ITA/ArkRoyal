using ArkRoyal
using DelimitedFiles

afl = Airfoil(NACA4(0012))

αs = collect(0.0:1.0:16.0)

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
	"n0012.plr",
	[αs polar]
)

res = solve(afl; α = 0.0)

writedlm(
	"n0012.cp",
	[afl.pts[:, 1] res.Cₚ res.θ res.H res.N res.Cf afl.pts[:, 2]]
)
