using ArkRoyal
using DelimitedFiles

afl = Airfoil("GRK1.dat")

αs = collect(0.0:1.0:14.0)

polar = vcat(
	[
		begin
			res = solve(afl; α = α)

			[res.CL res.CD res.Cm]
		end for α in αs
	]...
)

writedlm(
	"GRK1.plr",
	[αs polar]
)

res = solve(afl; α = 0.0)

writedlm(
	"GRK1.cp",
	[afl.pts[:, 1] res.Cₚ res.θ res.H res.N res.Cf afl.pts[:, 2]]
)
