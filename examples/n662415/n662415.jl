using ArkRoyal
using DelimitedFiles

afl = Airfoil("n662415.dat")

αs = collect(0.0:1.0:12.0)

polar = vcat(
	[
		begin
			res = solve(afl; α = α)

			[res.CL res.CD res.Cm]
		end for α in αs
	]...
)

writedlm(
	"n662415.plr",
	[αs polar]
)

res = solve(afl; α = 0.0)

writedlm(
	"n662415.cp",
	[afl.pts[:, 1] res.Cₚ res.θ res.H res.N res.Cf afl.pts[:, 2]]
)
