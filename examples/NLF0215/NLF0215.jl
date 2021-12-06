using ArkRoyal
using DelimitedFiles

afl = Airfoil("NLF0215.dat")

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
	"NLF0215.plr",
	[αs polar]
)

res = solve(afl; α = 0.0)

writedlm(
	"NLF0215.cp",
	[afl.pts[:, 1] res.Cₚ res.θ res.H res.N res.Cf afl.pts[:, 2]]
)
