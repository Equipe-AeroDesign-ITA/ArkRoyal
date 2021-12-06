using ArkRoyal
using DelimitedFiles

afl = Airfoil("GAW1.dat")

αs = collect(0.0:1.0:15.0)

polar = vcat(
	[
		begin
			res = solve(afl; α = α, n_iter = 15)

			[res.CL res.CD res.Cm]
		end for α in αs
	]...
)

writedlm(
	"GAW1.plr",
	[αs polar]
)

res = solve(afl; α = 0.0, n_iter = 15)

writedlm(
	"GAW1.cp",
	[afl.pts[:, 1] res.Cₚ res.θ res.H res.N res.Cf afl.pts[:, 2]]
)
