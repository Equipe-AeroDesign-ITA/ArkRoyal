using ArkRoyal

using ReverseDiff
using ForwardDiff

function drag(x::AbstractVector)

	m, p, yt = x

	afl = Airfoil(
		NACA4(m, p, yt)
	)

	u = solve_inviscid(afl)

	return sum(u)

end

x = [0.03, 0.6, 0.12]
@time ReverseDiff.gradient(drag, x)
@time ReverseDiff.gradient(drag, x)
@time ReverseDiff.gradient(drag, x)
@time ReverseDiff.gradient(drag, x)
@time ReverseDiff.gradient(drag, x)

@time ForwardDiff.gradient(drag, x)
@time ForwardDiff.gradient(drag, x)
@time ForwardDiff.gradient(drag, x)
@time ForwardDiff.gradient(drag, x)
@time ForwardDiff.gradient(drag, x)