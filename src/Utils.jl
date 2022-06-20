using ForwardDiff

"""
Bound function to limits and linearly extrapolate past boundaries
"""
macro bound_to_limits(fs::Symbol, L::Real, H::Real, fn::Symbol)

	f = eval(fs)

	fL = f(L)
	fH = f(H)
	fLp = ForwardDiff.derivative(f, L)
	fHp = ForwardDiff.derivative(f, H)

	q = :($fn(x::Real) = (x < $L ? ($fL + (x - $L) * $fLp) : (x > $H ? ($fH + (x - $H) * $fHp) : $fs(x))))

	return q

end
