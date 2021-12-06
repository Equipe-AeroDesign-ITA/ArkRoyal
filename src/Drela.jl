module Drela

	include("Utils.jl")

	export dn!dReθ, Reθ0

	H_validity_limits = (1.1, 4.0)

	_dn!dReθ(H::Real) = 0.01 * sqrt(
		(
			2.4 * H - 3.7 + 2.5 * tanh(1.5 * H - 4.65)
		) ^ 2 + 0.25
	)
	@bound_to_limits _dn!dReθ 2.0 4.0 dn!dReθ

	_Reθ0(H::Real) = 10.0 ^ min(
		(1.415 / (H - 1.0) - 0.485) * tanh(
			20.0 / (H - 1.0) - 12.9
		) + 3.295 / (H - 1.0) + 0.44,
		4.0
	)
	@bound_to_limits _Reθ0 2.0 4.0 Reθ0

	const turbulent_dn!dReθ = dn!dReθ(3.5)

end
