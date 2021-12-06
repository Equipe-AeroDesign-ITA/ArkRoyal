module Thwaites

	include("Utils.jl")

	export Cf, H, Λ

	λ_validity_limits = (- 0.1, 0.1)

	_H(λ::Real) = 2.61 - 4.1 * λ + 14.0 * λ ^ 3 + 0.56 * λ ^ 2 / (λ + 0.18) ^ 2
	@bound_to_limits _H -0.1 0.15 H

	_T(λ::Real) = 0.22 + 1.52 * λ - 5.0 * λ ^ 3 - 0.072 * λ ^ 2 / (λ + 0.18) ^ 2
	@bound_to_limits _T -0.1 0.15 T

	Cf(Λ::Real, Reθ::Real) = 2 * T(Λ * Reθ) / (Reθ + 1.0)
	H(Λ::Real, Reθ::Real) = H(Λ * Reθ)

	_λ(H::Real) = - 0.0715 * H ^ 3 + 0.7696 * H ^ 2 - 2.7866 * H + 3.3028
	@bound_to_limits _λ 2.0 4.0 λ

end
