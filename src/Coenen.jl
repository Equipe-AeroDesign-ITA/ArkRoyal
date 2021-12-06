module Coenen

	include("Utils.jl")

	export H₁, CE, laminar_H₁

	H_validity_limits = (1.1, 2.7)

	_H₁(H::Real) = 2.0 * H / (H - 1.0)
	@bound_to_limits _H₁ 1.1 2.7 H₁

	_CE(H₁::Real) = 0.0306 * (H₁ - 3.0) ^ (- 0.6169)
	@bound_to_limits _CE 3.1 12.0 CE

	Cf(H::Real, Reθ::Real) = 0.3 * exp(- 1.33 * H) / (
		(log10(Reθ + 1.0)) ^ (1.74 + 0.31 * H)
	)

	_laminar_H₁(H::Real) = 1.2024 * H ^ 4 - 15.5874 * H ^ 3 + 75.6505 * H ^ 2 - 163.5063 * H + 142.7506
	@bound_to_limits _laminar_H₁ 2.0 4.0 laminar_H₁

end
