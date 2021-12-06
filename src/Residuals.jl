σ(x::Real) = 1.0 / (exp(clamp(- x, -700.0, 700.0)) + 1.0)

"""
Get residuals for a combination of two successive 
mesh nodes
"""
function residual(
	Qm::AbstractVector,
	um::Real,
	Q::AbstractVector,
	u::Real,
	δx::Real;
	Ncr::Real = 9.0,
	Re::Real = 1e6,
	M∞::Real = 0.0,
	γ::Real = 1.4
)

	θm, Hm, Nm = Qm
	θ, H, N = Q

	lam, turb = begin
		t = σ(N - Ncr)

		1.0 - t, t
	end
	lam_m, turb_m = begin
		t = σ(Nm - Ncr)

		1.0 - t, t
	end

	δst = H * θ

	Λ = θ * (u - um) / ((abs(u) + 1e-7) * δx)

	Reθ = max(10.0, abs(u * θ * Re))

	Cf = (
		lam * Thwaites.Cf(Thwaites.λ(H) / Reθ, Reθ) + 
		turb * Coenen.Cf(H, Reθ)
	)

	H₁m = (
		lam_m * Coenen.laminar_H₁(Hm) + 
		turb_m * Coenen.H₁(Hm)
	)
	H₁ = (
		lam * Coenen.laminar_H₁(H) + 
		turb * Coenen.H₁(H)
	)

	CE = Coenen.CE(H₁)

	Mem = abs(M∞ * um)
	Me = abs(M∞ * u)

	ρm = ((1.0 + (γ - 1.0) * M∞ ^ 2 / 2) / (1.0 + (γ - 1.0) * Mem ^ 2 / 2)) ^ (1.0 / (γ - 1.0))
	ρ = ((1.0 + (γ - 1.0) * M∞ ^ 2 / 2) / (1.0 + (γ - 1.0) * Me ^ 2 / 2)) ^ (1.0 / (γ - 1.0))

	τ = ρ * u * abs(u) * Cf / 2

	rmomx = (ρ * θ * u ^ 2 - ρm * θm * um ^ 2) / δx + δst * ρ * u * (u - um) / δx - τ
	rentr = (
		(u * θ * H₁ - um * θm * H₁m) / δx - abs(u) * CE
	) * turb + (Thwaites.H(Λ, Reθ) - H) * lam

	∂Reθ = (
		max(abs(θ * u * Re), 10.0) - max(abs(θm * um * Re), 10.0)
	) / δx

	Reθ0 = (
		N > 1.0 ? 0.0 : Drela.Reθ0(H)
	)
	dn!dr = (
		Reθ > Reθ0 ? lam * Drela.dn!dReθ(H) + turb * Drela.turbulent_dn!dReθ : 0.0
	)

	rTS = dn!dr * ∂Reθ - (N - Nm) / δx

	Hl = 2.0 * lam + 1.1 * turb
	Hu = 4.0 * lam + 3.0 * turb

	[rmomx, rentr, rTS], [Hl, Hu]

end

"""
March variables from a station to another
"""
function march(
	qim1::AbstractVector,
	q::AbstractVector,
	uim1::Real,
	u::Real,
	δx::Real;
	kwargs...
)

	f = x -> begin
		r, _ = residual(
			qim1,
			uim1,
			x,
			u,
			δx;
			kwargs...
		)

		r
	end

	J = ForwardDiff.jacobian(f, q)
	r, i = residual(
		qim1,
		uim1,
		q,
		u,
		δx;
		kwargs...
	)

	q = q - J \ r

	θ, H, N = q
	θm, Hm, Nm = qim1

	if θ < θm
		θ = θm
	end

	if H < i[1]
		H = i[1]
	elseif H > i[2]
		H = i[2]
	end

	if N < Nm
		N = Nm
	end

	[θ, H, N]

end
