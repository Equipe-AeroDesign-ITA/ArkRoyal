export solve, AirfoilSolution

struct AirfoilSolution
	afl::Airfoil
	u::AbstractVector
	Cₚ::AbstractVector
	θ::AbstractVector
	H::AbstractVector
	N::AbstractVector
	Cf::AbstractVector
	CL::Real
	CD::Real
	Cm::Real
	converged::Bool
end

#=
"""
Obtain maximum variation in `Cₚ` between upper and lower surface
"""
function _max_dCp(
	afl::Airfoil,
	Cp::AbstractVector
)

	mdCp = 0.0

	j = size(afl.pts, 1)
	for i = 2:size(afl.pts, 1)
		if afl.pts[i + 1, 1] > afl.pts[i, 1] || afl.pts[j - 1, 1] > afl.pts[j, 1]
			return mdCp
		end

		while !(afl.pts[j, 1] >= afl.pts[i, 1] && afl.pts[j - 1, 1] <= afl.pts[i, 1])
			if afl.pts[j, 1] >= afl.pts[i, 1]
				j -= 1
			else
				j += 1
			end
		end

		η = (afl.pts[j, 1] - afl.pts[i, 1]) / (afl.pts[j, 1] - afl.pts[j - 1, 1])

		mdCp = max(
			mdCp,
			abs(
				Cp[j] * (1.0 - η) + Cp[j - 1] * η - Cp[i]
			)
		)
	end

	mdCp

end
=#

"""
Solve airfoil flow
"""
function solve(
	afl::Airfoil;
	α::Real = 0.0,
	Re::Real = 1e6,
	M∞::Real = 0.0,
	Ncr::Real = 9.0,
	n_iter::Int64 = 50,
	tol::Real = 0.0,
	ω::Real = 0.2,
	stall_correction::Symbol = :Eppler,
	kwargs...
)

	u = solve_inviscid(afl; α = α)

	q0 = [0.0, 2.5, 0.0]

	i1 = i2 = 0
	η = 0.0
	for i = 1:(length(u) - 1)
		if u[i] * u[i + 1] < 0.0
			i1 = i
			i2 = i + 1

			η = abs(u[i1]) / (abs(u[i1]) + abs(u[i2]))
		end
	end

	vs = Matrix{eltype(u)}(undef, 3, length(u))

	converged = true
	dr = 0.0

	# i1
	vs[:, i1] .= march(
		q0, q0, 0.0, u[i1], (afl.sV[i1] - afl.sV[i2]) * η;
		M∞ = M∞,
		Re = Re,
		Ncr = Ncr,
		ω = ω,
		kwargs...
	)
	for nit = 2:n_iter
		θo = vs[1, i1]

		vs[:, i1] .= march(
			q0, vs[:, i1], 0.0, u[i1], (afl.sV[i1] - afl.sV[i2]) * η;
			M∞ = M∞,
			Re = Re,
			Ncr = Ncr,
			ω = ω,
			kwargs...
		)

		if abs(vs[1, i1] - θo) / (θo + 1e-12) < tol
			break
		end

		if nit == n_iter
			converged = false
		end
	end
	
	for i = (i1-1):-1:1
		vs[:, i] .= march(
			vs[:, i + 1], vs[:, i + 1], u[i + 1], u[i], afl.sV[i] - afl.sV[i + 1];
			M∞ = M∞,
			Re = Re,
			Ncr = Ncr,
			ω = ω,
			kwargs...
		)
		for nit = 2:n_iter
			θo = vs[1, i]

			vs[:, i] .= march(
				vs[:, i + 1], vs[:, i], u[i + 1], u[i], afl.sV[i] - afl.sV[i + 1];
				M∞ = M∞,
				Re = Re,
				Ncr = Ncr,
				ω = ω,
				kwargs...
			)

			if abs(vs[1, i] - θo) / (θo + 1e-12) < tol
				break
			end

			if nit == n_iter
				converged = false
			end
		end
	end

	# i2
	vs[:, i2] .= march(
		q0, q0, 0.0, u[i2], (afl.sV[i2] - afl.sV[i1]) * (1.0 - η);
		M∞ = M∞,
		Re = Re,
		Ncr = Ncr,
		ω = ω,
		kwargs...
	)
	for nit = 2:n_iter
		θo = vs[1, i2]

		vs[:, i2] .= march(
			q0, vs[:, i2], 0.0, u[i2], (afl.sV[i2] - afl.sV[i1]) * (1.0 - η);
			M∞ = M∞,
			Re = Re,
			Ncr = Ncr,
			ω = ω,
			kwargs...
		)

		if abs(vs[1, i2] - θo) / (θo + 1e-12) < tol
			break
		end

		if nit == n_iter
			converged = false
		end
	end
	
	for i = (i2+1):length(u)
		vs[:, i] .= march(
			vs[:, i - 1], vs[:, i - 1], u[i - 1], u[i], afl.sV[i] - afl.sV[i - 1];
			M∞ = M∞,
			Re = Re,
			Ncr = Ncr,
			ω = ω,
			kwargs...
		)
		for nit = 2:n_iter
			θo = vs[1, i]

			vs[:, i] .= march(
				vs[:, i - 1], vs[:, i], u[i - 1], u[i], afl.sV[i] - afl.sV[i - 1];
				M∞ = M∞,
				Re = Re,
				Ncr = Ncr,
				ω = ω,
				kwargs...
			)

			if abs(vs[1, i] - θo) / (θo + 1e-12) < tol
				break
			end

			if nit == n_iter
				converged = false
			end
		end
	end

	# compressibility corrections from Karman-Tsien's correction
	u = @. u * (1.0 - M∞ ^ 2 / (2.0 - M∞ ^ 2)) / (1.0 - M∞ ^ 2 * u ^ 2 / (2.0 - M∞ ^ 2))

	Cₚ = @. 1.0 - u ^ 2

	θ = vs[1, :]
	H = vs[2, :]
	N = vs[3, :]

	Cf = [
		begin
			t = σ(
				N[i] - Ncr
			)

			Reθ = max(10.0, abs(θ[i] * u[i] * Re))

			(1.0 - t) * Thwaites.Cf(Thwaites.λ(H[i]) / Reθ, Reθ) +
			t * Coenen.Cf(H[i], Reθ)
		end for i = 1:length(u)
	]

	θTE = (θ[1] + θ[end])
	uTE = (abs(u[1]) + abs(u[end])) / 2
	HTE = (H[1] + H[end]) / 2
	Havg = (HTE + 1.0) / 2

	CL = 0.0
	CD = 2 * θTE * uTE ^ (Havg + 2.0)
	Cm = 0.0

	ŷ = [
		(- sind(α)), cosd(α)
	]
	x̂ = [cosd(α), sind(α)]

	for i = 1:size(afl.pts, 1)
		s = (
			i == 1 ?
			(afl.pts[2, :] .- afl.pts[1, :]) ./ 2 :
			(
				i == size(afl.pts, 1) ?
				(afl.pts[end, :] .- afl.pts[end - 1, :]) ./ 2 :
				(
					afl.pts[i + 1, :] .- afl.pts[i - 1, :]
				) ./ 2
			)
		)
		n = [- s[2], s[1]]

		F = s .* Cf[i] .* (u[i] > 0.0 ? 1.0 : - 1.0) .+ n .* Cₚ[i]

		CL += F ⋅ ŷ
		# CD += F ⋅ x̂
		Cm += afl.pts[i, 2] * F[1] - (afl.pts[i, 1] - 0.25) * F[2]
	end

	xdet_upper = afl.pts[1, 1]
	ydet_upper = afl.pts[1, 2]
	xdet_lower = afl.pts[end, 1]
	ydet_lower = afl.pts[end, 2]

	for i = 1:i1
		if N[i] > Ncr
			if H[i] >= 2.8 && H[i + 1] < 2.8
				η = (H[i] - 2.8) / (H[i] - H[i + 1])

				xdet_upper = afl.pts[i + 1, 1] * η + afl.pts[i, 1] * (1.0 - η)
				ydet_upper = afl.pts[i + 1, 2] * η + afl.pts[i, 2] * (1.0 - η)
			end
		end
	end

	for i = length(N):-1:i2
		if N[i] > Ncr
			if H[i] >= 2.8 && H[i - 1] < 2.8
				η = (H[i] - 2.8) / (H[i] - H[i - 1])

				xdet_lower = afl.pts[i - 1, 1] * η + afl.pts[i, 1] * (1.0 - η)
				ydet_lower = afl.pts[i - 1, 2] * η + afl.pts[i, 2] * (1.0 - η)
			end
		end
	end

	αr = deg2rad(α)

	if α > 0.0
		if stall_correction == :Eppler
			θTE = (
				xdet_upper < afl.pts[1, 1] ?
				atan((ydet_upper - afl.pts[1, 2]), - (xdet_upper - afl.pts[1, 1])) :
				0.0
			)

			CD += 0.2 * sin(αr + θTE) * (afl.pts[1, 1] - xdet_upper) ^ 2

			ΔCL = CL * (αr + θTE) * π * (afl.pts[1, 1] - xdet_upper)

			if ΔCL > 0.0
				CL = CL - ΔCL
			else
				CL *= (1.0 - sin(αr) * (afl.pts[1, 1] - xdet_upper))
			end

			Cm -= sin(αr) * (afl.pts[1, 1] - xdet_upper) * (
				0.5 * (1.0 + xdet_upper) - 0.25
			)
		elseif stall_correction == :CalcFoil
			CD += abs(
				sin(αr) ^ 2 * (afl.pts[1, 1] - xdet_upper) ^ 2 + 0.025 * cos(αr) * (afl.pts[1, 1] - xdet_upper) ^ 2
			)
			CL *= (1.0 - 0.2 * (afl.pts[1, 1] - xdet_upper))
			Cm *= 0.9 * xdet_lower ^ 2 * xdet_upper ^ 2
		else
			@error "Uknown stall correction mode $stall_correction"
		end
	else
		if stall_correction == :Eppler
			θTE = (
				xdet_lower < afl.pts[end, 1] ?
				atan((ydet_lower - afl.pts[end, 2]), - (xdet_lower - afl.pts[end, 1])) :
				0.0
			)

			CD -= 0.2 * sin(αr + θTE) * (afl.pts[end, 1] - xdet_lower) ^ 2

			ΔCL = CL * (αr + θTE) * π * (afl.pts[end, 1] - xdet_lower)

			if ΔCL < 0.0
				CL = CL - ΔCL
			else
				CL *= (1.0 - sin(αr) * (afl.pts[end, 1] - xdet_lower))
			end

			Cm -= sin(αr) * (afl.pts[end, 1] - xdet_lower) * (
				0.5 * (1.0 + xdet_lower) - 0.25
			)
		elseif stall_correction == :CalcFoil
			CD += abs(
				sin(αr) ^ 2 * (afl.pts[end, 1] - xdet_lower) ^ 2 + 0.025 * cos(αr) * (afl.pts[end, 1] - xdet_lower) ^ 2
			)
			CL *= (1.0 - 0.2 * (afl.pts[end, 1] - xdet_lower))
			Cm *= 0.9 * xdet_lower ^ 2 * xdet_upper ^ 2
		else
			@error "Uknown stall correction mode $stall_correction"
		end
	end

	A = (
		stall_correction == :Eppler ?
		30.0 :
		20.0
	)

	dCp = (maximum(Cₚ) - minimum(Cₚ)) # _max_dCp(afl, Cₚ)
	CL *= (1.0 / ((dCp / A) ^ 2 + 1.0))

	AirfoilSolution(
		afl,
		u,
		Cₚ,
		θ,
		H,
		N,
		Cf,
		CL,
		CD,
		Cm,
		converged || (tol == 0.0)
	)

end
