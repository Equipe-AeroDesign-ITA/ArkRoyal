export NACA4

const default_discretization = begin
	θ = collect(
		LinRange(0.0, π, 200)
	)

	@. 1.0 - sin(θ)
end

"""
NACA 4 digit airfoil generation, given camber, position of max. camber,
and thickness percentage
"""
function NACA4(m::Real, p::Real, yt::Real; xs = default_discretization, kwargs...)

	a0 = 0.2969
	a1 = -0.1260
	a2 = -0.3516
	a3 = 0.2843
    a4 = -0.1015

	yts = @. (yt / 0.2) * (a0 * sqrt(xs + 1e-10) + a1 * xs + a2 * xs ^ 2 + a3 * xs ^ 3 + a4 * xs ^ 4)

	pts = Matrix{eltype(xs)}(undef, length(xs), 2)

	yc = x -> (
		x > p ?
		m * ((1.0 - 2 * p) + 2 * p * x - x ^ 2) / ((1.0 - p) ^ 2 + 1e-10) :
		m * (2.0 * p * x - x ^ 2) / (p ^ 2 + 1e-10)
	)
	ycp = x -> ForwardDiff.derivative(yc, x)

	for i = 1:length(xs)
		is_extra = (
			i == 1 ?
			true :
			(
				i == length(xs) ?
				false :
				(
					xs[i + 1] - xs[i - 1] <= 0.0
				)
			)
		)

		θ = atan(ycp(xs[i]))

		pts[i, :] .= (
			is_extra ?
			[xs[i] - yts[i] * sin(θ), yc(xs[i]) + yts[i] * cos(θ)] :
			[xs[i] + yts[i] * sin(θ), yc(xs[i]) - yts[i] * cos(θ)]
		)
	end

	pts

end

"""
NACA 4 digit airfoil generation
"""
function NACA4(nacadigs::Int; xs = default_discretization)

	m = (nacadigs ÷ 1000) / 100.0
	p = ((nacadigs ÷ 100) % 10) / 10.0
	yt = (nacadigs % 100) / 100.0

	NACA4(m, p, yt; xs = xs)

end
