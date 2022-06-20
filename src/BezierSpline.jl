function Bezier_spline(
	pts::AbstractMatrix,
	η::Real
)

	if size(pts, 1) == 1
		return pts[1, :]
	end

	(1.0 - η) * Bezier_spline(pts[1:(end - 1), :], η) + η * Bezier_spline(pts[2:end, :], η)

end

Bezier_spline(
	pts::AbstractMatrix,
	η::AbstractVector
) = permutedims(
	hcat(
		[
			Bezier_spline(pts, p) for p in η
		]...
	)
)
