export Airfoil, solve_inviscid

function influence(
	panel_pts::AbstractMatrix,
	collocation::AbstractVector,
	γs::AbstractVector,
	is_same::Bool
)
	
	x̂ = panel_pts[:, 2] .- panel_pts[:, 1]
	
	L = norm(x̂)
	x̂ ./= L
	
	ŷ = [- x̂[2], x̂[1]]
	
	temp = collocation .- panel_pts[:, 1]
	x = dot(temp, x̂)
	y = dot(temp, ŷ) + 1e-10
	
	temp = collocation .- panel_pts[:, 2]
	x2 = dot(temp, x̂)
	y2 = dot(temp, ŷ) + 1e-10
	
	rj = sqrt(x ^ 2 + y ^ 2)
	rjp1 = sqrt(x2 ^ 2 + y2 ^ 2)
	
	θj = atan(y, x)
	θjp1 = (
		is_same ?
		pi :
		atan(y2, x2)
	)
	
	δ = (γs[2] - γs[1])
	
	u = (y * δ * (log(rjp1) - log(rj)) + 
		(γs[1] * L + δ * x) * (θjp1 - θj)) / (2 * π * L)
	
	w = (- y * δ * (θjp1 - θj) +
		(γs[1] * L + δ * x) * (log(rjp1) - log(rj))) / (2 * π * L) +
		δ / (2 * π)
	
	u * x̂ + w * ŷ
	
end

mutable struct Airfoil
	pts::AbstractMatrix
	colpts::AbstractMatrix
	n⃗::AbstractMatrix
	n̂::AbstractMatrix
	s⃗::AbstractVector
	sV::AbstractVector
	A::AbstractMatrix
	Ainv::AbstractMatrix
end

function A_gen(pts::AbstractMatrix, n̂::AbstractMatrix)
	
	Fg = eltype(pts)

	ndofs = size(pts, 1) - 1
	
	A = zeros(Fg, ndofs + 1, ndofs + 1)
	
	for i = 1:ndofs
		for j = 1:ndofs
			iptn = size(pts, 1) - j
			ipt = size(pts, 1) - j + 1
			
			panel_points = [
				pts[ipt, 1] pts[iptn, 1];
				pts[ipt, 2] pts[iptn, 2]
			]
			
			iptn = size(pts, 1) - i
			ipt = size(pts, 1) - i + 1
			
			collocation = [
				(pts[ipt, 1] + pts[iptn, 1]) / 2,
				(pts[ipt, 2] + pts[iptn, 2]) / 2
			]
			
			is_same = (i == j)
			
			u, w = influence(
				panel_points,
				collocation,
				[1.0, 0.0],
				is_same
			)
			
			A[i, j] += u * n̂[i, 1] + w * n̂[i, 2]
			
			u, w = influence(
				panel_points,
				collocation,
				[0.0, 1.0],
				is_same
			)

			A[i, j + 1] += u * n̂[i, 1] + w * n̂[i, 2]
		end
	end
	
	A[end, 1] = 1.0
	A[end, end] = 1.0
	
	A, inv(A)
	
end

function Airfoil(
		pts::AbstractMatrix
	)
	
	Fg = eltype(pts)

	colpts = pts[1:(end - 1), :]
	
	for i = 1:size(colpts, 1)
		colpts[i, :] .= (pts[i + 1, :] .+ pts[i, :]) ./ 2
	end
	
	n⃗ = vcat(
		[
			begin
				inxt = i % size(pts, 1) + 1
				
				x⃗ = pts[i, :] .- pts[inxt, :]
				
				[x⃗[2] (- x⃗[1])]
			end for i = (size(pts, 1) - 1):-1:1
		]...
	)
	
	
	nns = sqrt.(n⃗[:, 1] .^ 2 .+ n⃗[:, 2] .^ 2)
	
	n̂ = hcat(n⃗[:, 1] ./ nns, n⃗[:, 2] ./ nns)
	
	A, Ainv = A_gen(pts, n̂)
	
	currs = 0.0
	s⃗ = [
		begin
			inxt = i + 1
			
			L = sqrt((pts[inxt, 1] - pts[i, 1]) ^ 2 + 
				(pts[inxt, 2] - pts[i, 2]) ^ 2)
			
			rets = currs + L / 2
			currs += L
			
			rets
		end for i = 1:(size(pts, 1) - 1)
	]

	temp = 0.0
	sV = [
		0.0;
		[
			begin
				L = norm(pts[i + 1, :] .- pts[i, :])

				temp += L / 2

				s = temp

				temp += L / 2

				s
			end for i = 1:(size(pts, 1) - 1)
		]
	]
	
	Airfoil(pts, colpts, n⃗, n̂, s⃗, sV, A, Ainv)
	
end
# Airfoil(pts::Matrix{Fg}; kwargs...) where Fg = Airfoil{Fg}(pts; kwargs...)
Airfoil(filename::String; kwargs...) = Airfoil(
	readdlm(filename; skipstart = 1); kwargs...
)

function solve_inviscid(
		afl::Airfoil; 
		α::Real = 0.0
	)
	
	x̂ = [cosd(α), sind(α)]
	
	b = - x̂[1] * afl.n̂[:, 1] - x̂[2] * afl.n̂[:, 2]
	
	# push!(b, 0.0)
	b = [b; 0.0]

	γs = afl.Ainv * b

	γx = (γs[2] - γs[3]) / (afl.sV[2] - afl.sV[3])
	γ1 = γs[2] + γx * (afl.sV[1] - afl.sV[2])

	γx = (γs[end - 1] - γs[end - 2]) / (afl.sV[end - 1] - afl.sV[end - 2])
	γe = γs[end - 1] + γx * (afl.sV[end] - afl.sV[end - 1])

	γm = (γ1 - γe) / 2

	γs = [
		γm;
		(- γs[(end - 1):-1:2]);
		(- γm)
	]
	
end
