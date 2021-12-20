# ArkRoyal

Código para análise de aerofólios do mestrado do Frustrado (T-22).

Para instruções de utilização, veja os exemplos na pasta examples!

As funções principais (Airfoil() e solve()) têm os seguintes protótipos:

```

function Airfoil(
		pts::AbstractMatrix # pontos em form. Selig
	)

function solve(
	afl::Airfoil;
	α::Real = 0.0, # graus
	Re::Real = 1e6,
	M∞::Real = 0.0,
	Ncr::Real = 9.0,
	n_iter::Int64 = 5,
	stall_correction::Symbol = :Eppler,
	kwargs...
)

# solve retorna:
struct AirfoilSolution
	afl::Airfoil
	u::AbstractVector # distribuições pelos pontos
	Cₚ::AbstractVector
	θ::AbstractVector
	H::AbstractVector
	N::AbstractVector
	Cf::AbstractVector
	CL::Real
	CD::Real
	Cm::Real
end
```
