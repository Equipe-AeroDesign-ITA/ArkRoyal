# ArkRoyal

My code for airfoil analysis using integral boundary layer methods.

It employs:

* Coles' relationship for turbulent closure, for skin friction;
* Thwaites' method for laminar closure and laminar skin friction;
* Coenen/Head/Veldman's entrainment equation for turbulent closure;
* Drela's e-N envelope method for transition prediction.

```

function Airfoil(
		pts::AbstractMatrix # points in Selig format
	)

function Airfoil(
                fname::String # file name, with first line as the airfoil name
	)

function solve(
	afl::Airfoil;
	α::Real = 0.0, # degrees
	Re::Real = 1e6,
	M∞::Real = 0.0,
	Ncr::Real = 9.0, # as per Drela's e-N envelope method
	N_bypass::Real = 0.0, # initial value for bypass turbulence
	n_iter::Int64 = 50,
	tol::Real = 0.0, # tolerance for pointwise Newton method
	ω::Real = 0.2, # relaxation factor for pointwise Newton method
	kwargs...
)

# solve returns:
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
	converged::Bool # indicates whether the pointwise Newton method has succeeded
end
```
