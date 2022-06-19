module ArkRoyal

	using LinearAlgebra
	using ForwardDiff

	using DelimitedFiles

	using Interpolations
	
	include("VortexPanels.jl")

	include("Coenen.jl")
	include("Drela.jl")
	include("Thwaites.jl")

	include("Residuals.jl")
	include("Solve.jl")

	include("NACA.jl")
	include("Spline.jl")

end # module
