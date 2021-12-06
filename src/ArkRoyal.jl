module ArkRoyal

	using LinearAlgebra
	using ForwardDiff

	using DelimitedFiles
	
	include("VortexPanels.jl")

	include("Coenen.jl")
	include("Drela.jl")
	include("Thwaites.jl")

	include("Residuals.jl")
	include("Solve.jl")

	include("NACA.jl")

end # module
