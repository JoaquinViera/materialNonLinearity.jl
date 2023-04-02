module materialNonLinearity

# ===============================================

# Dependencies
include("deps.jl")

# Structs definition
include("init.jl")

# Export functions
include("exports.jl")

# Utils functions
include("Utils/utils.jl")

# Element functions
# -----------------------------------------------
# Frame element
include("Frame/frame_element.jl")

# Solver 
# -----------------------------------------------
# Internal variables 
include("Solver/initial_defs.jl")
# Main function
include("Solver/solver.jl")
# Tanget stiffness matrix and internal forces computation
include("Solver/finte_KT_int.jl")
# Assembler
include("Solver/assembler.jl")
# Material models
include("Solver/constitutive_model.jl")
# Store solution
include("Solver/store_sol.jl")

# Numerical methods - algorithms
# -----------------------------------------------
# Newton Raphson
include("Algorithms/NR.jl")
# Arc Length
include("Algorithms/AL_DD.jl")
include("Algorithms/AL_Cylindrical.jl")
# Convergence check
include("Algorithms/convergence_check.jl")
# Compute Fext
include("Algorithms/compute_Fext.jl")

end # module
