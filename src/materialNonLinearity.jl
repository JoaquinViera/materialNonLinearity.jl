module materialNonLinearity

# ===============================================

export hello, domath

"""
    hello(who::String)
Return "Hello, `who`".
"""
hello(who::String) = "Hello, $who"

"""
    domath(x::Number)
Return `x + 5`.
"""
domath(x::Number) = x + 5


# ===============================================

# Dependencies
include("deps.jl")

# Structs definition
include("init.jl")

# Solver initialization
# -----------------------------------------------
# Internal variables 
include("ini_defs.jl")

# General Functions
include("nodes2dofs.jl")
include("elemGeom.jl")

# Main function
include("solver.jl")
include("assembler.jl")
include("convergenceCheck.jl")

# Numerical methods
include("NR.jl")

# Core functions
include("finte_KT_int.jl")
include("constitutiveModel.jl")


# Export functions
include("exports.jl")

end # module
