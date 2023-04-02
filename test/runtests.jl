using Test, materialNonLinearity

# Checks bending moment, displacement, rotation and stiffness matrix
include("linearElastic.jl")

# Checks bending moment with NR and AL algortihms with Elastic Plastic Perfectly material
include("isotropicBiLinear_EPP.jl")

# User defined constitutive model
include("userModel.jl")

println("All tests passed!")