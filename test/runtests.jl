using Test, materialNonLinearity

# Checks bending moment, displacement, rotation and stiffness matrix
include("linearElastic.jl")

# Checks bending moment with NR and AL algortihms with Elastic PErfectly material
# todo add AL
include("isotropicBiLinear_EPP.jl")