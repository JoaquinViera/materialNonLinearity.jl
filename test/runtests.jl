using Test, materialNonLinearity

@test hello("Julia") == "Hello, Julia"
@test domath(2.0) ≈ 7.0

include("linearElastic.jl")