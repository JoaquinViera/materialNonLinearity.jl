# Declaration

# ============================================================================
# Material 
# ============================================================================

abstract type MaterialModel end

struct LinearElastic <: MaterialModel
    E::Real
end

function LinearElastic(; E::Real)
    LinearElastic(E)
end

struct IsotropicBiLinear <: MaterialModel
    E::Real
    σY0::Real
    K::Real
end

function IsotropicBiLinear(; E::Real, σY0::Real, K::Real)
    IsotropicBiLinear(E, σY0, K)
end

# ============================================================================
# Section
# ============================================================================

struct Section{Real}
    name::String
    params::Vector{Real}
    A::Float64
    Iy::Float64
end

function Section(name, params)
    if cmp(name, "rectangle") == 0
        length(params) != 2 ? error("Declare only width and height.") : nothing
        A = params[1] * params[2]
        Iy = params[1] * params[2]^3 / 12
    else
        error("To be implemented.")
    end
    return Section(name, params, A, Iy)
end

#=
secName = "rectangle"
secParams = [0.1, 0.6]
secStr = Section(secName, secParams)
=#

# ============================================================================
# Mesh
# ============================================================================

struct Mesh
    nodesMat::Matrix{Float64}
    conecMat::Matrix
end

#=
a = [1, 2, 3, 4]
b = [4, 3, 2, 1]
meshNodes = Nodes(a, b)
=#

# ============================================================================
# Boundary conditions
# ============================================================================

struct BoundaryConds
    suppMatrix::Matrix
    nodalForceMatrix::Matrix
end

# ============================================================================
# Analysis settings
# ============================================================================

struct AnalysisSettings
    tolk::Float64
    tolu::Float64
    tolf::Float64
    loadFactors::Vector{Float64}
end


# ============================================================================
# Plots settings
# ============================================================================

struct PlotSettings
    lw::Int64
    ms::Int64
    color::String
end


# ============================================================================
# Store Solution
# ============================================================================

mutable struct ModelSol
    Uk::Vector{Float64}
    convδu::Array{Float64}
    Fextk::Vector{Float64}
    Fintk::Vector{Float64}
    matUk::Array{Float64}
    matFext::Array{Float64}
    matFint::Array{Float64}
    freeDofs::Vector{Int64}
    loadFactors::Array{Float64}
end

# ============================================================================
# Iter parameters
# ============================================================================

mutable struct IterParams
    nTimes::Int64
    stopCrit::Array{Int64}
end