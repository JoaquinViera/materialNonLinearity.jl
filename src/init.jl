# Declaration

# ============================================================================
# Material 
# ============================================================================

abstract type MaterialModel end

struct LinearElastic <: MaterialModel
    E::Float64
    ne::Int64
    ns::Int64
end

function LinearElastic(; E::Float64, ne::Int64, ns::Int64)
    LinearElastic(E, ne, ns)
end

struct IsotropicBiLinear <: MaterialModel
    E::Float64
    σY0::Float64
    K::Float64
    ne::Int64
    ns::Int64
end

function IsotropicBiLinear(; E::Float64, σY0::Float64, K::Float64, ne::Int64, ns::Int64)
    IsotropicBiLinear(E, σY0, K, ne, ns)
end

struct UserModel <: MaterialModel
    ne::Int64
    ns::Int64
end

function UserModel(; ne::Int64, ns::Int64)
    UserModel(ne, ns)
end


# ============================================================================
# Section
# ============================================================================

abstract type Section end

mutable struct Rectangle <: Section
    b::Real
    h::Real
    A::Real
    Iy::Real
end

function Rectangle(; b, h)
    A = b * h
    Iy = b * h^3 / 12
    return Rectangle(b, h, A, Iy)
end

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
# Algorithm
# ============================================================================

abstract type AbstractAlgorithm end

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
    convδu::Vector{Vector{Float64}}
    Fextk::Vector{Float64}
    Fintk::Vector{Float64}
    matUk::Vector{Vector{Float64}}
    matFext::Vector{Vector{Float64}}
    matFint::Array{Float64}
    freeDofs::Vector{Int64}
    loadFactors::Vector{Float64}
end

# ============================================================================
# Iter parameters
# ============================================================================

mutable struct IterParams
    nTimes::Int64
    stopCrit::Vector{Int64}
end
