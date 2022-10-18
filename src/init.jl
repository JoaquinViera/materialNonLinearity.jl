# Declaration

# ============================================================================
# Material 
# ============================================================================
struct MaterialModel{Real}
    name::String
    params::Vector{Real}
    E::Real
    σY0::Real
    K::Real
end
# Constructor
function MaterialModel(name, params::Vector)
    if cmp(name, "linearElastic") == 1 && cmp(name, "isotropicBiLinear") == 1
        error("Current models are linearElastic & isotropicBiLinear.")
    else
        if cmp(name, "linearElastic") == 0
            length(params) != 1 ? error("Define only Young Modulus.") : nothing
            E = params[1]
            σY0 = 0.0
            K = 0.0
        elseif cmp(name, "isotropicBiLinear") == 0
            length(params) != 2 && length(params) != 3 ? error("Define at least 2 parameters.") : nothing
            E = params[1]
            σY0 = params[2]
            length(params) == 2 ? K = 0.0 : K = params[3]
        end
        return MaterialModel(name, params, E, σY0, K)
    end
end

#=
matName = "linearElastic"
matName2 = "isotropicBiLinear"
matParams = [10.0]
matParams2 = [10.0, 1.0, 7]
mat1 = MaterialModel(matName, matParams)
mat2 = MaterialModel(matName2, matParams2)
=#

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