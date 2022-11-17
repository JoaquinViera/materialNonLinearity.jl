# Declaration

# ============================================================================
# Material 
# ============================================================================

"""
    MaterialModel 

Abstract type to define material constitutive model.

"""
abstract type MaterialModel end


"""
    LinearElastic

Isotropic Linear elastic constitutive model. \\
Stress is proportional to strain given by Elasticity Modulus `E`.

Fields:

* `E::Float64`: Young Elasticity modulus.  
* `ne::Int64`: Gauss integration points along the element length.
* `ns::Int64`: Gauss integration points along the section height. 

"""
struct LinearElastic <: MaterialModel
    E::Float64
    ne::Int64
    ns::Int64
end

function LinearElastic(; E::Float64, ne::Int64, ns::Int64)
    LinearElastic(E, ne, ns)
end

"""
    IsotropicBiLinear

Isotropic BiLinear constitutive model i.e. Linear softening or Linear Hardening. \\
Stress is proportional to strain given by either the Elasticity Modulus `E` and the Elastoplastic tangent modulus `C`.\\
Change in proportionality occurs after the stress reaches a given value `σY0`.

Fields:

* `E::Float64`: Young Elasticity modulus.  
* `σY0::Float64`: Initial yield stress.
* `K::Float64`: Slope of the hardening function.
* `ne::Int64`: Gauss integration points along the element length.
* `ns::Int64`: Gauss integration points along the section height. 

"""
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

"""
    UserModel

User defined constitutive model. \\
The user shall use ``` import materialNonLinearity: constitutive_model ``` and define the function `constitutive_model` as follows.\\


```julia
function constitutive_model(ElemMaterialModel::UserModel, εₖ)
    ...
    return σ, ∂σ∂ε
end
```


Fields:

* `ne::Int64`: Gauss integration points along the element length.
* `ns::Int64`: Gauss integration points along the section height. 

"""
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

"""
    Section

Abstract type to define element sections.

"""
abstract type Section end

"""
    Rectangle

Rectangular section.\\


Fields:

* `b::Float64`: Section width.
* `h::Float64`: Section height.
* `A::Float64`: Section area.
* `Iy::Float64`: Second moment of inertia. 

Note:\\
`A` and `Iy` are computed internally, only provide fields `b` and `h`.

"""
mutable struct Rectangle <: Section
    b::Float64
    h::Float64
    A::Float64
    Iy::Float64
end

function Rectangle(; b, h)
    A = b * h
    Iy = b * h^3 / 12
    return Rectangle(b, h, A, Iy)
end

# ============================================================================
# Mesh
# ============================================================================

"""
    Mesh

Contains the mesh node coordinates and conectivity matrix.\\


Fields:

* `nodesMat::Matrix{Float64}`: Matrix containing the (``` xᵢ,zᵢ```) coordinates of the i-th node. 
* `conecMat::Matrix`: Matrix whose i-th row contains the information of the i-th element according the structure: [```Material, Section, (nodeᵢ, nodeⱼ)```].


"""
struct Mesh
    nodesMat::Matrix{Float64}
    conecMat::Matrix
end

# ============================================================================
# Boundary conditions
# ============================================================================

"""
    BoundaryConds

Contains the boundary conditions of the structure.\\
Includes the support definition and the external forces.\\

Fields:

* `suppMatrix::Matrix`: Matrix whose i-th row contains the fixed degrees of freedom of the fixed i-th node according to the structure: [```Nodeᵢ uz θy```].  
* `nodalForceMatrix::Matrix`: Matrix whose i-th row contains the external nodal forces applied to the loaded i-th node according the structure: [```Nodeᵢ Fz My```].


"""
struct BoundaryConds
    suppMatrix::Matrix
    nodalForceMatrix::Matrix
end

# ============================================================================
# Algorithm
# ============================================================================

"""
    AbstractAlgorithm

Abstract type to define numerical method algorithms.

"""
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

"""
    ModelSol

Stores the solution of the problem.\\

"""
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

"""
    IterPArams

Store the criterion stop at each load increment.\\

"""
mutable struct IterParams
    nTimes::Int64
    stopCrit::Vector{Int64}
end
