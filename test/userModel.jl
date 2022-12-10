# ===============================================
# Cantilever with user defined function 
# ===============================================

# Load solver module
using materialNonLinearity, LinearAlgebra

# example name
problemName = "cubicFunction"

# Define material model
# =======================================
E = 210e6
σY = 250e3
ne = 16
ns = 16

import materialNonLinearity: constitutive_model

# Materials struct
StrMaterialModels = UserModel(ne, ns)

function constitutive_model(ElemMaterialModel::UserModel, εₖ)

    E = 210e6
    σY = 250e3
    εY = σY / E

    a = -σY / (2 * εY^3)
    b = 3 * σY / (2 * εY)

    if abs(εₖ) >= sqrt(-b / a)
        σ = 0
        ∂σ∂ε = 0
    else
        σ = a * εₖ^3 + b * εₖ
        ∂σ∂ε = 3 * a * εₖ^2 + b
    end

    return σ, ∂σ∂ε

end


# Define section
# =======================================
b = 0.1
h = 0.1

# Section struct
StrSections = Rectangle(; b, h)

# Define Mesh
# =======================================

# Nodes
L = 1
nnodes = 31
xcoords = collect(LinRange(0, L, nnodes))
ycoords = zeros(length(xcoords))
Nodes = hcat(xcoords, ycoords)

# Conec
elemConec = []
for i in 1:(nnodes-1)
    global elemConec = vcat(elemConec, (i, i + 1))
end

nelems = nnodes - 1
matVec = ones(nelems)
secVec = ones(nelems)
Conec = hcat(matVec, secVec, elemConec)

# Mesh struct
StrMesh = Mesh(Nodes, Conec)

# Boundary conditions
# =======================================

# Define Supports
supps = [1 Inf Inf Inf]

# Define applied external loads
Fx = 0
Fz = -1
My = 0
nod = nnodes
nodalForces = [nod Fx Fz My]

# BoundaryConds struct
StrBoundaryConds = BoundaryConds(supps, nodalForces)

# Numerical method parameters
# =======================================

tolk = 50 # number of iters
tolu = 1e-4 # Tolerance of converged disps
tolf = 1e-6 # Tolerance of internal forces
initialDeltaLambda = 1e-3 #
arcLengthIncrem = vcat(ones(21) * 4e-4) #
nLoadSteps = length(arcLengthIncrem) # Number of load increments
controlDofs = [nod * 3 - 1] #
scalingProjection = 1 #

# Numerical method settings struct
StrAnalysisSettings = ArcLength(tolk, tolu, tolf, nLoadSteps, initialDeltaLambda, arcLengthIncrem, controlDofs, scalingProjection)

# Plot parameters
# =======================================
lw = 3
ms = 2
color = "black"

strPlots = PlotSettings(lw, ms, color)

# Stress Array
# =======================================
elems = []
xG_Rel_Ind = 0

StrStressArray = StressArraySets(elems, xG_Rel_Ind)

# ===============================================
# Process model parameters
# ===============================================

sol, time, IterData = solver(StrSections, StrMaterialModels, StrMesh, StrBoundaryConds, StrAnalysisSettings, problemName, StrStressArray)

# Post process
# --------------------------------
P = abs(Fz)
Iy = StrSections.Iy

# Numerical solution
matFint = sol.matFint
matUk = sol.matUk

# Clamped node
nod = 1
dofM = nod * 3

# Reaction Bending moment 
mVec = matFint[dofM, :]

# Compute curvatures
# --------------------------------
kappaHistElem = frame_curvature(nelems, StrMesh, nLoadSteps, matUk)

# Analytical solution M-κ
# --------------------------------

Mana = zeros(nLoadSteps)
εY = σY / E
ca = -σY / (2 * εY^3);
cb = 3 * σY / (2 * εY);

elem = 1
for i in 1:nLoadSteps
    κₖ = abs(kappaHistElem[elem, i])
    Mana[i] = κₖ * b * (ca * κₖ^2 * h^5 / 80 + cb * h^3 / 12)
end

err = (abs.(mVec[2:end]) - Mana[2:end]) ./ Mana[2:end] * 100

@test (maximum(abs.(abs.(mVec[2:end]) - Mana[2:end]) ./ Mana[2:end])) <= 1e-2
