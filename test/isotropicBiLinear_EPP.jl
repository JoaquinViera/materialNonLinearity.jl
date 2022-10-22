# ===============================================
# Cantilever with elasto plastic material model 
# ===============================================

# Load solver module
using materialNonLinearity, LinearAlgebra

# example name
problemName = "isotropicBiLinear_EPP"


# Define material model
# =======================================
E = 210e6
σY0 = 250e3
K = 0

# Materials struct
StrMaterialModels = IsotropicBiLinear(E, σY0, K)
# Define section
# =======================================
b = 0.1
h = 0.1
secName = "rectangle"
secParams = [b, h]

# Section struct
StrSections = Section(secName, secParams)

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
supps = [1 Inf Inf]

# Define applied external loads
Fy = -1
Mz = 0
nod = nnodes
nodalForces = [nod Fy Mz]

# BoundaryConds struct
StrBoundaryConds = BoundaryConds(supps, nodalForces)

# Numerical method parameters
# =======================================

tolk = 50 # number of iters
tolu = 1e-4 # Tolerance of converged disps
tolf = 1e-6 # Tolerance of internal forces
nLoadStε = 63 # Number of load increments
loadFactorsVec = ones(nLoadStε) # Load scaling factors

# Numerical method settings struct
StrAnalysisSettings = AnalysisSettings(tolk, tolu, tolf, loadFactorsVec)

# Plot parameters
# =======================================
lw = 3
ms = 2
color = "black"

strPlots = PlotSettings(lw, ms, color)

# ===============================================
# Process model parameters
# ===============================================

sol, time, IterData = solver(StrSections, StrMaterialModels, StrMesh, StrBoundaryConds, StrAnalysisSettings)

# Auxiliar
# --------------------------------
P = abs(Fy)
Iy = StrSections.Iy
κe = 2 * σY0 / (E * h)

# Numerical solution
matFint = sol.matFint
matUk = sol.matUk

nod = 1
dofM = nod * 2

mVec = matFint[dofM, :]

# Computes curvatures
# --------------------------------

kappaHistElem = zeros(nelems, nLoadStε)

rotXYXZ = Diagonal(ones(4, 4))
rotXYXZ[2, 2] = -1
rotXYXZ[4, 4] = -1

for j in 1:nelems
    nodeselem = StrMesh.conecMat[j, 3]
    elemdofs = nodes2dofs(nodeselem[:], 2)
    local R, l = element_geometry(StrMesh.nodesMat[nodeselem[1], :], StrMesh.nodesMat[nodeselem[2], :], 2)
    UkeL = R' * matUk[elemdofs, 1:end]
    Be = intern_function(0, l) * rotXYXZ
    kappaelem = Be * UkeL
    kappaHistElem[j, :] = abs.(kappaelem)
end

# Analytical solution M-κ
# --------------------------------

Mana = zeros(nLoadStε)
C = E * K / (E + K)
εY = σY0 / E
ε⃰ = εY - σY0 / C
κ⃰ = 2 * ε⃰ / h
elem = 1
for i in 1:nLoadStε
    κₖ = kappaHistElem[elem, i]
    if κₖ <= κe
        Mana[i] = E * StrSections.Iy * κₖ
    else
        Mana[i] = σY0 * b * h^2 / 12 * (3 - κe^2 / κₖ^2 + κₖ / κe * C / E * (2 - 3 * κe / κₖ + κe^3 / κₖ^3))
    end
end

@test (maximum(abs.(abs.(mVec[2:end]) - Mana[2:end]) ./ Mana[2:end])) <= 1e-2
