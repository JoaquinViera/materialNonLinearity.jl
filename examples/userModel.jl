# ===============================================
# Cantilever with elasto plastic material model 
# ===============================================

# Load solver module
using materialNonLinearity, Plots, LinearAlgebra

# example name
problemName = "CantileverEPP"

# Define material model
# =======================================
E = 210e6
σY = 250e3
K = -E / 100

# Materials struct
StrMaterialModels = UserModel()

# Define section
# =======================================
b = 0.3
h = 0.3

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
supps = [1 Inf Inf]

# Define applied external loads
Fy = -0.1
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
nLoadSteps = 40 # Number of load increments
initialDeltaLambda = 1e-5 #
arcLengthIncrem = 1e-5 #
controlDofs = [10] #
scalingProjection = 1 #

# Numerical method settings struct
StrAnalysisSettings = ArcLength(tolk, tolu, tolf, nLoadSteps, initialDeltaLambda, arcLengthIncrem, controlDofs, scalingProjection)

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

# Post process
# --------------------------------
P = abs(Fy)
Iy = StrSections.Iy

# Numerical solution
matFint = sol.matFint
matUk = sol.matUk

# Clamped node
nod = 1
dofM = nod * 2

# Loaded node
dofD = nnodes * 2 - 1
dofT = nnodes * 2

# Reaction Bending moment 
mVec = matFint[dofM, :]
Mnum = mVec[2]

# Displacements at loaded node
deltaNum = matUk[dofD, 2]
thetaNum = matUk[dofT, 2]

dVec = abs.(matUk[dofD, :])

pVec = abs.(mVec / L)

figPdelta = plot(dVec, pVec, markershape=:circle, lw=lw, ms=ms, title="P-δ", label="FEM", minorgrid=1, draw_arrow=1)
xlabel!("δ")
ylabel!("M")

# Compute curvatures
# --------------------------------

kappaHistElem = zeros(nelems, nLoadSteps)

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

# M-κ plot  
# --------------------------------
elem = 1
fig = plot(kappaHistElem[elem, :], abs.(mVec), markershape=:circle, lw=lw, ms=ms, title="M-κ", label="FEM", minorgrid=1, draw_arrow=1)
xlabel!("κ")
ylabel!("M")

# P-δ plot  
# --------------------------------
figPdelta = plot(dVec, pVec, markershape=:circle, lw=lw, ms=ms, title="P-δ", label="FEM", minorgrid=1, draw_arrow=1)
xlabel!("δ")
ylabel!("P")
