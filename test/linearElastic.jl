# ===============================================
# Cantilever with linear elastic material model 
# ===============================================

# Load solver module
using materialNonLinearity, LinearAlgebra

# example name
problemName = "linearElastic"

# Define material model
# =======================================
E = 210e6
ne = 2
ns = 2

# Materials struct
StrMaterialModels = LinearElastic(E, ne, ns)

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
nnodes = 11
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
Fz = -1
My = 0
nod = nnodes
nodalForces = [nod Fz My]

# BoundaryConds struct
StrBoundaryConds = BoundaryConds(supps, nodalForces)

# Numerical method parameters
# =======================================

tolk = 50 # number of iters
tolu = 1e-4 # Tolerance of converged disps
tolf = 1e-6 # Tolerance of internal forces
nLoadSteps = 3 # Number of load increments
loadFactorsVec = ones(nLoadSteps) # Load scaling factors

# Numerical method settings struct
StrAnalysisSettings = NewtonRaphson(tolk, tolu, tolf, loadFactorsVec)

# ===============================================
# Process model parameters
# ===============================================

sol, time, IterData = solver(StrSections, StrMaterialModels, StrMesh, StrBoundaryConds, StrAnalysisSettings, problemName)

# Check KTe
Iy = StrSections.Iy
Uke = zeros(4)
l = 1
rotXYXZ = Diagonal(ones(4, 4))
rotXYXZ[2, 2] = -1
rotXYXZ[4, 4] = -1

Finte, KTe = finte_KT_int(StrMaterialModels, l, [b, h], Uke, 1)
Kana = rotXYXZ * E * Iy / l^3 * [12 6l -12 6l; 6l 4l^2 -6l 2l^2; -12 -6l 12 -6l; 6l 2l^2 -6l 4l^2] * rotXYXZ

# Check First step
# --------------------------------
P = abs(Fz)
# Analytical solution
Man = -P * L
δan = -P * L^3 / (3 * E * Iy)
θan = P * L^2 / (2 * E * Iy)

# Numerical solution
matFint = sol.matFint
matUk = sol.matUk

nod = 1
dofM = nod * 2
dofD = nnodes * 2 - 1
dofT = nnodes * 2

mVec = matFint[dofM, :]

Mnum = mVec[2]
δNum = matUk[2][dofD]
θNum = matUk[2][dofT]

@test abs(Mnum - Man) <= tolf
@test abs(δNum - δan) <= tolu
@test abs(θNum - θan) <= tolu
@test norm(KTe - Kana) <= 1e-6