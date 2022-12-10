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
nnodes = 3
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

# Define Supports [ node ux uz thetay ]
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

tolk = 3 # number of iters
tolu = 1e-4 # Tolerance of converged disps
tolf = 1e-6 # Tolerance of internal forces
nLoadSteps = 10 # Number of load increments
loadFactorsVec = collect(1:nLoadSteps) # Load scaling factors

# Numerical method settings struct
StrAnalysisSettings = NewtonRaphson(tolk, tolu, tolf, loadFactorsVec)

# Stress Array
# =======================================
elems = []
xG_Rel_Ind = 0

StrStressArray = StressArraySets(elems, xG_Rel_Ind)

# ===============================================
# Process model parameters
# ===============================================

sol, time, IterData = solver(StrSections, StrMaterialModels, StrMesh, StrBoundaryConds, StrAnalysisSettings, problemName, StrStressArray)

# Check 
P = abs(Fz)
Iy = StrSections.Iy
A = StrSections.A
l = 1

# Check First step
# --------------------------------

# Analytical solution
Man = -P * L
δan = -P * L^3 / (3 * E * Iy)
θan = P * L^2 / (2 * E * Iy)

# Numerical solution
matFint = sol.matFint
matUk = sol.matUk

nod = 1
dofM = nod * 3
dofD = nnodes * 3 - 1
dofT = nnodes * 3

mVec = matFint[dofM, :]

Mnum = mVec[2]
δNum = matUk[2][dofD]
θNum = matUk[2][dofT]

# Check stiffness matrix
Uke = zeros(6)
rotXYXZ = Diagonal(ones(4, 4))
rotXYXZ[2, 2] = -1
rotXYXZ[4, 4] = -1

Finteb, Fintea, σArr, KTeb, KTea = finte_KT_int(StrMaterialModels, l, [b, h], Uke, 1, [], 1, 1, StrStressArray)

Kbending = rotXYXZ * E * Iy / l^3 * [12 6l -12 6l; 6l 4l^2 -6l 2l^2; -12 -6l 12 -6l; 6l 2l^2 -6l 4l^2] * rotXYXZ
Kaxial = E * A / l * [1 -1; -1 1]

@test abs(Mnum - Man) <= tolf
@test abs(δNum - δan) <= tolu
@test abs(θNum - θan) <= tolu
@test norm(KTeb - Kbending) <= 1e-6 * norm(Kbending)
@test norm(KTea - Kaxial) <= 1e-6 * norm(Kaxial)

# =======================================
# AL test
# =======================================

# Numerical method parameters
# =======================================

tolk = 15 # number of iters
tolu = 1e-10 # Tolerance of converged disps
tolf = 1e-6 # Tolerance of internal forces
initialDeltaLambda = 2e-3 #
arcLengthIncrem = vcat(ones(10) * 1e-4) #
nLoadSteps = length(arcLengthIncrem) # Number of load increments
controlDofs = [6] #
scalingProjection = 1 #

StrAnalysisSettings = ArcLength(tolk, tolu, tolf, nLoadSteps, initialDeltaLambda, arcLengthIncrem, controlDofs, scalingProjection)


# ===============================================
# Process model parameters
# ===============================================

sol, time, IterData = solver(StrSections, StrMaterialModels, StrMesh, StrBoundaryConds, StrAnalysisSettings, problemName, StrStressArray)

# Check 
P = abs(Fz) * abs(sol.loadFactors[2])
Iy = StrSections.Iy
A = StrSections.A
l = 1

# Check First step
# --------------------------------

# Analytical solution
Man = -P * L
δan = -P * L^3 / (3 * E * Iy)
θan = P * L^2 / (2 * E * Iy)

matFint = sol.matFint
matUk = sol.matUk

nod = 1
dofM = nod * 3
dofD = nnodes * 3 - 1
dofT = nnodes * 3

mVec = matFint[dofM, :]

Mnum = mVec[2]
δNum = matUk[2][dofD]
θNum = matUk[2][dofT]

# Check stiffness matrix
matUk = sol.matUk
elem = nnodes - 1
nodeselem = Conec[elem, 3]
dofselem = nodes2dofs(nodeselem, 3)
Uke = matUk[end][dofselem]

rotXYXZ = Diagonal(ones(4, 4))
rotXYXZ[2, 2] = -1
rotXYXZ[4, 4] = -1

Finteb, Fintea, σArr, KTeb, KTea = finte_KT_int(StrMaterialModels, l, [b, h], Uke, 1, [], 1, 1, StrStressArray)

Kbending = rotXYXZ * E * Iy / l^3 * [12 6l -12 6l; 6l 4l^2 -6l 2l^2; -12 -6l 12 -6l; 6l 2l^2 -6l 4l^2] * rotXYXZ
Kaxial = E * A / l * [1 -1; -1 1]

@test abs(Mnum - Man) <= tolf
@test abs(δNum - δan) <= tolu
@test abs(θNum - θan) <= tolu
@test norm(KTeb - Kbending) <= 1e-6 * norm(Kbending)
@test norm(KTea - Kaxial) <= 1e-6 * norm(Kaxial)

