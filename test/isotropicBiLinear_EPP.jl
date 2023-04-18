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
K = E / 10
ne = 12
ns = 12

# Materials struct
StrMaterialModels = IsotropicBiLinear(E, σY0, K, ne, ns)

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
nnodes = 21
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

# Stress Array
# =======================================
elems = []
xG_Rel_Ind = zeros(ne)

StrStressArray = StressArraySets(elems, xG_Rel_Ind)


# =======================================
# NR test
# =======================================

# Numerical method parameters
# =======================================

tolk = 15 # number of iters
tolu = 1e-4 # Tolerance of converged disps
tolf = 1e-6 # Tolerance of internal forces
nLoadSteps = 63 # Number of load increments
loadFactorsVec = collect(1:nLoadSteps) # Load scaling factors

# Numerical method settings struct
StrAnalysisSettings = NewtonRaphson(tolk, tolu, tolf, loadFactorsVec)

# ===============================================
# Process model parameters
# ===============================================

sol, time, IterData = solver(StrSections, StrMaterialModels, StrMesh, StrBoundaryConds, StrAnalysisSettings, problemName, StrStressArray)


# Auxiliar
# --------------------------------
P = abs(Fz)
Iy = StrSections.Iy
κe = 2 * σY0 / (E * h)

# Numerical solution
matFint = sol.matFint
matUk = sol.matUk

elem = 1
dofM = 3
mVec = hcat([i[dofM] for i in matFint[elem]])

# Computes curvatures
# --------------------------------
xrel = zeros(nelems)
kappaHistElem = frame_curvature(nelems, StrMesh, nLoadSteps, matUk, xrel)

# Analytical solution M-κ
# --------------------------------

Mana = zeros(nLoadSteps)
C = E * K / (E + K)
εY = σY0 / E
ε⃰ = εY - σY0 / C
κ⃰ = 2 * ε⃰ / h
for i in 1:nLoadSteps
    κₖ = abs(kappaHistElem[elem, i])
    if κₖ <= κe
        Mana[i] = E * StrSections.Iy * κₖ
    else
        Mana[i] = σY0 * b * h^2 / 12 * (3 - κe^2 / κₖ^2 + κₖ / κe * C / E * (2 - 3 * κe / κₖ + κe^3 / κₖ^3))
    end
end

# @test (maximum(abs.(abs.(mVec[2:end]) - Mana[2:end]) ./ Mana[2:end])) <= 1e-2

# =======================================
# AL test - dominant dof
# =======================================

# Numerical method parameters
# =======================================

# tolk = 50 # number of iters
# tolu = 1e-4 # Tolerance of converged disps
# tolf = 1e-6 # Tolerance of internal forces
# initialDeltaLambda = 1e-3 #
# arcLengthIncrem = vcat(ones(40) * 1e-3) #
# nLoadSteps = length(arcLengthIncrem) # Number of load increments
# controlDofs = [6] #
# scalingProjection = 1 #

# # Numerical method settings struct
# StrAnalysisSettings = ArcLength(tolk, tolu, tolf, nLoadSteps, initialDeltaLambda, arcLengthIncrem, controlDofs, scalingProjection)

# # ===============================================
# # Process model parameters
# # ===============================================

# sol, time, IterData = solver(StrSections, StrMaterialModels, StrMesh, StrBoundaryConds, StrAnalysisSettings, problemName, StrStressArray)

# # Auxiliar
# # --------------------------------
# P = abs(Fz)
# Iy = StrSections.Iy
# κe = 2 * σY0 / (E * h)

# # Numerical solution
# matFint = sol.matFint
# matUk = sol.matUk

# elem = 1
# dofM = 3
# mVec = hcat([i[dofM] for i in matFint[elem]])

# # Computes curvatures
# # --------------------------------
# xrel = zeros(nelems)
# kappaHistElem = frame_curvature(nelems, StrMesh, nLoadSteps, matUk, xrel)

# # Analytical solution M-κ
# # --------------------------------

# Mana = zeros(nLoadSteps)
# C = E * K / (E + K)
# εY = σY0 / E
# ε⃰ = εY - σY0 / C
# κ⃰ = 2 * ε⃰ / h
# for i in 1:nLoadSteps
#     κₖ = abs(kappaHistElem[elem, i])
#     if κₖ <= κe
#         Mana[i] = E * StrSections.Iy * κₖ
#     else
#         Mana[i] = σY0 * b * h^2 / 12 * (3 - κe^2 / κₖ^2 + κₖ / κe * C / E * (2 - 3 * κe / κₖ + κe^3 / κₖ^3))
#     end
# end

# @test (maximum(abs.(abs.(mVec[2:end]) - Mana[2:end]) ./ Mana[2:end])) <= 1e-2


# =======================================
# AL test - cylindrical constraint
# =======================================

# Numerical method parameters
# =======================================

tolk = 50 # number of iters
tolu = 1e-8 # Tolerance of converged disps
tolf = 1e-6 # Tolerance of internal forces
initialDeltaLambda = 1e-3 #
arcLengthIncrem = vcat(ones(45) * 1e-3) #
arcLengthIncrem = vcat(ones(20) * 1e-3, ones(60) * 5e-4) #
# arcLengthIncrem = vcat(ones(70) * 1e-3, ones(5) * 5e-4) #
nLoadSteps = length(arcLengthIncrem) # Number of load increments
controlDofs = [6] #

# Numerical method settings struct
StrAnalysisSettings = ArcLength_Cylindrical(tolk, tolu, tolf, nLoadSteps, initialDeltaLambda, arcLengthIncrem, controlDofs)

# ===============================================
# Process model parameters
# ===============================================

sol, time, IterData = solver(StrSections, StrMaterialModels, StrMesh, StrBoundaryConds, StrAnalysisSettings, problemName, StrStressArray)

# Auxiliar
# --------------------------------
P = abs(Fz)
Iy = StrSections.Iy
κe = 2 * σY0 / (E * h)

# Numerical solution
matFint = sol.matFint
matUk = sol.matUk

elem = 1
dofM = 3
mVec = hcat([i[dofM] for i in matFint[elem]])

# Computes curvatures
# --------------------------------
xrel = zeros(nelems)
kappaHistElem = frame_curvature(nelems, StrMesh, nLoadSteps, matUk, xrel)

# Analytical solution M-κ
# --------------------------------

Mana = zeros(nLoadSteps)
C = E * K / (E + K)
εY = σY0 / E
ε⃰ = εY - σY0 / C
κ⃰ = 2 * ε⃰ / h
for i in 1:nLoadSteps
    κₖ = abs(kappaHistElem[elem, i])
    if κₖ <= κe
        Mana[i] = E * StrSections.Iy * κₖ
    else
        Mana[i] = σY0 * b * h^2 / 12 * (3 - κe^2 / κₖ^2 + κₖ / κe * C / E * (2 - 3 * κe / κₖ + κe^3 / κₖ^3))
    end
end

fig = plot(abs.(kappaHistElem[elem, :]), abs.(mVec), markershape=:circle, legend=:false)

@test (maximum(abs.(abs.(mVec[2:end]) - Mana[2:end]) ./ Mana[2:end])) <= 1e-2

println("All tests passed for problem: $problemName !")

fig