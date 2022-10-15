# ===============================================
# Cantilever with elasto plastic material model 
# ===============================================

# Load solver module
using materialNonLinearity, Plots

# example name
problemName = "CantileverEPP"


# Define material model
# =======================================
E = 210e6
σY = 250e3
K = E / 20
matName = "isotropicBiLinear"
matParams = [E, σY, K]

# Materials struct
strMaterialModels = materialModel(matName, matParams)

# Define section
# =======================================
b = 0.1
h = 0.1
secName = "rectangle"
secParams = [b, h]

# Section struct
strSections = section(secName, secParams)

# Define Mesh
# =======================================

# Nodes
L = 1
nnodes = 101
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
strMesh = mesh(Nodes, Conec)

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
strBC = boundaryConds(supps, nodalForces)

# Numerical method parameters
# =======================================

tolk = 50 # number of iters
tolu = 1e-4 # Tolerance of converged disps
tolf = 1e-6 # Tolerance of internal forces
nLoadSteps = 80 # Number of load increments
loadFactorsVec = ones(nLoadSteps) # Load scaling factors

# Numerical method settings struct
strAnalysisSets = analysisSettings(tolk, tolu, tolf, nLoadSteps, loadFactorsVec)

# Plot parameters
# =======================================
lw = 3
ms = 2
color = "black"

strPlots = plotSettings(lw, ms, color)

# ===============================================
# Process model parameters
# ===============================================

sol, time = solver(strSections, strMaterialModels, strMesh, strBC, strAnalysisSets)



# Check First step
# --------------------------------
P = abs(Fy)
Iy = strSections.Iy
# Analytical solution
Man = P * L
fan = P * L^3 / (3 * E * Iy)
theta_an = P * L^2 / (2 * E * Iy)
kappae = 2 * σY / (E * h)

# Numerical solution
matFint = sol.matFint
matUk = sol.matUk
nod = 1

dofM = nod * 2
dofD = nnodes * 2 - 1
dofT = nnodes * 2

mVec = matFint[dofM, :]

Mnum = mVec[2]
deltaNum = matUk[dofD, 2]
thetaNum = matUk[dofT, 2]

# Computes curvatures
# --------------------------------

kappaHistElem = zeros(nelems, nLoadSteps)

for j in 1:nelems
    nodeselem = strMesh.conecMat[j, 3]
    elemdofs = nodes2dofs(nodeselem[:], 2)
    local R, l = elemGeom(strMesh.nodesMat[nodeselem[1], :], strMesh.nodesMat[nodeselem[2], :], 2)
    UkeL = R' * matUk[elemdofs, 1:end]
    Be = internFunction(0, l)
    kappaelem = Be * UkeL
    kappaHistElem[j, :] = abs.(kappaelem)
end

# Analytical solution M-κ
# --------------------------------

Mana = zeros(nLoadSteps)
C = E * K / (E + K)
elem = 1
for i in 1:nLoadSteps
    kappak = kappaHistElem[elem, i]
    if kappak <= kappae
        Mana[i] = E * strSections.Iy * kappak
    else
        Mana[i] = σY * b * h^2 / 12 * (3 - kappae^2 / kappak^2 + kappak / kappae * C / E * (2 - 3 * kappae / kappak + kappae^3 / kappak^3))
    end
end

fig = plot(kappaHistElem[elem, :], Mana, markershape=:circle, lw=lw, ms=ms)
plot!(fig, kappaHistElem[elem, :], mVec, markershape=:rect, lw=lw, ms=ms)



# Check KTe
Uke = zeros(4)
l = 1
Iy = strSections.Iy
E = strMaterialModels.E

Finte, KTe = finte_KT_int(strMaterialModels, l, strSections.params, Uke, 1)
Kana = E * Iy / l^3 * [12 6l -12 6l; 6l 4l^2 -6l 2l^2; -12 -6l 12 -6l; 6l 2l^2 -6l 4l^2]