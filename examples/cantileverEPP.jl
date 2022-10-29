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
matName = "isotropicBiLinear"
matParams = [E, σY, K]

# Materials struct
StrMaterialModels = IsotropicBiLinear(E, σY, K)

# Define section
# =======================================
b = 0.1
h = 0.1
secName = "rectangle"
secParams = [b, h]

# Section struct
#StrSections = Section(secName, secParams)
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
nLoadSteps = 100 # Number of load increments
loadFactorsVec = ones(nLoadSteps) # Load scaling factors

# Numerical method settings struct
StrAnalysisSettings = ArcLength(tolk, tolu, tolf, loadFactorsVec)

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



# Check First step
# --------------------------------
P = abs(Fy)
Iy = StrSections.Iy
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

Mana = zeros(nLoadSteps)
C = E * K / (E + K)
epsY = σY / E
eps_ast = epsY - σY / C
kappa_ast = 2 * eps_ast / h
elem = 1
for i in 1:nLoadSteps
    kappak = kappaHistElem[elem, i]
    if kappak <= kappae
        Mana[i] = E * StrSections.Iy * kappak
    elseif kappak <= kappa_ast
        #else
        Mana[i] = σY * b * h^2 / 12 * (3 - kappae^2 / kappak^2 + kappak / kappae * C / E * (2 - 3 * kappae / kappak + kappae^3 / kappak^3))
    else
        zy = epsY / kappak
        z0 = eps_ast / kappak
        zs = z0 - zy
        Mana[i] = 2 * zy^2 * σY * b / 3 + zs * σY * b * (zy + (z0 - zy) / 3)
    end
end

My = σY * b * h^2 / 6

fig = plot(kappaHistElem[elem, :], Mana, markershape=:circle, lw=lw, ms=ms)
plot!(fig, kappaHistElem[elem, :], abs.(mVec), markershape=:rect, lw=lw, ms=ms)

err = (abs.(mVec[2:end]) - Mana[2:end]) ./ Mana[2:end] * 100
maxErr = maximum(err)

# Check KTe
Uke = zeros(4)
l = 1
Iy = StrSections.Iy
E = StrMaterialModels.E

Finte, KTe = finte_KT_int(StrMaterialModels, l, [b, h], Uke, 1)
Kana = rotXYXZ * E * Iy / l^3 * [12 6l -12 6l; 6l 4l^2 -6l 2l^2; -12 -6l 12 -6l; 6l 2l^2 -6l 4l^2] * rotXYXZ

mNum = "mNum.txt"
f = open(mNum, "w")

for i in mVec
    println(f, abs(i))
end


close(f)

kappa = "kappa.txt"
f = open(kappa, "w")

for i in kappaHistElem[elem, :]
    println(f, i)
end

close(f)
