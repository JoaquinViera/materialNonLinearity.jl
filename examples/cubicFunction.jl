# ===============================================
# Cantilever with elasto plastic material model 
# ===============================================

# Load solver module
using materialNonLinearity, Plots, LinearAlgebra

# example name
problemName = "cubicFunction"

# Define material model
# =======================================
E = 210e6
σY = 250e3

import materialNonLinearity: constitutive_model


# Materials struct
StrMaterialModels = UserModel()

function constitutive_model(ElemMaterialModel::UserModel, εₖ)

    E = 210e6
    sigmaY = 250e3
    epsY = sigmaY / E

    a = -sigmaY / (2 * epsY^3)
    b = 3 * sigmaY / (2 * epsY)


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
initialDeltaLambda = 1e-3 #
#arcLengthIncrem = vcat(ones(30) * 1e-4, ones(5) * 5e-5, ones(15) * 1e-5) #
arcLengthIncrem = vcat(ones(40) * 4e-4, ones(1) * 1e-5) #
nLoadSteps = length(arcLengthIncrem) # Number of load increments
controlDofs = [10, 16, 24] #
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

#tVec = abs.(matUk[dofT, :])
#dVec = abs.(matUk[dofD, :])
pVec = abs.(mVec / L)

# Compute curvatures
# --------------------------------

kappaHistElem = zeros(nelems, nLoadSteps)

rotXYXZ = Diagonal(ones(4, 4))
rotXYXZ[2, 2] = -1
rotXYXZ[4, 4] = -1

for j in 1:nelems
    nodeselem = StrMesh.conecMat[j, 3]
    local elemdofs = nodes2dofs(nodeselem[:], 2)
    local R, l = element_geometry(StrMesh.nodesMat[nodeselem[1], :], StrMesh.nodesMat[nodeselem[2], :], 2)
    Be = intern_function(0, l) * rotXYXZ
    for i in 1:nLoadSteps
        #elemdofs
        UkeL = R' * matUk[i][elemdofs]
        kappaelem = Be * UkeL

        kappaHistElem[j, i] = abs.(kappaelem[1])
    end
end

# Analytical solution M-κ
# --------------------------------

Mana = zeros(nLoadSteps)
epsY = σY / E
ca = -σY / (2 * epsY^3);
cb = 3 * σY / (2 * epsY);

elem = 1
for i in 1:nLoadSteps
    kappak = kappaHistElem[elem, i]
    Mana[i] = kappak * b * (ca * kappak^2 * h^5 / 80 + cb * h^3 / 12)
end

err = (abs.(mVec[2:end]) - Mana[2:end]) ./ Mana[2:end] * 100
maxErrMk = maximum(err)
println(maxErrMk)

# M-κ plot  
# --------------------------------
elem = 1
fig = plot(kappaHistElem[elem, :], abs.(mVec), markershape=:circle, lw=lw, ms=ms, title="M-κ", label="FEM", minorgrid=1, draw_arrow=1, legend=:bottomright)
plot!(fig, kappaHistElem[elem, :], Mana, markershape=:rect, lw=lw, ms=ms, label="Analytic")
xlabel!("κ")
ylabel!("M")

# P-δ plot  
# --------------------------------
figPdelta = plot(dVec, pVec, markershape=:circle, lw=lw, ms=ms, title="P-δ", label="FEM", minorgrid=1, draw_arrow=1, legend=:bottomright)
xlabel!("δ")
ylabel!("P")

figPdelta2 = plot(dVec[end-15:end], pVec[end-15:end], markershape=:circle, lw=lw, ms=ms, title="P-δ", label="FEM", minorgrid=1, draw_arrow=1, legend=:bottomright)
xlabel!("δ")
ylabel!("P")

#plotlyjs()
fig3 = plot3d(dVec, tVec, pVec, markershape=:circle, lw=lw, ms=ms, title="(δ,θ,P)", label="FEM", minorgrid=1, draw_arrow=1, legend=:bottomright)
