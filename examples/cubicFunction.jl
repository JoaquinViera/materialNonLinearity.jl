# ===============================================
# Cantilever with elasto plastic material model 
# ===============================================

# Load solver module
using materialNonLinearity, Plots, LinearAlgebra, FastGaussQuadrature

# example name
problemName = "cubicFunction"

# Define material model
# =======================================
E = 210e6
σY = 250e3
ne = 20
ns = 20

import materialNonLinearity: constitutive_model


# Materials struct
StrMaterialModels = UserModel(ne, ns)

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

# Numerical method parameters
# =======================================

tolk = 50 # number of iters
tolu = 1e-7 # Tolerance of converged disps
tolf = 1e-6 # Tolerance of internal forces
initialDeltaLambda = 1e-5 #
# arcLengthIncrem = vcat(ones(10) * 4e-4) #
#arcLengthIncrem = vcat(ones(15) * 1e-4, ones(5) * 1e-5, ones(15) * 2e-6, ones(619) * 1e-6, ones(2) * 5e-7) #
arcLengthIncrem = vcat(ones(15) * 1e-4, ones(5) * 1e-5, ones(16) * 2e-6, ones(619) * 2e-6) #
nLoadSteps = length(arcLengthIncrem) # Number of load increments
controlDofs = [6] #
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
elems = [1, nnodes - 1]
xG_Rel_Ind = 0

StrStressArray = StressArraySets(elems, xG_Rel_Ind)

# ===============================================
# Process model parameters
# ===============================================

sol, time, IterData, σArr = solver(StrSections, StrMaterialModels, StrMesh, StrBoundaryConds, StrAnalysisSettings, problemName, StrStressArray)

println(IterData.stopCrit)

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

# Loaded node
dofD = nnodes * 3 - 1
dofT = nnodes * 3

# Reaction Bending moment 
mVec = matFint[dofM, :]

# Reaction Bending moment 
mVec = matFint[dofM, :]
println(mVec[end])
dVec = hcat([i[dofD] for i in matUk])
tVec = hcat([i[dofT] for i in matUk])
pVec = abs.(mVec / L)

# Compute curvatures
# --------------------------------

kappaHistElem = zeros(nelems, nLoadSteps)

rotXYXZ = Diagonal(ones(4, 4))
rotXYXZ[2, 2] = -1
rotXYXZ[4, 4] = -1
dofsbe = [2, 3, 5, 6]

for j in 1:nelems
    nodeselem = StrMesh.conecMat[j, 3]
    local elemdofs = nodes2dofs(nodeselem[:], 3)
    local R, l = element_geometry(StrMesh.nodesMat[nodeselem[1], :], StrMesh.nodesMat[nodeselem[2], :], 3)
    Be = intern_function(0, l) * rotXYXZ
    for i in 1:nLoadSteps
        #elemdofs
        UkeL = R[dofsbe, dofsbe]' * matUk[i][elemdofs[dofsbe]]
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

p, w = gausslegendre(ns)
figspath = "..\\paper_matnonliniden\\tex\\2_Informe\\figs\\"

elem = 1
fig = plot(kappaHistElem[elem, :], abs.(mVec), markershape=:circle, lw=lw, ms=ms, title="M-κ", label="FEM", minorgrid=1, draw_arrow=1, legend=:bottomright)
#plot!(fig, kappaHistElem[elem, :], Mana, markershape=:rect, lw=lw, ms=ms, label="Analytic")
xlabel!("κ")
ylabel!("M")

savefig(fig, "$(figspath)ejemplo5M-k.png")

# P-δ plot  
# --------------------------------
fig2 = plot(abs.(dVec), pVec, markershape=:circle, lw=lw, ms=ms, title="P-δ", label="FEM", minorgrid=1, draw_arrow=1, legend=:bottomright)
xlabel!("δ")
ylabel!("P")

savefig(fig2, "$(figspath)ejemplo5P-d.png")

#plotlyjs()
fig3 = plot3d(abs.(dVec)[:, 1], tVec[:, 1], pVec, markershape=:circle, lw=lw, ms=ms, title="(δ,θ,P)", label="FEM", minorgrid=1, draw_arrow=1, legend=:bottomright)

sfig = plot(σArr[1][convert(Int, ceil(nLoadSteps / 10))], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[convert(Int, ceil(nLoadSteps / 5))]), minorgrid=1, draw_arrow=1, legend=:bottomright)
plot!(sfig, σArr[1][convert(Int, ceil(2 * nLoadSteps / 10))], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[convert(Int, ceil(2 * nLoadSteps / 10))]), minorgrid=1, draw_arrow=1)
plot!(sfig, σArr[1][convert(Int, ceil(3 * nLoadSteps / 10))], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[convert(Int, ceil(3 * nLoadSteps / 10))]), minorgrid=1, draw_arrow=1)
plot!(sfig, σArr[1][convert(Int, ceil(4 * nLoadSteps / 10))], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[convert(Int, ceil(4 * nLoadSteps / 10))]), minorgrid=1, draw_arrow=1)
plot!(sfig, σArr[1][convert(Int, ceil(5 * nLoadSteps / 10))], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[convert(Int, ceil(5 * nLoadSteps / 10))]), minorgrid=1, draw_arrow=1)
plot!(sfig, σArr[1][convert(Int, ceil(6 * nLoadSteps / 10))], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[convert(Int, ceil(6 * nLoadSteps / 10))]), minorgrid=1, draw_arrow=1)
plot!(sfig, σArr[1][convert(Int, ceil(7 * nLoadSteps / 10))], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[convert(Int, ceil(7 * nLoadSteps / 10))]), minorgrid=1, draw_arrow=1)
plot!(sfig, σArr[1][convert(Int, ceil(8 * nLoadSteps / 10))], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[convert(Int, ceil(8 * nLoadSteps / 10))]), minorgrid=1, draw_arrow=1)
plot!(sfig, σArr[1][convert(Int, ceil(9 * nLoadSteps / 10))], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[convert(Int, ceil(9 * nLoadSteps / 10))]), minorgrid=1, draw_arrow=1)
plot!(sfig, σArr[1][end], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[end]), minorgrid=1, draw_arrow=1)
plot!(sfig, zeros(length(p)), p * h / 2, lw=lw, ms=ms, label="", color=:"black")

savefig(sfig, "$(figspath)ejemplo5stress1.png")