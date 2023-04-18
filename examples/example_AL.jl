# ===============================================
# Cantilever with elasto plastic material model 
# ===============================================

# Load solver module
using materialNonLinearity, Plots, LinearAlgebra, Polynomials

# example name
problemName = "CantileverEPP"

# Define material model
# =======================================
E = 210e6
σY = 250e3
K = E / 10
ne = 12
ns = 12

# Materials struct
StrMaterialModels = IsotropicBiLinear(E, σY, K, ne, ns)

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
tolu = 1e-8 # Tolerance of converged disps
tolf = 1e-6 # Tolerance of internal forces
initialDeltaLambda = 1e-5 #

arcLengthIncrem = vcat(ones(45) * 1e-3)
arcLengthIncrem = vcat(ones(20) * 1e-3, ones(80) * 5e-4) # norma disps
# arcLengthIncrem = vcat(ones(2) * 4.36e-2, ones(15) * 3.615e-3) # norma todo igual a 1
nLoadSteps = length(arcLengthIncrem) # Number of load increments

controlDofs = [6] #

# Numerical method settings struct
StrAnalysisSettings = ArcLength_Cylindrical(tolk, tolu, tolf, nLoadSteps, initialDeltaLambda, arcLengthIncrem, controlDofs)

# arcLengthIncrem = vcat(ones(2) * 1.81e-2, ones(10) * 1e-3)
# nLoadSteps = length(arcLengthIncrem) # Number of load increments
# scalingProjection = 1
# StrAnalysisSettings = ArcLength(tolk, tolu, tolf, nLoadSteps, initialDeltaLambda, arcLengthIncrem, controlDofs, scalingProjection)

# Stress Array
# =======================================
elems = [1]
xG_Rel_Ind = collect(1:ne)

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
elem = 1
dofM = 3

# Loaded node
dofD = nnodes * 3 - 1
dofT = nnodes * 3

# Applied loads
pVec = sol.loadFactors * P

# Reaction Bending moment 
mVec = hcat([i[dofM] for i in matFint[elem]])

# Displacements at loaded node
dVec = hcat([i[dofD] for i in matUk])

# Compute curvatures
# --------------------------------
xrel = zeros(nelems)
kappaHistElem = frame_curvature(nelems, StrMesh, nLoadSteps, matUk, xrel)

# Analytical solution M-κ
# --------------------------------
Mana = zeros(nLoadSteps)
Iy = StrSections.Iy
κe = 2 * σY / (E * h)

C = E * K / (E + K)
εY = σY / E
ε⃰ = εY - σY / C
κ⃰ = 2 * ε⃰ / h
for i in 1:nLoadSteps
    κₖ = abs(kappaHistElem[elem, i])
    if κₖ <= κe
        Mana[i] = E * StrSections.Iy * κₖ
    else
        Mana[i] = σY * b * h^2 / 12 * (3 - κe^2 / κₖ^2 + κₖ / κe * C / E * (2 - 3 * κe / κₖ + κe^3 / κₖ^3))
    end
end

My = σY * b * h^2 / 6

# Plot parameters
# =======================================
include("../src/Utils/plots.jl")

lw = 3
ms = 2
color = "black"
minorGridBool = 1
legend_pos = :topright

StrPlots = PlotSettings(lw, ms, color, minorGridBool, legend_pos)

# Deformed shape plot
ndivs = 2
timesPlot = [1, nLoadSteps]

include("../src/Utils/plots.jl")

figsD = DeformedShapePlot(timesPlot, StrMesh, StrPlots, matUk)

# Constitutive model plot
SEfig = ConstitutiveModelPlot(StrMaterialModels, [-epsY * 3, epsY * 3], 50, 1000.0, 1e-3)

# M-κ plot
# --------------------------------
fig = plot(abs.(kappaHistElem[elem, :]), Mana, markershape=:circle, lw=lw, ms=ms, title="M-κ", label="Analytic", minorgrid=1, draw_arrow=1, legend=:bottomright)
plot!(fig, abs.(kappaHistElem[elem, :]), abs.(mVec), markershape=:rect, lw=lw, ms=ms, label="FEM")
xlabel!("κ")
ylabel!("M")

err = (abs.(mVec[2:end]) - Mana[2:end]) ./ Mana[2:end] * 100
maxErrMk = maximum(err)

# P-δ plot  
# --------------------------------
fig2 = plot(abs.(dVec), pVec, markershape=:circle, lw=lw, ms=ms, title="P-δ", label="FEM", minorgrid=1, draw_arrow=1, legend=:bottomright)
xlabel!("δ")
ylabel!("P")

# Stress plot  
# --------------------------------
p, w = gausslegendre(ns)

tf = length(arcLengthIncrem)

sfig = plot(σArr[1][tf][1], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[end]), minorgrid=1, draw_arrow=1, legend=:topright)
plot!(sfig, σArr[1][tf][10], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[end]), minorgrid=1, draw_arrow=1)
plot!(sfig, σArr[1][tf][end], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[end]), minorgrid=1, draw_arrow=1)
plot!(sfig, zeros(length(p)), p * h / 2, lw=lw, ms=ms, label="", color=:"black")

sfig1 = plot(σArr[1][tf-1][1], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tf-1]), minorgrid=1, draw_arrow=1, legend=:topright)
plot!(sfig1, σArr[1][tf-1][10], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tf-1]), minorgrid=1, draw_arrow=1)
plot!(sfig1, σArr[1][tf-1][end], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tf-1]), minorgrid=1, draw_arrow=1)
plot!(sfig1, zeros(length(p)), p * h / 2, lw=lw, ms=ms, label="", color=:"black")


fig