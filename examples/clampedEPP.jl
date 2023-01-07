# ===============================================
# Cantilever with elasto plastic material model 
# ===============================================

# Load solver module
using materialNonLinearity, Plots, LinearAlgebra, Polynomials

# example name
problemName = "ClampedEPP"

# Define material model
# =======================================
E = 210e6
σY = 250e3
K = 0

# Gauss points
ne = 16 # element
ns = 16 # section

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
supps = [1 Inf Inf Inf; nnodes Inf Inf Inf]

# Define applied external loads
Fx = 0
Fz = -1
My = 0
nod = (nnodes + 1) / 2
nodalForces = [nod Fx Fz My]

# BoundaryConds struct
StrBoundaryConds = BoundaryConds(supps, nodalForces)

# Numerical method parameters
# =======================================

tolk = 50 # number of iters
tolu = 1e-7 # Tolerance of converged disps
tolf = 1e-6 # Tolerance of internal forces

initialDeltaLambda = 1e-2 #

arcLengthIncrem = vcat(ones(18) * 1e-3, ones(80) * 1e-4) # 21 node
arcLengthIncrem = vcat(ones(18) * 1e-3) # 21 node
nLoadSteps = length(arcLengthIncrem) # Number of load increments
controlDofs = [6] #
scalingProjection = 1 #

# Numerical method settings struct
StrAnalysisSettings = ArcLength(tolk, tolu, tolf, nLoadSteps, initialDeltaLambda, arcLengthIncrem, controlDofs, scalingProjection)

# Stress Array
# =======================================
elems = []
xG_Rel_Ind = 0

StrStressArray = StressArraySets(elems, xG_Rel_Ind)

# ===============================================
# Process model parameters
# ===============================================

sol, time, IterData = solver(StrSections, StrMaterialModels, StrMesh, StrBoundaryConds, StrAnalysisSettings, problemName, StrStressArray)

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
kappaHistElem = frame_curvature(nelems, StrMesh, nLoadSteps, matUk)

# Analytical solution M-κ
# --------------------------------
Mana = zeros(nLoadSteps)
κₑ = 2 * σY / (E * h)
C = E * K / (E + K)
epsY = σY / E
eps_ast = epsY - σY / C
kappa_ast = 2 * eps_ast / h
elem = 1
for i in 1:nLoadSteps
    κₖ = abs(kappaHistElem[elem, i])
    if κₖ <= κₑ
        Mana[i] = E * StrSections.Iy * κₖ
    elseif κₖ <= kappa_ast
        Mana[i] = σY * b * h^2 / 12 * (3 - κₑ^2 / κₖ^2 + κₖ / κₑ * C / E * (2 - 3 * κₑ / κₖ + κₑ^3 / κₖ^3))
    else
        zy = epsY / κₖ
        z0 = eps_ast / κₖ
        zs = z0 - zy
        Mana[i] = 2 * zy^2 * σY * b / 3 + zs * σY * b * (zy + (z0 - zy) / 3)
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

figspath = "..\\paper_matnonliniden\\tex\\2_Informe\\figs\\"

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

savefig(fig, "$(figspath)ejemplo3M-k.png")

err = (abs.(mVec[2:end]) - Mana[2:end]) ./ Mana[2:end] * 100
maxErrMk = maximum(err)

# P-δ plot  
# --------------------------------
fig2 = plot(abs.(dVec), pVec, markershape=:circle, lw=lw, ms=ms, title="P-δ", label="FEM", minorgrid=1, draw_arrow=1, legend=:bottomright)
xlabel!("δ")
ylabel!("P")

savefig(fig2, "$(figspath)ejemplo3P-d.png")

# Bending moment plot
# --------------------------------
ndivs = 2
timesPlot = [1, nLoadSteps]

include("../src/Utils/plots.jl")

figsM = BendingMomentPlot(timesPlot, StrMesh, StrPlots, matFint)
