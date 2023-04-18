# ===============================================
# Cantilever with elasto plastic material model 
# ===============================================

# Load solver module
using materialNonLinearity, Plots, LinearAlgebra, FastGaussQuadrature, Printf

# example name
problemName = "FRConcrete"

# Define material model
# =======================================
# Tension
# ---------------
# Tramo 1
fctd = 1.97 * 1000 # kN/m2
epsf = 0.06 / 1000
# Tramo 2
fctR1d = 0.75 * 1000 # kN/m2
eps1 = 0.16 / 1000
# Tramo 3
fctR3d = 0.52 * 1000 # kN/m2
eps2 = 12.5 / 1000
# Tramo 4
fctlim = 0.38 * 1000 # kN/m2
epslim = 20 / 1000

E = fctd / epsf # kN/m2

# Gauss points
ne = 20
ns = 50

import materialNonLinearity: constitutive_model

# Materials struct
StrMaterialModels = UserModel(ne, ns)

function constitutive_model(ElemMaterialModel::UserModel, εₖ)

    # Tension
    # ---------------
    # Tramo 1
    fctd = 1.97 * 1000 # kN/m2
    epsf = 0.06 / 1000
    # Tramo 2
    fctR1d = 0.75 * 1000 # kN/m2
    eps1 = 0.16 / 1000
    # Tramo 3
    fctR3d = 0.52 * 1000 # kN/m2
    eps2 = 12.5 / 1000
    # Tramo 4
    fctlim = 0.38 * 1000 # kN/m2
    epslim = 20 / 1000

    E = fctd / epsf # kN/m2

    # Tension
    if εₖ >= 0.0
        if εₖ <= epsf # Tramo 1
            σ = E * εₖ
            ∂σ∂ε = E
        elseif εₖ <= eps1 # Tramo 2
            ∂σ∂ε = (fctR1d - fctd) / (eps1 - epsf)
            σ = ∂σ∂ε * (εₖ - epsf) + fctd
        elseif εₖ <= eps2 # Tramo 3
            ∂σ∂ε = (fctR3d - fctR1d) / (eps2 - eps1)
            σ = ∂σ∂ε * (εₖ - eps1) + fctR1d
        elseif εₖ <= epslim # Tramo 4
            ∂σ∂ε = (fctlim - fctR3d) / (epslim - eps2)
            σ = ∂σ∂ε * (εₖ - eps2) + fctR3d
        elseif εₖ > epslim # Rotura
            error("Collapse")
        end
    else #Compression
        σ = E * εₖ
        ∂σ∂ε = E
    end

    return σ, ∂σ∂ε

end

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

tolk = 75 # number of iters
tolu = 1e-10 # Tolerance of converged disps
tolf = 1e-6 # Tolerance of internal forces
initialDeltaLambda = 1e-7 #
arcLengthIncrem = vcat(ones(30) * 1e-5, ones(100) * 3e-6)
nLoadSteps = length(arcLengthIncrem)
controlDofs = [6] #
scalingProjection = 1 #

# Numerical method settings struct
StrAnalysisSettings = ArcLength(tolk, tolu, tolf, nLoadSteps, initialDeltaLambda, arcLengthIncrem, controlDofs, scalingProjection)

# Stress Array
# =======================================
elems = [1, 2, 3, 4, 5]
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
σY = fctd
Mfis = σY * Iy / (h / 2)
println(Mfis)

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

# Plot parameters
# =======================================
include("../../src/Utils/plots.jl")
lw = 3
ms = 2
color = "black"
minorGridBool = 1
legend_pos = :topright

StrPlots = PlotSettings(lw, ms, color, minorGridBool, legend_pos)

figspath = "..\\..\\paper_matnonliniden\\tex\\2_Informe\\figs\\"

# Constitutive model plot

SEfig = ConstitutiveModelPlot(StrMaterialModels, [-epslim / 300, epslim], 1000, 1000.0, 1e-3)

# savefig(SEfig, "$(figspath)ejemplo6sigma-epsilon.png")


# stop
# M-κ plot  
# --------------------------------
elem = 1
fig = plot(abs.(kappaHistElem[elem, :]), abs.(mVec), markershape=:circle, lw=lw, ms=ms, title="M-κ", label="FEM", minorgrid=1, draw_arrow=1, legend=:bottomright)
xlabel!("κ")
ylabel!("M")

# savefig(fig, "$(figspath)ejemplo6M-k.png")

# P-δ plot  
# --------------------------------
fig2 = plot(abs.(dVec), pVec, markershape=:circle, lw=lw, ms=ms, title="P-δ", label="FEM", minorgrid=1, draw_arrow=1, legend=:bottomright)
xlabel!("δ")
ylabel!("P")

# savefig(fig2, "$(figspath)ejemplo6P-d.png")

# Stress plot  
# --------------------------------
p, w = gausslegendre(ns)

sfig = plot(σArr[1][convert(Int, ceil(nLoadSteps / 5))], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[convert(Int, ceil(nLoadSteps / 5))]), minorgrid=1, draw_arrow=1, legend=:topleft)
plot!(sfig, σArr[1][convert(Int, ceil(2 * nLoadSteps / 5))], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[convert(Int, ceil(2 * nLoadSteps / 5))]), minorgrid=1, draw_arrow=1)
plot!(sfig, σArr[1][convert(Int, ceil(3 * nLoadSteps / 5))], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[convert(Int, ceil(3 * nLoadSteps / 5))]), minorgrid=1, draw_arrow=1)
plot!(sfig, σArr[1][convert(Int, ceil(4 * nLoadSteps / 5))], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[convert(Int, ceil(4 * nLoadSteps / 5))]), minorgrid=1, draw_arrow=1)
plot!(sfig, σArr[1][end], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[end]), minorgrid=1, draw_arrow=1)
plot!(sfig, zeros(length(p)), p * h / 2, lw=lw, ms=ms, label="", color=:"black")

# savefig(sfig, "$(figspath)ejemplo6stress1.png")

tp = 42
tf = nLoadSteps
# prueba
sfig2 = plot(σArr[2][tp][1], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tp]), minorgrid=1, draw_arrow=1, legend=:bottomright)
plot!(sfig2, σArr[2][tp-1][10], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tp-1]), minorgrid=1, draw_arrow=1)
plot!(sfig2, σArr[2][tp-2][10], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tp-2]), minorgrid=1, draw_arrow=1)
plot!(sfig2, σArr[2][tp-4][10], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tp-4]), minorgrid=1, draw_arrow=1)
plot!(sfig2, σArr[2][tp-6][10], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tp-6]), minorgrid=1, draw_arrow=1)

sfig22 = plot(σArr[2][tp][1], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tp]), minorgrid=1, draw_arrow=1, legend=:bottomright)
plot!(sfig22, σArr[2][tp+1][10], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tp+1]), minorgrid=1, draw_arrow=1)
plot!(sfig22, σArr[2][tp+2][10], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tp+2]), minorgrid=1, draw_arrow=1)
plot!(sfig22, σArr[2][tp+3][10], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tp+3]), minorgrid=1, draw_arrow=1)
plot!(sfig22, σArr[2][tp+4][10], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tp+4]), minorgrid=1, draw_arrow=1)
plot!(sfig22, σArr[2][tp+5][10], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tp+5]), minorgrid=1, draw_arrow=1)
plot!(sfig22, σArr[2][tp+6][10], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tp+6]), minorgrid=1, draw_arrow=1)
# plot!(sfig22, σArr[2][tp+7][10], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tp+7]), minorgrid=1, draw_arrow=1)
plot!(sfig22, zeros(length(p)), p * h / 2, lw=lw, ms=ms, label="", color=:"black")

sfig33 = plot(σArr[2][tp][1], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tp]), minorgrid=1, draw_arrow=1, legend=:bottomright)
plot!(sfig33, σArr[2][tp+1][10], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tp+1]), minorgrid=1, draw_arrow=1)
plot!(sfig33, σArr[2][tp+2][10], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tp+2]), minorgrid=1, draw_arrow=1)
plot!(sfig33, σArr[2][tp+3][10], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tp+3]), minorgrid=1, draw_arrow=1)
plot!(sfig33, σArr[2][tp+4][10], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tp+4]), minorgrid=1, draw_arrow=1)
plot!(sfig33, σArr[2][tp+5][10], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tp+5]), minorgrid=1, draw_arrow=1)
plot!(sfig33, σArr[2][tp+6][10], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tp+6]), minorgrid=1, draw_arrow=1)
# plot!(sfig22, σArr[2][tp+7][10], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tp+7]), minorgrid=1, draw_arrow=1)
plot!(sfig33, zeros(length(p)), p * h / 2, lw=lw, ms=ms, label="", color=:"black")

sfig44 = plot(σArr[3][tp][1], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tp]), minorgrid=1, draw_arrow=1, legend=:bottomright)
plot!(sfig44, σArr[3][tp+1][10], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tp+1]), minorgrid=1, draw_arrow=1)
plot!(sfig44, σArr[3][tp+2][10], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tp+2]), minorgrid=1, draw_arrow=1)
plot!(sfig44, σArr[3][tp+3][10], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tp+3]), minorgrid=1, draw_arrow=1)
plot!(sfig44, σArr[3][tp+4][10], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tp+4]), minorgrid=1, draw_arrow=1)
plot!(sfig44, σArr[3][tp+5][10], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tp+5]), minorgrid=1, draw_arrow=1)
plot!(sfig44, σArr[3][tp+6][10], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tp+6]), minorgrid=1, draw_arrow=1)
# plot!(sfig22, σArr[2][tp+7][10], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tp+7]), minorgrid=1, draw_arrow=1)
plot!(sfig44, zeros(length(p)), p * h / 2, lw=lw, ms=ms, label="", color=:"black")

sfig5 = plot(σArr[4][tf][1], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tp]), minorgrid=1, draw_arrow=1, legend=:bottomright)
plot!(sfig5, σArr[4][tf-1][10], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tp+1]), minorgrid=1, draw_arrow=1)
plot!(sfig5, σArr[4][tf-2][10], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tp+2]), minorgrid=1, draw_arrow=1)
plot!(sfig5, σArr[4][tf-3][10], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tp+3]), minorgrid=1, draw_arrow=1)
plot!(sfig5, σArr[4][tf-4][10], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tp+4]), minorgrid=1, draw_arrow=1)
plot!(sfig5, σArr[4][tf-5][10], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tp+5]), minorgrid=1, draw_arrow=1)
plot!(sfig5, σArr[4][tf-6][10], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tp+6]), minorgrid=1, draw_arrow=1)
# plot!(sfig22, σArr[2][tp+7][10], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tp+7]), minorgrid=1, draw_arrow=1)
plot!(sfig5, zeros(length(p)), p * h / 2, lw=lw, ms=ms, label="", color=:"black")

sfig55 = plot(σArr[4][tp][1], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tp]), minorgrid=1, draw_arrow=1, legend=:bottomright)
plot!(sfig55, σArr[4][tp+1][10], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tp+1]), minorgrid=1, draw_arrow=1)
plot!(sfig55, σArr[4][tp+2][10], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tp+2]), minorgrid=1, draw_arrow=1)
plot!(sfig55, σArr[4][tp+3][10], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tp+3]), minorgrid=1, draw_arrow=1)
plot!(sfig55, σArr[4][tp+4][10], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tp+4]), minorgrid=1, draw_arrow=1)
plot!(sfig55, σArr[4][tp+5][10], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tp+5]), minorgrid=1, draw_arrow=1)
plot!(sfig55, σArr[4][tp+6][10], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tp+6]), minorgrid=1, draw_arrow=1)
# plot!(sfig22, σArr[2][tp+7][10], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[tp+7]), minorgrid=1, draw_arrow=1)
plot!(sfig55, zeros(length(p)), p * h / 2, lw=lw, ms=ms, label="", color=:"black")


#
sfig2 = plot(σArr[end][convert(Int, ceil(nLoadSteps / 5))], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[convert(Int, ceil(nLoadSteps / 5))]), minorgrid=1, draw_arrow=1, legend=:topleft)
plot!(sfig2, σArr[end][convert(Int, ceil(2 * nLoadSteps / 5))], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[convert(Int, ceil(2 * nLoadSteps / 5))]), minorgrid=1, draw_arrow=1)
plot!(sfig2, σArr[end][convert(Int, ceil(3 * nLoadSteps / 5))], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[convert(Int, ceil(3 * nLoadSteps / 5))]), minorgrid=1, draw_arrow=1)
plot!(sfig2, σArr[end][convert(Int, ceil(4 * nLoadSteps / 5))], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[convert(Int, ceil(4 * nLoadSteps / 5))]), minorgrid=1, draw_arrow=1)
plot!(sfig2, σArr[end][end], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[end]), minorgrid=1, draw_arrow=1)
plot!(sfig2, zeros(length(p)), p * h / 2, lw=lw, ms=ms, label="", color=:"black")

# savefig(sfig2, "$(figspath)ejemplo6stress2.png")

# Bending moment plot
# --------------------------------
ndivs = 2
timesPlot = [1, 130]

figsM = BendingMomentPlot(timesPlot, StrMesh, StrPlots, matFint)

# savefig(figsM[end], "$(figspath)ejemplo6bending.png")

# Deformed shape plot
# --------------------------------
ndivs = 2
timesPlot = [1, 15, 30, 45, 60, 75, 90, 105, 130]

figsD = DeformedShapePlot(timesPlot, StrMesh, StrPlots, matUk)

figsDefs = plot(figsD[end], figsD[end-1], title="Deformed shapes")

# savefig(figsDefs, "$(figspath)ejemplo6deformed.png")
