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
K = -E / 100
ne = 16
ns = 16

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
nnodes = 101
xcoords = collect(LinRange(0, L, nnodes))
#t = 0:1/(nnodes-1):1
#xcoords = t .^ 1.5 * L
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

initialDeltaLambda = 1e-2 #

arcLengthIncrem = vcat(ones(18) * 1e-3, ones(80) * 1e-4) # 21 nodes
arcLengthIncrem = vcat(ones(28) * 1e-3, ones(23) * 1e-4, ones(200) * 1e-5) # 51 nodes
arcLengthIncrem = vcat(ones(39) * 1e-3, ones(21) * 1e-4, ones(30) * 1e-5) # 101 nodes
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
lw = 3
ms = 2
color = "black"
minorGridBool = 1
legend_pos = :topright

StrPlots = PlotSettings(lw, ms, color, minorGridBool, legend_pos)

figspath = "..\\paper_matnonliniden\\tex\\2_Informe\\figs\\"

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

#=
# convergence analysis
κₚ = 0.4 # 1/m
δκ = 4e-2
kappaVecElem = kappaHistElem[elem, :]
idx = maximum(findall(x -> x <= δκ, (abs.(kappaVecElem .- κₚ))))

kappa_i = kappaVecElem[idx]
M_i = abs(mVec[idx])

if kappa_i <= κₑ
    M_a = E * StrSections.Iy * kappa_i
elseif kappa_i <= kappa_ast
    M_a = σY * b * h^2 / 12 * (3 - κₑ^2 / kappa_i^2 + kappa_i / κₑ * C / E * (2 - 3 * κₑ / kappa_i + κₑ^3 / kappa_i^3))
else
    zy = epsY / kappa_i
    z0 = eps_ast / kappa_i
    zs = z0 - zy
    M_a = 2 * zy^2 * σY * b / 3 + zs * σY * b * (zy + (z0 - zy) / 3)
end

err = (M_i - M_a) / M_a * 100

println(kappa_i)
println(M_i)
println(M_a)
println(err)
=#


# mNum & κ

#=
mNum = "mNum10.txt"
f = open(mNum, "w")

for i in mVec
    println(f, abs(i))
end

close(f)

κ = "kappa10.txt"
f = open(κ, "w")

for i in kappaHistElem[elem, :]
    println(f, i)
end

close(f)

errVec = [5.65, 4.8, 4.26, 4.13]
ne = [20, 40, 80, 100]

#xs = range(0, 150, step=1)
#p = fit(ne, errVec, 3)

figErr = plot(ne, errVec, markershape=:circle, lw=lw, ms=ms, title="err-nₑ", label="err-nₑ", minorgrid=1, draw_arrow=1)
#plot!(figErr, xs, p.(xs), markershape=:circle, lw=lw, ms=ms, label="Fit")
xlabel!("nₑ")
ylabel!("err")

savefig(figErr, "ejemplo3err.png")


figConv = plot(kappa20, M20, markershape=:circle, lw=0.5 * lw, ms=ms, title="M-κ", label="nₑ=20", minorgrid=1, draw_arrow=1, legend=:bottomright)
plot!(kappa40[1:end-18], M40[1:end-18], markershape=:circle, lw=0.5 * lw, ms=ms, title="M-κ", label="nₑ=40", minorgrid=1, draw_arrow=1)
plot!(kappa80[1:end-10], M80[1:end-10], markershape=:circle, lw=0.5 * lw, ms=ms, title="M-κ", label="nₑ=80", minorgrid=1, draw_arrow=1)
plot!(kappa100[1:end-6], M100[1:end-6], markershape=:circle, lw=0.5 * lw, ms=ms, title="M-κ", label="nₑ=100", minorgrid=1, draw_arrow=1)
xlabel!("κ")
ylabel!("M")
savefig(figConv, "ejemplo3conv.png")
=#