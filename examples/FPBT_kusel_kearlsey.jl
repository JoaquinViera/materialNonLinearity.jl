# ===============================================
# Cantilever with elasto plastic material model 
# ===============================================

# Load solver module
using materialNonLinearity, Plots, LinearAlgebra, FastGaussQuadrature, Printf

# example name
problemName = "FRC_kusel_kearlsey"

# Define material model
# =======================================
# Tension
# ---------------
# Tramo 1
# fctd = 6.44 * 1000 # kN/m2
# epsf = 1.74e-4
# # Tramo 2
# ft1 = 2.8 * 1000 # kN/m2
# eps1 = 1.8e-4
# # Tramo 3
# ft2 = 3.4 * 1000 # kN/m2
# eps2 = 9.00e-3
# # Tramo 4
# ft3 = 1.35 * 1000 # kN/m2
# eps3 = 2.50e-2
# # Tramo 5
# fctlim = 0 # kN/m2
# epslim = 1.00e-1

# E = fctd / epsf # kN/m2
E = 6.44 * 1000 / 1.74e-4 # kN/m2

# Gauss points
ne = 20
ns = 50

import materialNonLinearity: constitutive_model

# Materials struct
StrMaterialModels = UserModel(ne, ns)

function constitutive_model(ElemMaterialModel::UserModel, εₖ)

    # E = fctd / epsf # kN/m2
    E = 6.44 * 1000 / 1.74e-4 # kN/m2

    # Tension - simplificada
    # ---------------
    # Tramo 1
    fctd = 2.8 * 1000 # kN/m2
    epsf = fctd / E
    # epsf = 1.74e-4
    # Tramo 2
    ft1 = 3.4 * 1000 # kN/m2
    eps1 = 9.00e-3
    # Tramo 3
    ft2 = 1.35 * 1000 # kN/m2
    eps2 = 2.50e-2
    # Tramo 4
    fctlim = 0 # kN/m2
    epslim = 1.00e-1

    # Tension - paper
    # ---------------
    # Tramo 1
    fctd = 6.44 * 1000 # kN/m2
    epsf = 1.74e-4
    # Tramo 2
    ft1 = 2.8 * 1000 # kN/m2
    eps1 = 1.74e-4 * 1.05
    # Tramo 3
    ft2 = 3.4 * 1000 # kN/m2
    eps2 = 9.00e-3
    # Tramo 4
    ft3 = 1.35 * 1000 # kN/m2
    eps3 = 2.50e-2
    # Tramo 5
    fctlim = 0 # kN/m2
    epslim = 1.00e-1


    # Tension
    if εₖ >= 0.0
        # if εₖ <= epsf # Tramo 1
        #     σ = E * εₖ
        #     ∂σ∂ε = E
        # elseif εₖ <= eps1 # Tramo 2
        #     ∂σ∂ε = (ft1 - fctd) / (eps1 - epsf)
        #     σ = ∂σ∂ε * (εₖ - epsf) + fctd
        # elseif εₖ <= eps2 # Tramo 3
        #     ∂σ∂ε = (ft2 - ft1) / (eps2 - eps1)
        #     σ = ∂σ∂ε * (εₖ - eps1) + ft1
        # elseif εₖ <= epslim # Tramo 4
        #     ∂σ∂ε = (fctlim - ft2) / (epslim - eps2)
        #     σ = ∂σ∂ε * (εₖ - eps2) + ft2
        # elseif εₖ > epslim # Rotura
        #     error("Collapse")
        # end
        if εₖ <= epsf # Tramo 1
            σ = E * εₖ
            ∂σ∂ε = E
        elseif εₖ <= eps1 # Tramo 2
            ∂σ∂ε = (ft1 - fctd) / (eps1 - epsf)
            σ = ∂σ∂ε * (εₖ - epsf) + fctd
        elseif εₖ <= eps2 # Tramo 3
            ∂σ∂ε = (ft2 - ft1) / (eps2 - eps1)
            σ = ∂σ∂ε * (εₖ - eps1) + ft1
        elseif εₖ <= eps3 # Tramo 4
            ∂σ∂ε = (ft3 - ft2) / (eps3 - eps2)
            σ = ∂σ∂ε * (εₖ - eps2) + ft2
        elseif εₖ <= epslim # Tramo 5
            ∂σ∂ε = (fctlim - ft3) / (epslim - eps3)
            σ = ∂σ∂ε * (εₖ - eps3) + ft3
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
b = 0.1
h = 0.1

# Section struct
StrSections = Rectangle(; b, h)

# Define Mesh
# =======================================

# Nodes
# L = 4.5
L = 0.45
# nnodes = 55
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
# mid = findall(x -> x == L / 2, xcoords)
supps = [1 Inf Inf 0; nnodes 0 Inf 0]

# Define applied external loads
# load_coord1 = 2.25 - 1
load_coord1 = L / 4
# load_coord2 = 2.25 + 1
load_coord2 = 3L / 4

delta_x = L / (nnodes - 1)
n1 = convert(Int64, round(load_coord1 / delta_x))
n2 = convert(Int64, round(load_coord2 / delta_x))

x1 = xcoords[n1+1]
x2 = xcoords[n2+1]

Fx = 0
Fz = -1
My = 0
nodalForces = [(n1+1) Fx Fz My; (n2+1) Fx Fz My]


# BoundaryConds struct
StrBoundaryConds = BoundaryConds(supps, nodalForces)

# Numerical method parameters
# =======================================

tolk = 75 # number of iters
tolu = 1e-10 # Tolerance of converged disps
tolf = 1e-6 # Tolerance of internal forces
initialDeltaLambda = 1e-5 #

arcLengthIncrem = vcat(ones(25) * 1e-5, ones(125) * 3e-6) # dips

nLoadSteps = length(arcLengthIncrem)
dof1 = (n1 + 1) * 3 - 1
dof2 = (n2 + 1) * 3 - 1
controlDofs = [dof1, dof2] #
# controlDofs = [dof1] #
scalingProjection = -1 #

# Numerical method settings struct
StrAnalysisSettings = ArcLength(tolk, tolu, tolf, nLoadSteps, initialDeltaLambda, arcLengthIncrem, controlDofs, scalingProjection)

# Stress Array
# =======================================
elems = [n1, n1 + 1]
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
σY = fctd
Mfis = σY * Iy / (h / 2)
println(Mfis)

# Numerical solution
matFint = sol.matFint
matUk = sol.matUk

# Clamped node
nod = 1
elem = n1 + 1
dofM = 3

# Loaded node
dofD = (n1 + 1) * 3 - 1
dofT = (n1 + 1) * 3

# Applied loads
pVec = sol.loadFactors * P

# Reaction Bending moment 
mVec = hcat([i[dofM] for i in matFint[elem]])
# vVec = hcat([i[2] for i in matFint[elem]])
# Displacements at loaded node
dVec = hcat([i[dofD] for i in matUk])
# tVec = hcat([i[dofT] for i in matUk])
# Compute curvatures
# --------------------------------
kappaHistElem = frame_curvature(nelems, StrMesh, nLoadSteps, matUk)

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

# Constitutive model plot

SEfig = ConstitutiveModelPlot(StrMaterialModels, [-epslim / 300, epslim], 1000, 1000.0, 1e-3)

savefig(SEfig, "$(figspath)ejemplo8sigma-epsilon.png")


# stop
# M-κ plot  
# --------------------------------
elem = n1 + 1
fig = plot(abs.(kappaHistElem[elem, :]), abs.(mVec), markershape=:circle, lw=lw, ms=ms, title="M-κ", label="FEM", minorgrid=1, draw_arrow=1, legend=:bottomright)
xlabel!("κ")
ylabel!("M")

# savefig(fig, "$(figspath)ejemplo6M-k.png")

# P-δ plot  
# --------------------------------
fig2 = plot(abs.(dVec) * 1000, pVec, markershape=:circle, lw=lw, ms=ms, title="P-δ", label="FEM", minorgrid=1, draw_arrow=1, legend=:bottomright)
xlabel!("δ")
ylabel!("P")

# savefig(fig2, "$(figspath)ejemplo6P-d.png")

# Stress plot  
# --------------------------------
p, w = gausslegendre(ns)

sfig = plot(σArr[1][convert(Int, ceil(nLoadSteps / 5))], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[convert(Int, ceil(nLoadSteps / 5))]), minorgrid=1, draw_arrow=1, legend=:topright)
plot!(sfig, σArr[1][convert(Int, ceil(2 * nLoadSteps / 5))], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[convert(Int, ceil(2 * nLoadSteps / 5))]), minorgrid=1, draw_arrow=1)
plot!(sfig, σArr[1][convert(Int, ceil(3 * nLoadSteps / 5))], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[convert(Int, ceil(3 * nLoadSteps / 5))]), minorgrid=1, draw_arrow=1)
plot!(sfig, σArr[1][convert(Int, ceil(4 * nLoadSteps / 5))], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[convert(Int, ceil(4 * nLoadSteps / 5))]), minorgrid=1, draw_arrow=1)
plot!(sfig, σArr[1][end], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[end]), minorgrid=1, draw_arrow=1)
plot!(sfig, zeros(length(p)), p * h / 2, lw=lw, ms=ms, label="", color=:"black")

# savefig(sfig, "$(figspath)ejemplo6stress1.png")

sfig2 = plot(σArr[end][convert(Int, ceil(nLoadSteps / 5))], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[convert(Int, ceil(nLoadSteps / 5))]), minorgrid=1, draw_arrow=1, legend=:topright)
plot!(sfig2, σArr[end][convert(Int, ceil(2 * nLoadSteps / 5))], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[convert(Int, ceil(2 * nLoadSteps / 5))]), minorgrid=1, draw_arrow=1)
plot!(sfig2, σArr[end][convert(Int, ceil(3 * nLoadSteps / 5))], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[convert(Int, ceil(3 * nLoadSteps / 5))]), minorgrid=1, draw_arrow=1)
plot!(sfig2, σArr[end][convert(Int, ceil(4 * nLoadSteps / 5))], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[convert(Int, ceil(4 * nLoadSteps / 5))]), minorgrid=1, draw_arrow=1)
plot!(sfig2, σArr[end][end], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label=@sprintf("M = %0.2f", mVec[end]), minorgrid=1, draw_arrow=1)
plot!(sfig2, zeros(length(p)), p * h / 2, lw=lw, ms=ms, label="", color=:"black")

savefig(sfig2, "$(figspath)ejemplo8stress2.png")

# # Bending moment plot
# # --------------------------------
ndivs = 2
timesPlot = [1, length(pVec)]

figsM = BendingMomentPlot(timesPlot, StrMesh, StrPlots, matFint)

savefig(figsM[end], "$(figspath)ejemplo8bending.png")


# Deformed shape plot
# --------------------------------
ndivs = 2
timesPlot = [1, length(pVec)]

figsD = DeformedShapePlot(timesPlot, StrMesh, StrPlots, matUk)

figsDefs = plot(figsD[end], figsD[end-1], title="Deformed shapes")

savefig(figsDefs, "$(figspath)ejemplo8deformed.png")

# Paper results

dVec_paper =
    [0.000
        0.013
        0.026
        0.033
        0.046
        0.047
        0.053
        0.072
        0.073
        0.143
        0.175
        0.226
        0.277
        0.359
        0.417
        0.461
        0.519
        0.557
        0.582
        0.659
        0.793
        0.875
        0.932
        0.977
        1.009
        1.041
        1.066
        1.117
        1.155
        1.206
        1.244
        1.288
        1.409
        1.505
        1.594
        1.664
        1.727
        1.816
        1.918
        2.000
        2.128
        2.216
        2.312
        2.471
        2.630
        2.820
        2.979
        3.170
        3.348
        3.539
        3.787
        4.073
        4.372
        4.664
        4.721
        4.880
        4.982
    ]

pVec_paper =
    [
        0.000
        4.728
        7.728
        9.455
        12.455
        15.728
        17.546
        18.638
        20.183
        20.821
        21.912
        23.277
        24.460
        24.916
        25.462
        26.645
        27.737
        27.465
        27.829
        28.467
        29.742
        30.289
        29.199
        28.109
        28.018
        28.564
        28.838
        28.748
        27.748
        27.386
        26.841
        26.842
        26.753
        26.300
        26.120
        25.394
        24.668
        24.397
        24.308
        23.673
        22.584
        20.768
        19.406
        18.954
        18.321
        17.324
        16.145
        15.694
        14.879
        14.155
        13.432
        12.529
        11.534
        10.721
        10.449
        10.180
        10.091]

# plot!(fig2, dVec_paper[1:length(dVec)], pVec_paper[1:length(dVec)])
plot!(fig2, dVec_paper[1:10], pVec_paper[1:10], label="Kusel-Kearlsey")

savefig(fig2, "$(figspath)ejemplo8Pd.png")


kappaVec_paper =
    [
        0.000
        0.003
        0.004
        0.004
        0.008
        0.013
        0.023
        0.029
        0.041
        0.052
        0.054
        0.061
        0.066
        0.081
        0.092
        0.108
        0.123
        0.135
        0.149
        0.162
        0.171
        0.194
        0.204
        0.222
        0.229
        0.238
        0.259
        0.278
        0.295
        0.324
        0.353
        0.376
        0.416
        0.439
        0.446
        0.470
        0.491
        0.525
        0.552
        0.572
        0.594
        0.630
        0.657
        0.680
        0.704
        0.732
        0.763
        0.794
        0.828
        0.849
        0.886
        0.925
        0.961
        0.981
        0.997]

mVec_paper =
    [
        0.005
        0.503
        0.779
        1.010
        1.078
        1.069
        1.141
        1.218
        1.286
        1.268
        1.223
        1.273
        1.318
        1.350
        1.422
        1.458
        1.499
        1.481
        1.449
        1.395
        1.345
        1.313
        1.286
        1.245
        1.191
        1.146
        1.123
        1.078
        1.064
        1.005
        0.956
        0.915
        0.851
        0.815
        0.788
        0.752
        0.729
        0.679
        0.661
        0.652
        0.616
        0.571
        0.548
        0.534
        0.530
        0.503
        0.480
        0.453
        0.430
        0.426
        0.403
        0.389
        0.385
        0.371
        0.376
    ]

plot!(fig, kappaVec_paper[1:10], mVec_paper[1:10], label="Kusel-Kearlsey")

savefig(fig, "$(figspath)ejemplo8Mk.png")

fig