# ===============================================
# Cantilever with elasto plastic material model 
# ===============================================

# Load solver module
using materialNonLinearity, Plots, LinearAlgebra, FastGaussQuadrature

# example name
problemName = "Concrete"

# Define material model
# =======================================
E = 28e6
σY = 3e3
K = -E
ne = 20
ns = 20



import materialNonLinearity: constitutive_model

# Materials struct
StrMaterialModels = UserModel(ne, ns)


function constitutive_model(ElemMaterialModel::UserModel, εₖ)
    fck = 30 # MPa
    E = 28e6 # kN/m2
    fctmfl = 3e3 # kN/m2
    yc = 1.5
    fcd = fck / yc

    # Tension
    eps1 = fctmfl / E
    K = -E

    # Compression
    if fck <= 50
        epsc0 = 2e-3
        epscu = 3.5e-3
        n = 2
    else
        epsc0 = 2e-3 + 0.000085 * (fck - 50)^(0.50)
        epscu = 0.0026 + 0.0144 * ((100 - fck) / 100)^4
        n = 1.4 + 9.6 * ((100 - fck) / 100)^4
    end

    if εₖ >= 0.0
        # Tension
        if εₖ <= eps1
            σ = E * εₖ
            ∂σ∂ε = E
        else
            if εₖ >= eps1 * (1 - E / K)
                σ = 0.0
                ∂σ∂ε = 0.0
                # println("n")
            else
                σ = fctmfl + K * (εₖ - eps1)
                ∂σ∂ε = K
            end
        end
        #Compression
    else
        #epsc = abs(εₖ)
        #=
                if epsc <= epsc0
                    σ = -fcd * (1 - (1 - epsc / epsc0) .^ n) * 1000
                    ∂σ∂ε = fcd * n * (1 - epsc / epsc0) .^ (n - 1) / epsc0 * 1000
                else
                    σ = -fcd * 1000
                    ∂σ∂ε = 0
                end
        =#
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
initialDeltaLambda = 1e-5 #
#arcLengthIncrem = vcat(ones(6) * 1e-4, ones(17) * 1e-5, ones(75) * 1e-6)  # compresion lineal
arcLengthIncrem = vcat(ones(6) * 1e-4, ones(4) * 1e-5, ones(30) * 2e-6)  # compresion polinomica
nLoadSteps = length(arcLengthIncrem)
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

# ===============================================
# Process model parameters
# ===============================================

sol, time, IterData, σArr = solver(StrSections, StrMaterialModels, StrMesh, StrBoundaryConds, StrAnalysisSettings, problemName)

println(IterData.stopCrit)

# Post process
# --------------------------------
P = abs(Fz)
Iy = StrSections.Iy
Mfis = σY * Iy / (h / 2)
println(Mfis)

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
println(mVec[end])
dVec = hcat([i[dofD] for i in matUk])
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
#=
Mana = zeros(nLoadSteps)
C = K
epsY = σY / E
kappae = 2 * σY / (E * h)
#eps_ast = epsY - σY / C
#kappa_ast = 2 * eps_ast / h
elem = 1
for i in 1:nLoadSteps
    kappak = kappaHistElem[elem, i]
    if kappak <= kappae
        Mana[i] = E * StrSections.Iy * kappak
    else
        Msup = E * StrSections.Iy * kappak / 2
        Minf = σY * b * h^2 / 24 * (3 - kappae^2 / kappak^2) + σY * K * b * h^2 * kappak / kappae / (12 * E) * (1 - 3 / 2 * kappae / kappak + (kappae / kappak)^3 / 2)
        Mana[i] = Msup + Minf
    end
end

err = (abs.(mVec[2:end]) - Mana[2:end]) ./ Mana[2:end] * 100
maxErrMk = maximum(err)
println(maxErrMk)
=#
# M-κ plot  
# --------------------------------
elem = 1
fig = plot(kappaHistElem[elem, :], abs.(mVec), markershape=:circle, lw=lw, ms=ms, title="M-κ", label="FEM", minorgrid=1, draw_arrow=1, legend=:bottomright)
#plot!(fig, kappaHistElem[elem, :], Mana, markershape=:rect, lw=lw, ms=ms, label="Analytic")
xlabel!("κ")
ylabel!("M")

# P-δ plot  
# --------------------------------
fig2 = plot(abs.(dVec), pVec, markershape=:circle, lw=lw, ms=ms, title="P-δ", label="FEM", minorgrid=1, draw_arrow=1, legend=:bottomright)
xlabel!("δ")
ylabel!("P")

p, w = gausslegendre(ns)


sfig = plot(σArr[end], p * h / 2, markershape=:circle, lw=lw, ms=ms, title="stress", label="stress", minorgrid=1, draw_arrow=1)
plot!(sfig, zeros(length(p)), p * h / 2, lw=lw, ms=ms, label="z", color=:"black")
