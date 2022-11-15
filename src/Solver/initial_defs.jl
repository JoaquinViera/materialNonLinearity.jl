#
# Definition of initial variables to run
#

function initial_defs(Mesh, BoundaryConds, AnalysisSettings)

    println(" \n==================================================")
    println("Starting analysis.")
    println("================================================== \n")

    ndofs = 2 # degrees of freedom per node
    nnodes = size(Mesh.nodesMat, 1)

    # Iteration parameters
    nTimes = AnalysisSettings.nTimes
    stopCrit = Vector{Int64}(undef, nTimes)

    # varFext
    varFext = zeros(ndofs * nnodes)
    # Loaded nodes
    loadedNodes = BoundaryConds.nodalForceMatrix[:, 1]
    for i in 1:length(loadedNodes)
        dofs = nodes2dofs(loadedNodes[i], ndofs)
        varFext[dofs] = BoundaryConds.nodalForceMatrix[i, 2:3]
    end

    # Solution parameters
    Uₖ = zeros(ndofs * nnodes)
    Fextk = zeros(ndofs * nnodes)
    Fintk = zeros(ndofs * nnodes)

    matUₖ = Vector{Vector{Float64}}(undef, nTimes) # Matrix to store disps
    matUₖ[1] = Uₖ

    #matFext = vcat(Fextk, []) # Matrix to store applied external forces
    matFext = Vector{Vector{Float64}}(undef, nTimes) # Matrix to store applied external forces
    matFext[1] = Fextk
    matFint = vcat(Fintk, []) # Matrix to store interal forces 



    # Supports
    fixed_dofs = []
    suppNodes = BoundaryConds.suppMatrix[:, 1]
    for i in 1:length(suppNodes)
        for j in 1:ndofs
            if BoundaryConds.suppMatrix[i, j+1] == Inf
                fixed_dofs = [fixed_dofs; nodes2dofs(suppNodes[i], ndofs)[j]]
            end
        end
    end

    # Degrees of freedom
    fixed_dofs = unique(fixed_dofs)
    dofs_vec = Vector(1:nnodes*ndofs)
    free_dofs = copy(dofs_vec)
    deleteat!(free_dofs, fixed_dofs)


    δUₖ = Vector{Vector{Float64}}(undef, nTimes)
    δUₖ[1] = zeros(length(free_dofs))


    loadFactors = Vector{Float64}(undef, nTimes)
    loadFactors[1] = 0.0

    λₖ, U, c = sets!(AnalysisSettings, nnodes, ndofs)

    # store struct
    ModelStore = ModelSol(Uₖ, δUₖ, Fextk, Fintk, matUₖ, matFext, matFint, free_dofs, loadFactors)
    IterData = IterParams(nTimes, stopCrit)

    return ModelStore, IterData, varFext, λₖ, U, c
end
