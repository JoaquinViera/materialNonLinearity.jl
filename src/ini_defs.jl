#
# Definition of initial variables to run
#

function ini_defs(mesh, boundaryConds, analysisSettings)

    ndofs = 2 # degrees of freedom per node
    nnodes = size(mesh.nodesMat, 1)

    # varFext
    varFext = zeros(ndofs * nnodes)
    # Loaded nodes
    loadedNodes = boundaryConds.nodalForceMatrix[:, 1]
    for i in 1:length(loadedNodes)
        dofs = nodes2dofs(loadedNodes[i], ndofs)
        varFext[dofs] = boundaryConds.nodalForceMatrix[i, 2:3]
    end

    # Solution parameters
    Uk = zeros(ndofs * nnodes)
    Fextk = zeros(ndofs * nnodes)
    Fintk = zeros(ndofs * nnodes)

    matUk = vcat(Uk, []) # Matrix to store disps
    matFext = vcat(Fextk, []) # Matrix to store applied external forces
    matFint = vcat(Fintk, []) # Matrix to store interal forces 

    # Iteration parameters
    nTimes = length(analysisSettings.loadFactors)
    stopCrit = []

    # Supports
    fixed_dofs = []
    suppNodes = boundaryConds.suppMatrix[:, 1]
    for i in 1:length(suppNodes)
        for j in 1:ndofs
            if boundaryConds.suppMatrix[i, j+1] == Inf
                fixed_dofs = [fixed_dofs; nodes2dofs(suppNodes[i], ndofs)[j]]
            end
        end
    end

    # Degrees of freedom
    fixed_dofs = unique(fixed_dofs)
    dofs_vec = Vector(1:nnodes*ndofs)
    free_dofs = copy(dofs_vec)
    deleteat!(free_dofs, fixed_dofs)


    δUk = zeros(length(free_dofs))

    # store struct
    modelStore = modelSol(Uk, δUk, Fextk, Fintk, matUk, matFext, matFint, free_dofs, [0.0])
    iterData = iterParams(nTimes, stopCrit)

    return modelStore, iterData, varFext
end




