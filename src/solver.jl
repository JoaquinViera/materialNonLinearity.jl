# =============================================
# Calls numerical method 
# =============================================

function solver(strSections, strMaterialModels, strMesh, strBC, strAnalysisSets)
    # Initialize required variables
    modelSol, iterData = ini_defs(strMesh, strBC, strAnalysisSets)

    # Counters
    time = 1
    nTimes = iterData.nTimes
    # Loaded nodes
    loadedNodes = strBC.nodalForceMatrix[:, 1]
    ndofs = 2
    while nTimes > time
        println(time)
        # Sets current disp Vector

        Uk = modelSol.matUk[:, time]

        # increment external force
        λk = strAnalysisSets.loadFactors[time]
        for i in 1:length(loadedNodes)
            dofs = nodes2dofs(loadedNodes[i], ndofs)
            modelSol.Fextk[dofs] = modelSol.Fextk[dofs] + λk * strBC.nodalForceMatrix[i, 2:3]
        end
        modelSol.matFext = hcat(modelSol.matFext, modelSol.Fextk)



        # NR iter
        dispIter = 0 # iter counter
        convIter = 0 # Convergence control parameter

        while convIter == 0
            Uk, Fintk, FintkL, cond, convIter = NR(strSections, strMaterialModels, strMesh, Uk, modelSol, time, strAnalysisSets, dispIter)
            dispIter += 1
            if convIter == 1
                # Stores results
                modelSol.matUk = hcat(modelSol.matUk, Uk)
                modelSol.matFint = hcat(modelSol.matFint, Fintk)
                iterData.stopCrit = vcat(iterData.stopCrit, cond)
            end
        end
        # Updates time
        time += 1

    end

    return modelSol, time
end