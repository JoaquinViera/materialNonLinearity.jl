# =============================================
# Calls numerical method 
# =============================================

function solver(strSections, strMaterialModels, strMesh, strBC, strAnalysisSets)
    # Initialize required variables
    modelSol, iterData, varFext = ini_defs(strMesh, strBC, strAnalysisSets)

    # Counters
    time = 1
    nTimes = iterData.nTimes

    while nTimes > time
        println(time)

        # Sets current disp Vector
        Uk = modelSol.matUk[:, time]
        convδu = modelSol.convδu[:, time]
        currδu = zeros(length(modelSol.freeDofs))
        #println(size(currδu))
        NRBool = 0
        # increment external force
        if NRBool == 1
            λk = strAnalysisSets.loadFactors[time]
        else # AL
            if time == 1
                λk = 0
            else
                λk = modelSol.loadFactors[time]
                println("lambda")
                #println(λk)
            end
        end

        modelSol.Fextk = modelSol.Fextk + λk * varFext
        modelSol.matFext = hcat(modelSol.matFext, modelSol.Fextk)

        # Iters
        dispIter = 0 # iter counter
        convIter = 0 # Convergence control parameter


        while convIter == 0
            # Computes Tangent Stiffness Matrix KTk & Internal Forces
            intBool = 1 # Boolean to control computations of variables
            Fintk, KTk = assembler(strSections, strMaterialModels, strMesh, Uk, modelSol, time, strAnalysisSets, dispIter, intBool)

            # Computes Uk & δUk

            if NRBool == 1
                Uk, δUk = NR(Uk, modelSol, KTk, Fintk)
                #λk = strAnalysisSets.loadFactors[time]    
            else
                Uk, δUk, λk = AL(Uk, modelSol, KTk, Fintk, time, analysisSettings, dispIter, varFext, currδu, convδu)
                currδu = δUk + currδu
                println(λk)
            end

            # Computes Fintk at computed Uk
            intBool = 0
            Fintk = assembler(strSections, strMaterialModels, strMesh, Uk, modelSol, time, strAnalysisSets, dispIter, intBool)

            # Check convergence
            cond, convIter = convergenceCheck(modelSol.freeDofs, Uk, δUk, modelSol.Fextk, Fintk, strAnalysisSets, dispIter)

            # Stores results if convergence
            if convIter == 1
                modelSol.loadFactors = hcat(modelSol.loadFactors, λk)
                modelSol.convδu = hcat(modelSol.convδu, δUk)
                modelSol.matUk = hcat(modelSol.matUk, Uk)
                modelSol.matFint = hcat(modelSol.matFint, Fintk)
                iterData.stopCrit = vcat(iterData.stopCrit, cond)
            end

            # Updates disp Iter
            dispIter += 1
        end
        # Updates time
        time += 1
    end

    return modelSol, time
end