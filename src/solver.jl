# =============================================
# Calls numerical method 
# =============================================

function solver(Section, MaterialModel, Mesh, BoundaryConds, AnalysisSettings)
    # Initialize required variables
    ModelSol, IterData, varFext = initial_defs(Mesh, BoundaryConds, AnalysisSettings)

    # Counters
    time = 1
    nTimes = IterData.nTimes

    while nTimes > time
        println(time)

        # Sets current disp Vector
        Uk = ModelSol.matUk[:, time]
        convδu = ModelSol.convδu[:, time]
        currδu = zeros(length(ModelSol.freeDofs))
        #println(size(currδu))
        NRBool = 0
        # increment external force
        if NRBool == 1
            λk = AnalysisSettings.loadFactors[time]
        else # AL
            if time == 1
                λk = 0
            else
                λk = ModelSol.loadFactors[time]
            end
        end

        ModelSol.Fextk = ModelSol.Fextk + λk * varFext
        ModelSol.matFext = hcat(ModelSol.matFext, ModelSol.Fextk)

        # Iters
        dispIter = 0 # iter counter
        convIter = 0 # Convergence control parameter


        while convIter == 0
            # Computes Tangent Stiffness Matrix KTk & Internal Forces
            intBool = 1 # Boolean to control computations of variables
            Fintk, KTk = assembler(Section, MaterialModel, Mesh, Uk, intBool)

            # Computes Uk & δUk

            if NRBool == 1
                Uk, δUk = NR(Uk, ModelSol, KTk, Fintk)
                #λk = AnalysisSettings.loadFactors[time]    
            else
                Uk, δUk, λk = AL(Uk, ModelSol, KTk, Fintk, time, AnalysisSettings, dispIter, varFext, currδu, convδu)
                currδu = δUk + currδu
            end

            # Computes Fintk at computed Uk
            intBool = 1
            Fintk, KTk = assembler(Section, MaterialModel, Mesh, Uk, intBool)

            # Check convergence
            cond, convIter = convergence_check(ModelSol.freeDofs, Uk, δUk, ModelSol.Fextk, Fintk, AnalysisSettings, dispIter, time)

            # Stores results if convergence
            if convIter == 1
                #println("cond $cond")
                #neigs = sum(eigvals(KTk) .< 0)
                #println("number of negative eigs $neigs")
                ModelSol.loadFactors = hcat(ModelSol.loadFactors, λk)
                ModelSol.convδu = hcat(ModelSol.convδu, δUk)
                ModelSol.matUk = hcat(ModelSol.matUk, Uk)
                ModelSol.matFint = hcat(ModelSol.matFint, Fintk)
                IterData.stopCrit = vcat(IterData.stopCrit, cond)
            end

            # Updates disp Iter
            dispIter += 1
        end
        # Updates time
        time += 1
    end

    return ModelSol, time, IterData
end