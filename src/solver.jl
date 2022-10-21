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
        Uₖ = ModelSol.matUk[:, time]
        convδu = ModelSol.convδu[:, time]
        currδu = zeros(length(ModelSol.freeDofs))
        #println(size(currδu))
        NRBool = 1
        # increment external force
        if NRBool == 1
            λₖ = AnalysisSettings.loadFactors[time]
        else # AL
            if time == 1
                λₖ = 0
            else
                λₖ = ModelSol.loadFactors[time]
            end
        end

        ModelSol.Fextk = ModelSol.Fextk + λₖ * varFext

        ModelSol.matFext = hcat(ModelSol.matFext, ModelSol.Fextk)

        # Iters
        dispIter = 0 # iter counter
        convIter = 0 # Convergence control parameter


        while convIter == 0
            # Computes Tangent Stiffness Matrix KTₖ & Internal Forces
            intBool = 1 # Boolean to control computations of variables
            Fintₖ, KTₖ = assembler(Section, MaterialModel, Mesh, Uₖ, intBool)

            # Computes Uₖ & δUₖ

            if NRBool == 1
                Uₖ, δUₖ = NR(Uₖ, ModelSol, KTₖ, Fintₖ)
                #λₖ = AnalysisSettings.loadFactors[time]    
            else
                Uₖ, δUₖ, λₖ = AL(Uₖ, ModelSol, KTₖ, Fintₖ, time, AnalysisSettings, dispIter, varFext, currδu, convδu)
                currδu = δUₖ + currδu
            end

            # Computes Fintₖ at computed Uₖ
            intBool = 1
            Fintₖ, KTₖ = assembler(Section, MaterialModel, Mesh, Uₖ, intBool)

            # Check convergence
            cond, convIter = convergence_check(ModelSol.freeDofs, Uₖ, δUₖ, ModelSol.Fextk, Fintₖ, AnalysisSettings, dispIter, time)

            # Stores results if convergence
            if convIter == 1
                neigs = sum(eigvals(KTₖ[ModelSol.freeDofs, ModelSol.freeDofs]) .< 0)
                println("number of negative eigs $neigs")
                ModelSol.loadFactors = hcat(ModelSol.loadFactors, λₖ)
                ModelSol.convδu = hcat(ModelSol.convδu, δUₖ)
                ModelSol.matUk = hcat(ModelSol.matUk, Uₖ)
                ModelSol.matFint = hcat(ModelSol.matFint, Fintₖ)
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