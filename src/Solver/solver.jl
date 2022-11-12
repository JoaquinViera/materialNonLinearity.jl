# =============================================
# Calls numerical method 
# =============================================



function solver(Section, MaterialModel, Mesh, BoundaryConds, AnalysisSettings)
    # Initialize required variables
    ModelSol, IterData, varFext, λₖ, U, c = initial_defs(Mesh, BoundaryConds, AnalysisSettings)

    # Counters
    time = 1
    nTimes = IterData.nTimes

    # Progress print
    progressFrame = round((nTimes - 1) / 5)
    if progressFrame == 0
        counter = 5
    else
        counter = 0
    end
    #println(ModelSol.matUk)
    while nTimes > time
        #println(time)

        if time == counter * progressFrame + 1
            println("$(counter * 20) %")
            counter = counter + 1
        end
        # Sets current disp Vector
        #Uₖ = ModelSol.matUk[:, time]
        Uₖ = ModelSol.matUk[time]

        convδu = ModelSol.convδu[time]
        #convδu = view(ModelSol.convδu, time)[1]

        currδu = zeros(length(ModelSol.freeDofs))

        # increment external force
        time > 1 ? view(ModelSol.loadFactors, time)[1] : nothing

        ModelSol.Fextk = ModelSol.Fextk + λₖ * varFext
        ModelSol.matFext = hcat(ModelSol.matFext, ModelSol.Fextk)

        # Iters
        dispIter = 0 # iter counter
        convIter = 0 # Convergence control parameter

        while convIter == 0
            # Computes Tangent Stiffness Matrix KTₖ & Internal Forces
            # intBool = 1 # Boolean to control computations of variables
            Fintₖ, KTₖ = assembler(Section, MaterialModel, Mesh, Uₖ, 1)

            # Computes Uₖ & δUₖ
            Uₖ, δUₖ, λₖ, currδu = step!(AnalysisSettings, Uₖ, ModelSol, KTₖ, Fintₖ, time, U, dispIter, varFext, currδu, convδu, c)
            # Computes Fintₖ at computed Uₖ
            #intBool = 1
            Fintₖ, KTₖ = assembler(Section, MaterialModel, Mesh, Uₖ, 1)

            # Check convergence
            cond, convIter = convergence_check(Uₖ[ModelSol.freeDofs], δUₖ, ModelSol.Fextk[ModelSol.freeDofs], Fintₖ[ModelSol.freeDofs], AnalysisSettings, dispIter, time)

            # Stores results if convergence
            if convIter == 1
                #rintln(ModelSol.matUk)
                neigs = sum(eigvals(KTₖ[ModelSol.freeDofs, ModelSol.freeDofs]) .< 0)
                #println("number of negative eigs $neigs")
                push!(ModelSol.loadFactors, λₖ)
                #ModelSol.convδu = hcat(ModelSol.convδu, δUₖ)
                #println(δUₖ)
                push!(ModelSol.convδu, δUₖ)
                #ModelSol.matUk = hcat(ModelSol.matUk, Uₖ)
                #println(Uₖ)
                #println(ModelSol.matUk)
                #println(ModelSol.matUk[time])
                #push!(ModelSol.matUk, Uₖ)
                #println(ModelSol.matUk[time])
                ModelSol.matUk[time+1] = Uₖ
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
