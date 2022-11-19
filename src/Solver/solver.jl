# =============================================
# Calls numerical method 
# =============================================



function solver(Section, MaterialModel, Mesh, BoundaryConds, AnalysisSettings, problemName)

    # Initialize required variables
    ModelSol, IterData, varFext, λₖ, U, c = initial_defs(Mesh, BoundaryConds, AnalysisSettings, problemName)

    # Counters
    time = 1
    nTimes = IterData.nTimes

    pbar = Progress(nTimes, dt=0.25, barglyphs=BarGlyphs("[=> ]"), barlen=35, color=:cyan)
    aux = zeros(length(ModelSol.freeDofs))

    while nTimes > time
        #println(time)
        # Sets current disp Vector
        Uₖ = ModelSol.matUk[time]
        convδu = ModelSol.convδu[time]
        currδu = aux

        # increment external force
        time > 1 ? λₖ = view(ModelSol.loadFactors, time)[1] : nothing

        ModelSol.Fextk = ModelSol.Fextk + λₖ * varFext
        ModelSol.matFext[time] = ModelSol.Fextk

        # Iters
        dispIter = 0 # iter counter
        convIter = 0 # Convergence control parameter

        while convIter == 0
            # Computes Tangent Stiffness Matrix KTₖ & Internal Forces
            Fintₖ, KTₖ = assembler(Section, MaterialModel, Mesh, Uₖ, 1)

            # Computes Uₖ & δUₖ
            Uₖ, δUₖ, λₖ, currδu = step!(AnalysisSettings, Uₖ, ModelSol, KTₖ, Fintₖ, time, U, dispIter, varFext, currδu, convδu, c)

            # Computes Fintₖ at computed Uₖ
            Fintₖ = assembler(Section, MaterialModel, Mesh, Uₖ, 0)

            # Check convergence
            cond, convIter = convergence_check(Uₖ[ModelSol.freeDofs], δUₖ, ModelSol.Fextk[ModelSol.freeDofs], Fintₖ[ModelSol.freeDofs], AnalysisSettings, dispIter, time)

            # Stores results if convergence
            if convIter == 1
                ModelSol, IterData = store_sol(time, ModelSol, IterData, Uₖ, δUₖ, Fintₖ, λₖ, cond)
            end

            # Updates displacement iter
            dispIter += 1
        end
        # Updates time
        time += 1
        next!(pbar)

    end

    println("\n")

    println("End.")
    println("==================================================")

    println("\n\n")

    return ModelSol, time, IterData
end
