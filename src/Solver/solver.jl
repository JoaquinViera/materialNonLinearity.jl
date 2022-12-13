# =============================================
# Calls numerical method 
# =============================================



function solver(Section, MaterialModel, Mesh, BoundaryConds, AnalysisSettings, problemName, StressArraySets)

    # Initialize required variables
    ModelSol, IterData, varFext, λₖ, U, c = initial_defs(Mesh, BoundaryConds, AnalysisSettings, problemName)

    # Counters
    time = 1
    nTimes = IterData.nTimes

    pbar = Progress(nTimes, dt=0.25, barglyphs=BarGlyphs("[=> ]"), barlen=35, color=:cyan)
    aux = zeros(length(ModelSol.freeDofs))

    σArr = [[zeros(MaterialModel.ns) for _ in 1:nTimes] for _ in StressArraySets.elems]
    loadFactor = zeros(nTimes)

    while nTimes > time
        # Sets current disp Vector
        Uₖ = ModelSol.matUk[time]
        convδu = ModelSol.convδu[time]
        currδu = aux

        # increment external force
        time > 1 ? λₖ = view(ModelSol.loadFactors, time)[1] : λₖ = 0.0
        ModelSol.matFext[time] = ModelSol.Fextk

        ModelSol.Fextk, loadFactor = compute_Fext!(AnalysisSettings, varFext, 0.0, time, ModelSol.Fextk, loadFactor)

        # Iters
        dispIter = 0 # iter counter
        convIter = 0 # Convergence control parameter

        while convIter == 0
            # Updates displacement iter
            dispIter += 1
            # Computes Tangent Stiffness Matrix KTₖ & Internal Forces
            Fintₖ, σArr, matFint, KTₖ = assembler(Section, MaterialModel, Mesh, Uₖ, 1, σArr, time, StressArraySets, ModelSol.matFint)

            # Computes Uₖ & δUₖ
            Uₖ, δUₖ, λₖ, currδu = step!(AnalysisSettings, Uₖ, ModelSol, KTₖ, Fintₖ, time, U, dispIter, varFext, currδu, convδu, c, λₖ)

            # Computes Fintₖ at computed Uₖ
            Fintₖ, σArr, matFint = assembler(Section, MaterialModel, Mesh, Uₖ, 0, σArr, time, StressArraySets, ModelSol.matFint)

            # Computes Fext
            ModelSol.Fextk, loadFactor = compute_Fext!(AnalysisSettings, varFext, λₖ, time, ModelSol.Fextk, loadFactor)

            # Check convergence
            cond, convIter = convergence_check(Uₖ[ModelSol.freeDofs], δUₖ, ModelSol.Fextk[ModelSol.freeDofs], Fintₖ[ModelSol.freeDofs], AnalysisSettings, dispIter, time)

            # Stores results if convergence
            if convIter == 1
                ModelSol, IterData = store_sol(time, ModelSol, IterData, Uₖ, δUₖ, λₖ, cond)
            end

        end
        # Updates time
        time += 1
        next!(pbar)
        # time > nTimes ? ModelSol.loadFactors = loadFactor : nothing
    end

    ModelSol.loadFactors = loadFactor

    println("\n")

    println("End.")
    println("==================================================")

    println("\n\n")

    return ModelSol, time, IterData, σArr
end
