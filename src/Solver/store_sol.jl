function store_sol(time, ModelSol, IterData, Uₖ, δUₖ, λₖ, cond)

    ModelSol.matUk[time+1] = Uₖ
    ModelSol.convδu[time+1] = δUₖ
    ModelSol.loadFactors[time+1] = λₖ
    IterData.stopCrit[time+1] = cond

    return ModelSol, IterData
end