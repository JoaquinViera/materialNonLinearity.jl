# Check convergence

function convergence_check(Uₖ, δUₖ, Fext_red, Fintk, AnalysisSettings, dispIter, time)

    # Disps stop
    normUₖ = norm(Uₖ)
    normδUₖ = norm(δUₖ)

    # Forces stop

    if dispIter >= AnalysisSettings.tolk
        cond = 1 # iters
        convIter = 0
        error("Non convergence at step $time")
    elseif (normδUₖ < AnalysisSettings.tolu * normUₖ) || (norm(Fintk - Fext_red) < AnalysisSettings.tolf * norm(Fext_red))
        if normδUₖ < AnalysisSettings.tolu * normUₖ
            cond = 2 # disps
        else
            cond = 3 # forces
        end
        convIter = 1
    else
        cond = 0
        convIter = 0
    end

    return cond, convIter
end