# Check convergence

function convergence_check(Uₖ, δUₖ, Fext_red, Fint_red, AnalysisSettings, dispIter, time)

    # Disps stop
    normUₖ = norm(Uₖ)
    normδUₖ = norm(δUₖ)

    # Forces stop
    if dispIter >= AnalysisSettings.tolk
        cond = 1 # iters
        convIter = 0
        error("Non convergence at step $time")
    elseif (normδUₖ <= AnalysisSettings.tolu * normUₖ) || (norm(Fint_red - Fext_red) <= AnalysisSettings.tolf * norm(Fext_red))
        if norm(Fint_red - Fext_red) <= AnalysisSettings.tolf * norm(Fext_red)
            cond = 3 # forces
        else
            cond = 2 # disps
        end
        convIter = 1
    else
        cond = 0
        convIter = 0
    end

    return cond, convIter
end