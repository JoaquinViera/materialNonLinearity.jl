# Check convergence

function convergence_check(freeDofs, Uₖ, δUₖ, Fextk, Fintk, AnalysisSettings, dispIter, time)

    Fext_red = Fextk[freeDofs]

    # Disps stop
    normUₖ = norm(Uₖ[freeDofs])
    normδUₖ = norm(δUₖ)

    # Forces stop
    normResidual = norm(Fintk[freeDofs] - Fext_red)
    normFext = norm(Fext_red)

    if dispIter >= AnalysisSettings.tolk
        cond = 1 # iters
        convIter = 0
        error("Non convergence at step $time")
    elseif (normδUₖ < AnalysisSettings.tolu * normUₖ) || (normResidual < AnalysisSettings.tolf * normFext)
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