# Check convergence

function convergence_check(freeDofs, Uk, δUk, Fextk, Fintk, analysisSettings, dispIter, time)

    Fext_red = Fextk[freeDofs]

    # Disps stop
    normUk = norm(Uk[freeDofs])
    normDeltaUk = norm(δUk)

    # Forces stop
    norm_r = norm(Fintk[freeDofs] - Fext_red)
    normFext = norm(Fext_red)

    if dispIter >= analysisSettings.tolk
        cond = 1 # iters
        convIter = 0
        error("Non convergence at step $time")
    elseif (normDeltaUk < analysisSettings.tolu * normUk) || (norm_r < analysisSettings.tolf * normFext)
        if normDeltaUk < analysisSettings.tolu * normUk
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