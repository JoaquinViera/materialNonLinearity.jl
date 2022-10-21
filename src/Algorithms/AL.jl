function AL(Uₖ, ModelSol, KTₖ, Fintk, time, analysisSettings, dispIter, varFext, currδu, convδu)


    Fextk = ModelSol.Fextk
    freeDofs = ModelSol.freeDofs

    Fext_red = Fextk[freeDofs]
    Fint_red = Fintk[freeDofs]
    varFext = varFext[freeDofs]
    KTₖ_red = KTₖ[freeDofs, freeDofs]

    r = [Fext_red - Fint_red varFext]
    deltas = KTₖ_red \ r

    δu⃰ = deltas[:, 1]
    δū = deltas[:, 2]

    incremArcLen = 1e-5

    initialDeltaLambda = 1e-2
    arcLengthNorm = zeros(length(freeDofs))
    arcLengthNorm[1:2:end] .= 1

    #if length(analysisSettings.incremArcLen) > 1
    if length(incremArcLen) > 1
        #Δl = (analysisSettings.incremArcLen)[time]
        Δl = incremArcLen[time]
    else
        #Δl = analysisSettings.incremArcLen
        Δl = incremArcLen
    end

    if dispIter == 1 # Predictor
        if norm(Uₖ) == 0
            #δλ = analysisSettings.initialDeltaLambda
            δλ = initialDeltaLambda
        else
            δλ = sign((convδu' * (arcLengthNorm .* δū))) * Δl / sqrt(δū' * (arcLengthNorm .* δū))
        end
    else # Jirasek method
        #controlDofs = freeDofs[end-1]
        controlDofs = 10
        scalingProjection = 1
        # controlDofs = analysisSettings.controlDofs # ajustar esto luego
        #scalingProjection = analysisSettings.sign
        c = zeros(length(Uₖ))
        c[controlDofs] = scalingProjection
        c = c[freeDofs]
        δλ = (Δl - c' * currδu - c' * δu⃰) / (c' * δū)
    end

    δUₖ = δu⃰ + δλ * δū

    Uₖ[freeDofs] = Uₖ[freeDofs] + δUₖ

    return Uₖ, δUₖ, δλ

end