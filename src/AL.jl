function AL(Uk, modelSol, KTk, Fintk, time, analysisSettings, dispIter, varFext, currδu, convδu)



    Fextk = modelSol.Fextk
    freeDofs = modelSol.freeDofs

    Fext_red = Fextk[freeDofs]
    Fint_red = Fintk[freeDofs]
    varFext = varFext[freeDofs]
    KTkred = KTk[freeDofs, freeDofs]

    r = [Fext_red - Fint_red varFext]
    deltas = KTkred \ r

    δu_ast = deltas[:, 1]
    δu_bar = deltas[:, 2]

    incremArcLen = 1e-4
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
        if norm(Uk) == 0
            #δλ = analysisSettings.initialDeltaLambda
            δλ = initialDeltaLambda
        else
            δλ = sign((convδu' * (arcLengthNorm .* δu_bar))) * Δl / sqrt(δu_bar' * (arcLengthNorm .* δu_bar))
        end
    else # Jirasek method
        #controlDofs = freeDofs[end-1]
        controlDofs = 10
        scalingProjection = -1
        # controlDofs = analysisSettings.controlDofs # ajustar esto luego
        #scalingProjection = analysisSettings.sign
        c = zeros(length(Uk))
        c[controlDofs] = scalingProjection
        c = c[freeDofs]
        δλ = (Δl - c' * currδu - c' * δu_ast) / (c' * δu_bar)
    end

    δUk = δu_ast + δλ * δu_bar

    Uk[freeDofs] = Uk[freeDofs] + δUk

    return Uk, δUk, δλ

end