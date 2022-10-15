function AL(section, material, mesh, Uk, modelSol, time, analysisSettings, dispIter)

    deltas = KT \ residual
    δu_ast = deltas[1]
    δu_bar = deltas[2]

    if length(analysisSettings.incremArcLen) > 1
        Δl = (analysisSettings.incremArcLen)(time)
    else
        Δl = analysisSettings.incremArcLen
    end

    if dispIter == 1 # Predictor
        if norm(Uk) == 0
            δλ = analysisSettings.initialDeltaLAmbda
        else
            δλ = sign((convUk' * (arcLengthNorm .* δu_bar))) * Δl / sqrt(δu_bar' * (arcLengthNorm .* δu_bar))
        end
    else # Jirasek method
        controlDofs = analysisSettings.controlDofs # ajustar esto luego
        scalingProjection = analysisSettings.sign
        c = zeros(size(Uk), 1)
        c[controlDofs] = scalingProjection
        δλ = (Δl - c' * currδu - c' * δu_ast) / (c' * δu_bar)
    end

    δUk = δu_ast + δλ * δu_bar

    Uk[freeDofs] = Uk[freeDofs] + δUk

    return Uk, Fintk, FintkL, cond, convParam

end