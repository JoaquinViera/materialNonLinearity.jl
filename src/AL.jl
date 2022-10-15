function AL(section, material, mesh, Uk, modelSol, time, analysisSettings, dispIter)

    #=
        deltas = KT \ residual
        δu_ast = deltas[1]
        δu_bar = deltas[2]

        if length(analysisSettings.incremArcLen) > 1
            Δl = (analysisSettings.incremArcLen)(time)
        else
            Δl = analysisSettings.incremArcLen
        end


        if dispIter == 1 # Predictor
            # todo...
        end

        # Jirasek method
        controlDofs = analysisSettings.controlDofs # ajustar esto luego
        projection = analysisSettings.sign
        c = zeros(size(Uk), 1)
        c[controlDofs] = projection
        δλ = (Δl - c'* )

        Uk = 0
        Fintk = 0
        FintkL = 0
        cond = 0
        convParam = 0

        return Uk, Fintk, FintkL, cond, convParam
    =#
end