# Arc-Length 
# =====================================

struct ArcLength <: AbstractAlgorithm
    tolk::Float64
    tolu::Float64
    tolf::Float64
    nTimes::Float64
    initialDeltaLambda::Float64
    arcLengthIncrem::Vector{Float64}
    #arcLengthIncrem::Float64
    controlDofs::Vector{Int64}
    scalingProjection::Float64
end

function step!(alg::ArcLength, Uₖ, ModelSol, KTₖ, Fintk, time, dispIter, varFext, currδu, convδu)

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

    incremArcLen = alg.arcLengthIncrem
    initialDeltaLambda = alg.initialDeltaLambda
    controlDofs = alg.controlDofs
    scalingProjection = alg.scalingProjection

    arcLengthNorm = zeros(length(freeDofs))
    arcLengthNorm[1:2:end] .= 1

    if length(incremArcLen) > 1
        Δl = incremArcLen[time]
    else
        Δl = incremArcLen[1]
    end

    if dispIter == 1 # Predictor
        if norm(Uₖ) == 0
            δλ = initialDeltaLambda
        else
            δλ = sign((convδu' * (arcLengthNorm .* δū))) * Δl / sqrt(δū' * (arcLengthNorm .* δū))
        end
    else # Jirasek method
        c = zeros(length(Uₖ))
        c[controlDofs] .= scalingProjection
        c = c[freeDofs]
        δλ = (Δl - c' * currδu - c' * δu⃰) / (c' * δū)
    end

    δUₖ = δu⃰ + δλ * δū

    Uₖ[freeDofs] = Uₖ[freeDofs] + δUₖ
    currδu = δUₖ + currδu

    return Uₖ, δUₖ, δλ, currδu

end