# Arc-Length 
# =====================================

struct ArcLength <: AbstractAlgorithm
    tolk::Float64
    tolu::Float64
    tolf::Float64
    nTimes::Int64
    initialDeltaLambda::Float64
    arcLengthIncrem::Vector{Float64}
    controlDofs::Vector{Int64}
    scalingProjection::Float64
end

function step!(alg::ArcLength, Uₖ, ModelSol, KTₖ, Fintk, time, U, dispIter, varFext, currδu, convδu, c, λₖ)

    Fextk = ModelSol.Fextk
    freeDofs = ModelSol.freeDofs

    Fext_red = view(Fextk, freeDofs)
    Fint_red = view(Fintk, freeDofs)
    varFext = view(varFext, freeDofs)
    KTₖ_red = view(KTₖ, freeDofs, freeDofs)

    r = [Fext_red - Fint_red varFext]

    deltas = KTₖ_red \ r

    δu⃰ = view(deltas, :, 1)
    δū = view(deltas, :, 2)


    arcLengthNorm = zeros(length(Uₖ))
    arcLengthNorm[1:3:end] .= 1
    arcLengthNorm[2:3:end] .= 1
    arcLengthNorm = view(arcLengthNorm, freeDofs)

    if length(alg.arcLengthIncrem) > 1
        Δl = view(alg.arcLengthIncrem, time)[1]
    else
        Δl = view(alg.arcLengthIncrem, 1)[1]
    end

    if dispIter == 1 # Predictor
        if norm(convδu) == 0
            δλ = alg.initialDeltaLambda
        else
            δλ = sign((convδu' * (arcLengthNorm .* δū))) * Δl / sqrt(δū' * (arcLengthNorm .* δū))
        end
    else # Jirasek method
        c = view(c, freeDofs)
        δλ = (Δl - c' * currδu - c' * δu⃰) / (c' * δū)
    end

    δUₖ = δu⃰ + δλ * δū
    currδu = δUₖ + currδu

    Uk_red = view(Uₖ, freeDofs) + δUₖ
    U[freeDofs] = Uk_red

    return copy(U), δUₖ, δλ, currδu

end

function sets!(alg::ArcLength, nnodes, ndofs, args...)
    λₖ = 0.0
    c = zeros(nnodes * ndofs)
    c[alg.controlDofs] .= alg.scalingProjection
    U = zeros(nnodes * ndofs)
    return λₖ, U, c
end