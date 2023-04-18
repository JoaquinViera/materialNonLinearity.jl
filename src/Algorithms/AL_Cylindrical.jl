# Arc-Length 
# =====================================

struct ArcLength_Cylindrical <: AbstractAlgorithm
    tolk::Float64
    tolu::Float64
    tolf::Float64
    nTimes::Int64
    initialDeltaLambda::Float64
    arcLengthIncrem::Vector{Float64}
    controlDofs::Vector{Int64} # not used ATM
end

function step!(alg::ArcLength_Cylindrical, Uₖ, ModelSol, KTₖ, Fintk, time, U, dispIter, varFext, currδu, convδu, c, λₖ)

    Fextk = ModelSol.Fextk
    freeDofs = ModelSol.freeDofs

    Fext_red = view(Fextk, freeDofs)
    Fint_red = view(Fintk, freeDofs)
    varFext = view(varFext, freeDofs)
    KTₖ_red = view(KTₖ, freeDofs, freeDofs)

    r = [(Fext_red - Fint_red) varFext]

    deltas = KTₖ_red \ r

    δu⃰ = view(deltas, :, 1)
    δū = view(deltas, :, 2)

    arcLengthNorm = zeros(length(Uₖ))
    # arcLengthNorm .= 1
    # arcLengthNorm[1:3:end] .= 1 # θy dofs
    arcLengthNorm[2:3:end] .= 1 # uz dofs
    # arcLengthNorm[6] = 1
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
            # println(" $time")
            # δλ = sign((convδu' * (arcLengthNorm .* δū))) * Δl / sqrt(δū' * (arcLengthNorm .* δū))
            δλ = Δl / sqrt(δū' * (arcLengthNorm .* δū))
            # δλ = sign(convδu' * δū) * Δl / sqrt(δū' * δū)
        end
    else # Cylindrical constraint method - De Souza Neto
        # println(dispIter)
        a = δū' * (arcLengthNorm .* δū) # eq 4.117
        b = 2 * (currδu + δu⃰)' * (arcLengthNorm .* δū)
        c = (currδu + δu⃰)' * (arcLengthNorm .* (currδu + δu⃰)) - Δl^2
        d = b^2 - 4 * a * c
        d < 0 ? error("Negative discriminant at time $time, reduce arc length increment.") : nothing
        λsol = -b / (2 * a) .+ sqrt(d) / (2 * a) * [-1; 1]
        verif = [(currδu + δu⃰ + λsol[1] * δū)' * (arcLengthNorm .* currδu)
            (currδu + δu⃰ + λsol[2] * δū)' * (arcLengthNorm .* currδu)]

        # aa = δū' * δū # eq 4.117
        # bb = 2 * (currδu + δu⃰)' * δū
        # cc = (currδu + δu⃰)' * (currδu + δu⃰) - Δl^2
        # dd = bb^2 - 4 * aa * cc

        # λsol = -bb / (2 * aa) .+ sqrt(dd) / (2 * aa) * [-1; 1]

        # val1 = (currδu + δu⃰ + λsol[1] * δū)' * currδu
        # val2 = (currδu + δu⃰ + λsol[2] * δū)' * currδu
        # verif = [val1; val2]

        δλ = λsol[findall(x -> x == maximum(verif), verif)[1]]

    end

    δUₖ = δu⃰ + δλ * δū
    currδu = δUₖ + currδu

    Uk_red = view(Uₖ, freeDofs) + δUₖ
    U[freeDofs] = Uk_red

    return copy(U), δUₖ, δλ, currδu

end

function sets!(alg::ArcLength_Cylindrical, nnodes, ndofs, args...)
    λₖ = 0.0
    U = zeros(nnodes * ndofs)
    return λₖ, U, nothing
end