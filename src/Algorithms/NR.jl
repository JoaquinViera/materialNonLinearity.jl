# Newton-Raphson 
# =====================================

struct NewtonRaphson <: AbstractAlgorithm
    tolk::Float64
    tolu::Float64
    tolf::Float64
    loadFactors::Vector{Float64}
    nTimes::Int64
    function NewtonRaphson(tolk, tolu, tolf, loadFactors)
        nTimes = length(loadFactors)
        return new(tolk, tolu, tolf, loadFactors, nTimes)
    end
end

function step!(alg::NewtonRaphson, Uk, ModelSol, KTk, Fintk, time, args...)

    Fextk = ModelSol.Fextk
    freeDofs = ModelSol.freeDofs

    # Solve system
    KTkred = KTk[freeDofs, freeDofs]
    Fext_red = Fextk[freeDofs]
    Fint_red = Fintk[freeDofs]

    r = Fext_red - Fint_red

    δUk = KTkred \ r

    # Computes Uk
    Uk[freeDofs] = Uk[freeDofs] + δUk

    # Computes load factor
    λₖ = alg.loadFactors[time]

    return Uk, δUk, λₖ, nothing

end
