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

function step!(alg::NewtonRaphson, Uk, ModelSol, KTk, Fintk, time, U, args...)

    Fextk = ModelSol.Fextk
    freeDofs = ModelSol.freeDofs

    # Solve system
    KTₖ_red = view(KTk, freeDofs, freeDofs)
    Fext_red = view(Fextk, freeDofs)
    Fint_red = view(Fintk, freeDofs)

    δUₖ = KTₖ_red \ (Fext_red - Fint_red)

    # Computes Uk
    Uk_red = view(Uk, freeDofs) + δUₖ
    U[freeDofs] = Uk_red

    # Computes load factor
    λₖ = view(alg.loadFactors, time)[1]

    return copy(U), δUₖ, λₖ, nothing

end

function sets!(alg::NewtonRaphson, nnodes, ndofs, args...)
    λₖ = alg.loadFactors[1]
    U = zeros(nnodes * ndofs)
    return λₖ, U, nothing
end