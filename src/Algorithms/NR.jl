# Newton-Raphson 
# =====================================

function NR(Uk, ModelSol, KTk, Fintk)

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

    return Uk, δUk

end