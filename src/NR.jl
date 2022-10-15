# Newton-Raphson 
# =====================================

function NR(Uk, modelSol, KTk, Fintk)

    Fextk = modelSol.Fextk
    freeDofs = modelSol.freeDofs

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