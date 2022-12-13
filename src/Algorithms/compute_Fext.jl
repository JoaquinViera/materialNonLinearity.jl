function compute_Fext!(alg::ArcLength, varFext, currλ, time, Fextk, loadFactor, args...)
    Fext = varFext * currλ + Fextk
    loadFactor[time+1] = currλ + loadFactor[time]
    return Fext, loadFactor
end

function compute_Fext!(alg::NewtonRaphson, varFext, currλ, time, Fextk, loadFactor, args...)

    Fext = varFext * alg.loadFactors[time]
    loadFactor[time+1] = alg.loadFactors[time]
    return Fext, loadFactor
end

