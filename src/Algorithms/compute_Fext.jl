function compute_Fext!(alg::ArcLength, varFext, currλ, time, Fextk, AL_factor, args...)
    Fext = varFext * currλ + Fextk
    loadFactor = AL_factor + currλ
    return Fext, loadFactor
end

function compute_Fext!(alg::NewtonRaphson, varFext, currλ, time, Fextk, args...)
    Fext = varFext * alg.loadFactors[time]
    return Fext, alg.loadFactors[time]
end

