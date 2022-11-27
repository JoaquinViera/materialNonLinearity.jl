function compute_Fext!(alg::ArcLength, varFext, currλ, time, Fextk, args...)
    Fext = varFext * currλ + Fextk
    return Fext
end

function compute_Fext!(alg::NewtonRaphson, varFext, currλ, time, Fextk, args...)

    Fext = varFext * alg.loadFactors[time]
    return Fext
end

