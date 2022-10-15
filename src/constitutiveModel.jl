# 
# Computes stress and tangent modulus
#

function constitutiveModel(material, epsk)

    E = material.E

    if cmp(material.name, "linearElastic") == 0
        dsigdeps = E
        sigma = E * epsk
    elseif cmp(material.name, "isotropicBiLinear") == 0
        sigmaY = material.σY0
        sigma_tr = abs(E * epsk)
        if sigma_tr >= sigmaY
            K = material.params[3]
            epsY = sigmaY / E
            dsigdeps = E * K / (E + K)
            sigma = sigmaY * sign(epsk) + dsigdeps * (epsk - epsY * sign(epsk))
            if dsigdeps < 0
                error("Softening not available yet")
            end
        else
            dsigdeps = E
            sigma = E * epsk
        end
    else
        error("Material model to be implemented")
    end

    return sigma, dsigdeps
end

