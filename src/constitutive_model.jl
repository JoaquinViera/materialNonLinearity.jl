# 
# Computes stress and tangent modulus
#

function constitutive_model(MaterialModel, εₖ)

    E = MaterialModel.E

    if cmp(MaterialModel.name, "linearElastic") == 0
        dsigdeps = E
        σ = E * εₖ
    elseif cmp(MaterialModel.name, "isotropicBiLinear") == 0
        σY = MaterialModel.σY0
        σ_tr = abs(E * εₖ)
        if σ_tr >= σY
            K = MaterialModel.params[3]
            εY = σY / E
            dsigdeps = E * K / (E + K)
            σ = σY * sign(εₖ) + dsigdeps * (εₖ - εY * sign(εₖ))
            if dsigdeps < 0 && sign(σ) != sign(σY * sign(εₖ))
                #error("Softening not available yet")
                σ = 0
                dsigdeps = 0
            end
        else
            dsigdeps = E
            σ = E * εₖ
        end
    else
        error("MaterialModel model to be implemented")
    end

    return σ, dsigdeps
end

