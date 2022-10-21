# 
# Computes stress and tangent modulus
#

function constitutive_model(MaterialModel, εₖ)

    E = MaterialModel.E

    if cmp(MaterialModel.name, "linearElastic") == 0
        ∂σ∂ε = E
        σ = E * εₖ
    elseif cmp(MaterialModel.name, "isotropicBiLinear") == 0
        σY = MaterialModel.σY0
        σ_ₜᵣ = abs(E * εₖ)
        if σ_ₜᵣ >= σY
            K = MaterialModel.params[3]
            εY = σY / E
            ∂σ∂ε = E * K / (E + K)
            σ = σY * sign(εₖ) + ∂σ∂ε * (εₖ - εY * sign(εₖ))

            if ∂σ∂ε < 0 && sign(σ) != sign(σY * sign(εₖ))
                σ = 0
                ∂σ∂ε = 0
            end

        else
            ∂σ∂ε = E
            σ = E * εₖ
        end
    else
        error("MaterialModel model to be implemented")
    end

    return σ, ∂σ∂ε
end

