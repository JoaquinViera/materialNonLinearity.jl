# 
# Computes stress and tangent modulus
#

function constitutive_model(ElemMaterialModel::LinearElastic, εₖ)

    E = ElemMaterialModel.E
    ∂σ∂ε = E
    σ = E * εₖ

    return σ, ∂σ∂ε
end



function constitutive_model(ElemMaterialModel::IsotropicBiLinear, εₖ)

    E = ElemMaterialModel.E
    σY = ElemMaterialModel.σY0
    σ_ₜᵣ = abs(E * εₖ)
    if σ_ₜᵣ >= σY

        K = ElemMaterialModel.K
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

    return σ, ∂σ∂ε
end
