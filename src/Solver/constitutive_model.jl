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


function constitutive_model(ElemMaterialModel::UserModel, εₖ)
    fck = 30 # MPa
    E = 28e6 # kN/m2
    fctmfl = 3e3 # kN/m2
    yc = 1.5
    fcd = fck / yc

    # Tension
    eps1 = fctmfl / E
    K = -E / 10

    # Compression
    if fck <= 50
        epsc0 = 2e-3
        epscu = 3.5e-3
        n = 2
    else
        epsc0 = 2e-3 + 0.000085 * (fck - 50)^(0.50)
        epscu = 0.0026 + 0.0144 * ((100 - fck) / 100)^4
        n = 1.4 + 9.6 * ((100 - fck) / 100)^4
    end

    if εₖ >= 0
        # Tension
        if εₖ <= eps1
            σ = E * εₖ
            ∂σ∂ε = E
        else
            if εₖ >= eps1 * (1 - E / K)
                σ = 0
                ∂σ∂ε = 0
            else
                σ = fctmfl + K * (εₖ - eps1)
                ∂σ∂ε = K
            end
        end
        #Compression
    else
        epsc = abs(εₖ)

        if epsc <= epsc0
            σ = -fcd * (1 - (1 - epsc / epsc0) .^ n) * 1000
            ∂σ∂ε = fcd * n * (1 - epsc / epsc0) .^ (n - 1) / epsc0 * 1000
        else
            σ = -fcd * 1000
            ∂σ∂ε = 0
        end

        #σ = E * εₖ
        #∂σ∂ε = E
    end

    return σ, ∂σ∂ε

end
