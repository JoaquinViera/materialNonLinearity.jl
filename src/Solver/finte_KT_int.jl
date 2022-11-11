# 
# Computes Internal force and tanget stiffness matrix
#

function finte_KT_int(ElemMaterialModel, l, secParams, Uke, intBool)

    rotXYXZ = Diagonal(ones(4, 4))
    rotXYXZ[2, 2] = -1
    rotXYXZ[4, 4] = -1

    KTe = zeros(4, 4)
    Finte = zeros(4)

    # Section
    b = secParams[1]
    h = secParams[2]

    # Gauss points
    ne = 10
    ns = 10
    xge, we = gausslegendre(ne)
    xgs, ws = gausslegendre(ns)

    pgeVec = l / 2 * xge .+ l / 2
    pgsVec = h / 2 * xgs

    for j in 1:length(we)

        secFinte = 0
        secKTe = 0

        pge = pgeVec[j]

        # Bending inter functions second derivative
        B = intern_function(pge, l) * rotXYXZ

        # Strain array
        εₖVec = -pgsVec * B * Uke

        for m in 1:length(ws)
            pgs = pgsVec[m]
            εₖ = εₖVec[m]

            σ, ∂σ∂ε = constitutive_model(ElemMaterialModel, εₖ)
            secFinte = h / 2 * (b * (-B') * pgs * σ * ws[m]) .+ secFinte

            if intBool == 1
                secKTe = l / 2 * (b * ∂σ∂ε * pgs^2 * ws[m]) + secKTe
            end

        end # endfor ws

        # Tangent stiffness matrix
        if intBool == 1
            KTe = h / 2 * (B' * secKTe * B * we[j]) + KTe
        end

        # Internal force
        Finte = l / 2 * we[j] * secFinte + Finte

    end

    return Finte, KTe
end

