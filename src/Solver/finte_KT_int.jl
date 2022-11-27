# 
# Computes Internal force and tanget stiffness matrix
#

function finte_KT_int(ElemMaterialModel, l, secParams, Uke, intBool, σArr, time, elem)

    rotXYXZ = Diagonal(ones(4, 4))
    rotXYXZ[2, 2] = -1
    rotXYXZ[4, 4] = -1

    # Axial
    KTea = zeros(2, 2)
    Fintea = zeros(2)


    # Bending
    KTeb = zeros(4, 4)
    Finteb = zeros(4)

    # Diagonal
    #KTeab = zeros(2, 4)

    # Section
    b = secParams[1]
    h = secParams[2]

    # Gauss points

    xge, we = gausslegendre(ElemMaterialModel.ne)
    xgs, ws = gausslegendre(ElemMaterialModel.ns)
    minXe = minimum(xge)
    pgeVec = l / 2 * xge .+ l / 2
    pgsVec = h / 2 * xgs

    # Disp vector
    Ua = view(Uke, [1, 4])
    Ub = view(Uke, [2, 3, 5, 6])

    for j in 1:length(we)

        secFintea = 0
        secKTea = 0

        secFinteb = 0
        secKTeb = 0

        #secKTeab = 0

        pge = view(pgeVec, j)[1]

        # Axial interpolation functions first derivative
        Ba = intern_function_a(pge, l)
        # Bending interpolation functions second derivative
        Bb = intern_function(pge, l) * rotXYXZ

        # Strain array
        εₖVec = -pgsVec * Bb * Ub .+ Ba * Ua

        for m in 1:length(ws)
            pgs = view(pgsVec, m)[1]
            εₖ = view(εₖVec, m)[1]

            σ, ∂σ∂ε = constitutive_model(ElemMaterialModel, εₖ)

            secFintea = h / 2 * (b * Ba' * σ * ws[m]) .+ secFintea
            secFinteb = h / 2 * (b * (-Bb') * pgs * σ * ws[m]) .+ secFinteb

            if intBool == 1
                secKTea = h / 2 * (b * ∂σ∂ε * ws[m]) + secKTea
                secKTeb = h / 2 * (b * ∂σ∂ε * pgs^2 * ws[m]) + secKTeb
                #secKTeab = h / 2 * (b * ∂σ∂ε * (-pgs) * ws[m]) + secKTeab
            else
                if minXe == xge[j] && elem == 1
                    σArr[time+1][m] = σ
                end
            end

        end # endfor ws

        # Tangent stiffness matrix
        if intBool == 1

            KTea = l / 2 * (Ba' * secKTea * Ba * we[j]) + KTea
            KTeb = l / 2 * (Bb' * secKTeb * Bb * we[j]) + KTeb

            #KTeab = l / 2 * (Ba' * secKTeab * Bb * we[j]) + KTeab
        end

        # Internal force
        Fintea = l / 2 * we[j] * secFintea + Fintea
        Finteb = l / 2 * we[j] * secFinteb + Finteb

    end

    return Finteb, Fintea, σArr, KTeb, KTea
end

