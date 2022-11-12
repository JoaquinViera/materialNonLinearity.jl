
function assembler(Section, MaterialModel, Mesh, Uₖ, intBool)

    ndofs = 2 # degrees of freedom per node
    nnodes = size(Mesh.nodesMat, 1)
    nelems = size(Mesh.conecMat, 1)

    KTₖ = zeros(ndofs * nnodes, ndofs * nnodes)
    FintₖL = zeros(nnodes * ndofs)
    Fintₖ = zeros(nnodes * ndofs, 1)

    nodes = view(Mesh.conecMat, (nelems*2+1):(nelems*3))

    for i in 1:nelems
        # Elem nodes, dofs, Material and geometry

        nodeselem = nodes[i]
        elemdofs = nodes2dofs(nodeselem, ndofs)

        R, l = element_geometry(Mesh.nodesMat[nodeselem[1], :], Mesh.nodesMat[nodeselem[2], :], ndofs)
        ElemSecParams = [Section.b Section.h]
        ElemMaterial = MaterialModel

        # Elem disps in local system
        UₖeL = R' * Uₖ[elemdofs]

        # Internal force & Tangent mat
        #intBool = 1
        #Finte, KTe = finte_KT_int(ElemMaterial, l, ElemSecParams, UₖeL, intBool)
        if intBool == 1
            Finte, KTe = finte_KT_int(ElemMaterial, l, ElemSecParams, UₖeL, intBool)
            # Assemble tangent stiffness matrix
            KTₖ[elemdofs, elemdofs] = KTₖ[elemdofs, elemdofs] + R * KTe * R'
        else
            Finte, ~ = finte_KT_int(ElemMaterial, l, ElemSecParams, UₖeL, intBool)
        end

        # Assemble internal force
        FintₖL[elemdofs] = Finte
        Fintₖ[elemdofs] = Fintₖ[elemdofs] + R * Finte

    end

    if intBool == 1
        return Fintₖ, KTₖ
    else
        return Fintₖ
    end
end