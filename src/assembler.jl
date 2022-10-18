
function assembler(Section, MaterialModel, Mesh, Uk, intBool)

    ndofs = 2 # degrees of freedom per node
    nnodes = size(Mesh.nodesMat, 1)
    nelems = size(Mesh.conecMat, 1)

    KTk = zeros(ndofs * nnodes, ndofs * nnodes)
    FintkL = zeros(nnodes * ndofs)
    Fintk = zeros(nnodes * ndofs, 1)

    for i in 1:nelems
        # Elem nodes, dofs, Material and geometry
        nodeselem = Mesh.conecMat[i, 3]
        elemdofs = nodes2dofs(nodeselem[:], ndofs)
        R, l = element_geometry(Mesh.nodesMat[nodeselem[1], :], Mesh.nodesMat[nodeselem[2], :], ndofs)
        elemSecParams = Section.params
        elemMaterial = MaterialModel

        # Elem disps in local system
        UkeL = R' * Uk[elemdofs]

        # Internal force & Tangent mat
        #intBool = 1
        #Finte, KTe = finte_KT_int(elemMaterial, l, elemSecParams, UkeL, intBool)
        if intBool == 1
            Finte, KTe = finte_KT_int(elemMaterial, l, elemSecParams, UkeL, intBool)
            # Assemble tangent stiffness matrix
            KTk[elemdofs, elemdofs] = KTk[elemdofs, elemdofs] + R * KTe * R'
        else
            Finte, ~ = finte_KT_int(elemMaterial, l, elemSecParams, UkeL, intBool)
        end

        # Assemble internal force
        FintkL[elemdofs] = Finte
        Fintk[elemdofs] = Fintk[elemdofs] + R * Finte

    end

    if intBool == 1
        return Fintk, KTk
    else
        return Fintk
    end
end