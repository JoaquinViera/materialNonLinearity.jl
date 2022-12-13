
function assembler(Section, MaterialModel, Mesh, Uₖ, intBool, σArr, time, StressArraySets, F)

    ndofs = 3 # degrees of freedom per node
    nnodes = size(Mesh.nodesMat, 1)
    nelems = size(Mesh.conecMat, 1)

    KTₖ = zeros(ndofs * nnodes, ndofs * nnodes)
    Fintₖ = zeros(nnodes * ndofs, 1)

    nodes = view(Mesh.conecMat, (nelems*2+1):(nelems*3))

    dofsa = [1, 4]
    dofsb = [2, 3, 5, 6]


    for i in 1:nelems
        # Elem nodes, dofs, Material and geometry

        nodeselem = nodes[i]
        elemdofs = nodes2dofs(nodeselem, ndofs)
        dofsbe = elemdofs[dofsb]
        dofsae = elemdofs[dofsa]

        R, l = element_geometry(view(Mesh.nodesMat, nodeselem[1], :), view(Mesh.nodesMat, nodeselem[2], :), ndofs)
        ElemSecParams = [Section.b Section.h]
        ElemMaterial = MaterialModel

        # Elem disps in local system
        UₖeL = R' * view(Uₖ, elemdofs)

        # Internal force & Tangent mat
        if intBool == 1
            Finteb, Fintea, σArr, KTeb, KTea = finte_KT_int(ElemMaterial, l, ElemSecParams, UₖeL, intBool, σArr, time, i, StressArraySets)
            # Assemble tangent stiffness matrix
            #KTₖ[elemdofs, elemdofs] = KTₖ[elemdofs, elemdofs] + R * KTe * R'
            #KTₖ[dofsb, dofsb] = KTₖ[dofsb, dofsb] + R[dofsb, dofsb] * KTeb * R[dofsb, dofsb]'

            KTₖ[dofsbe, dofsbe] = KTₖ[dofsbe, dofsbe] + R[dofsb, dofsb] * KTeb * R[dofsb, dofsb]'
            KTₖ[dofsae, dofsae] = KTₖ[dofsae, dofsae] + R[dofsa, dofsa] * KTea * R[dofsa, dofsa]'

        else
            Finteb, Fintea, σArr = finte_KT_int(ElemMaterial, l, ElemSecParams, UₖeL, intBool, σArr, time, i, StressArraySets)
        end

        # Assemble internal force
        Fintₖ[dofsbe] = Fintₖ[dofsbe] + R[dofsb, dofsb] * Finteb
        Fintₖ[dofsae] = Fintₖ[dofsae] + R[dofsa, dofsa] * Fintea
        F[i][time+1][dofsb] = Finteb
        F[i][time+1][dofsa] = Fintea
    end

    if intBool == 1
        return Fintₖ, σArr, F, KTₖ
    else
        return Fintₖ, σArr, F
    end
end