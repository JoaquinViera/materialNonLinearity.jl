## =======================================
## Computes nodal degrees of freedom 
## =======================================

function nodes2dofs(nodes, ndofs)
    n = length(nodes)
    gdl = zeros(Int64, n * ndofs)
    vec = Vector(1:ndofs)
    for i in 1:n
        gdl[(i-1)*ndofs.+vec] = (nodes[i] - 1) * ndofs .+ vec
    end
    return gdl
end

function frame_curvature(nelems, Mesh, len, matUk)

    ndofs = 3
    κHistElem = zeros(nelems, len)
    rotXYXZ = Diagonal(ones(4, 4))
    rotXYXZ[2, 2] = -1
    rotXYXZ[4, 4] = -1
    dofsbe = [2, 3, 5, 6]

    for j in 1:nelems
        nodeselem = Mesh.conecMat[j, ndofs]
        elemdofs = nodes2dofs(nodeselem[:], ndofs)
        R, l = element_geometry(view(Mesh.nodesMat, nodeselem[1], :), view(Mesh.nodesMat, nodeselem[2], :), ndofs)
        Bₑ = intern_function(0, l) * rotXYXZ
        for i in 1:len
            UₖₑL = R[dofsbe, dofsbe]' * matUk[i][elemdofs[dofsbe]]
            κHistElem[j, i] = (Bₑ*UₖₑL)[1]
        end
    end
    return κHistElem
end
