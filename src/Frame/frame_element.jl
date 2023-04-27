# 
# Frame element functions
#

## 
## Computes rotation matrix and element length
##

function element_geometry(cord1, cord2, ndofs)

    l = sqrt(sum((cord2 - cord1) .^ 2))

    Δx = cord2[1] - cord1[1]
    Δy = cord2[2] - cord1[2]

    ex = [Δx Δy] / l
    ey = [-Δy Δx] / l

    R = zeros(2 * ndofs, 2 * ndofs)
    aux = [ex' ey']

    R[[1, 4], [1, 4]] = aux
    R[[2, 3], [2, 3]] = aux
    R[[5, 6], [5, 6]] = aux

    return R, l
end

function frame_curvature(nelems, Mesh, len, matUk, xrel)

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
        Bₑ = intern_function(xrel[j], l, 2) * rotXYXZ
        for i in 1:len
            UₖₑL = R[dofsbe, dofsbe]' * matUk[i][elemdofs[dofsbe]]
            κHistElem[j, i] = (Bₑ*UₖₑL)[1]
        end
    end
    return κHistElem
end

function intern_function(x, l, deriv)
    if deriv == 0
        N1 = (2 * x .^ 3 - 3 * l * x .^ 2 + l .^ 3) / l^3
        N2 = (x .^ 3 - 2 * l * x .^ 2 + l .^ 2 * x) / l^2
        N3 = -(2 * x .^ 3 - 3 * l * x .^ 2) / l^3
        N4 = (x .^ 3 - l * x .^ 2) / l^2
    elseif deriv == 1
        N1 = (6 * x .^ 2 - 6 * x * l) / l^3
        N2 = (3 * x .^ 2 - 4 * l * x + l .^ 2) / l^2
        N3 = -(6 * x .^ 2 - 6 * x * l) / l^3
        N4 = (3 * x .^ 2 - 2 * l * x) / l^2
    elseif deriv == 2
        N1 = (12x - 6l) / l^3
        N2 = (6x - 4l) / l^2
        N3 = -(12x - 6l) / l^3
        N4 = (6x - 2l) / l^2
    end
    f = [N1 N2 N3 N4]
    return f
end

function intern_function_a(x, l)
    N1 = -1 / l
    N2 = 1 / l
    f = [N1 N2]
    return f
end