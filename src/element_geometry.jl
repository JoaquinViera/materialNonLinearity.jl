# 
# Computes rotation matrix and elem length
#

function element_geometry(cord1, cord2, ndofs)

    l = sqrt(sum((cord2 - cord1) .^ 2))

    Deltax = cord2[1] - cord1[1]
    Deltay = cord2[2] - cord1[2]

    ex = [Deltax Deltay] / l
    ey = [-Deltay Deltax] / l

    R = zeros(2 * ndofs, 2 * ndofs)
    aux = [ex' ey']

    R[1:2, 1:2] = aux
    R[3:4, 3:4] = aux

    return R, l
end