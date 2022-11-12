


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