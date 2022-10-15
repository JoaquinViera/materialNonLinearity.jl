# =======================================
# Computes nodal degrees of freedom 
# =======================================
function nodes2dofs(nodes, ndofs)

    n = length(nodes)
    gdl = zeros(Int64, n * ndofs)
    for i in 1:n
        gdl[(i-1)*ndofs.+Vector(1:ndofs)] = (nodes[i] - 1) * ndofs .+ Vector(1:ndofs)
    end
    return gdl
end