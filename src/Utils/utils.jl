


## =======================================
## Computes nodal degrees of freedom 
## =======================================

#=
function nodes2dofs(nodes, degreespernode)
    n = length(nodes)
    dofs = Vector{Int64}(undef, n*degreespernode)
    for i = 1:n
        for j = 1:degreespernode
            @inbounds dofs[(i-1)*degreespernode + j] = degreespernode*(nodes[i]-1) + j
        end
    end
    return dofs
end
=#

function nodes2dofs(nodes, ndofs)

    n = length(nodes)
    gdl = zeros(Int64, n * ndofs)
    x = Vector(1:ndofs)
    for i in 1:n
        gdl[((i-1)*ndofs).+x] = (nodes[i] - 1) * ndofs .+ x
    end
    return gdl
end
