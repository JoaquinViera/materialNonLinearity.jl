# isotropic Bilinear

# using LinearAlgebra, BenchmarkTools

#using ProgressMeter
#=
p = Progress(10000, dt=0.25, barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:cyan)

i = 1
function prueba(i)

    while i <= 10000
        ones(300) * ones(300)'
        next!(p)
        i = i + 1
        #println(i)
    end

end
=#
#=

function prueba1(t)
    Uk = zeros(t)
    aux = ones(t)
    matUk = vcat(Uk, [])
    for i = 1:t
        Uk = matUk[:, i]
        Uk = Uk + aux
        matUk = hcat(matUk, Uk)
    end
end

function prueba2(t)
    Uk = zeros(t)
    aux = ones(t)
    matUk = Vector{Vector{Float64}}()
    push!(matUk, Uk)
    for i = 1:t
        Uk = view(matUk, i)
        Uk = Uk[1] + aux
        push!(matUk, Uk)
    end
end

function prueba3(t)
    Uk = zeros(t)
    aux = ones(t)
    matUk = Vector{Vector{Float64}}()
    push!(matUk, Uk)
    for i = 1:t
        Uk = view(matUk, i)[1] + aux
        push!(matUk, Uk)
    end
end

function prueba4(t)
    Uk = zeros(t)
    aux = ones(t)
    matUk = Vector{Vector{Float64}}(undef, t)
    matUk[1] = Uk
    for i = 1:(t-1)
        Uk = view(matUk, i)[1]
        matUk[i+1] = Uk
    end
end

t = 100

@btime prueba1(t)
@btime prueba2(t)
@btime prueba3(t)
@btime prueba4(t)

=#

# without mods
#Memory estimate: 756.18 MiB, allocs estimate: 17773651.
#=
1.857 s (17774411 allocations: 756.21 MiB)

1.546 s (17774316 allocations: 756.16 MiB)

1.601 s (17774221 allocations: 756.11 MiB)

1.587 s (17774601 allocations: 756.08 MiB)

1.657 s (17774404 allocations: 756.21 MiB)

1.623 s (17774309 allocations: 756.20 MiB)

2.015 s (17774356 allocations: 753.64 MiB)

1.863 s (17774356 allocations: 753.57 MiB)

1.598 s (17757253 allocations: 752.27 MiB)

1.601 s (17755955 allocations: 752.20 MiB)

1.694 s (17755681 allocations: 751.24 MiB)

1.622 s (17755776 allocations: 751.24 MiB)

1.520 s (17754676 allocations: 750.22 MiB)

1.556 s (17754676 allocations: 750.22 MiB)

  1.640 s (17754771 allocations: 750.27 MiB)

 1.743 s (17755776 allocations: 751.29 MiB)

 1.616 s (17754832 allocations: 750.27 MiB)

 3.765 s (44637435 allocations: 1.76 GiB)
 3.678 s (44637337 allocations: 1.76 GiB)
 3.575 s (44637241 allocations: 1.76 GiB)

 12.345 s (163687389 allocations: 6.34 GiB)
 14.176 s (163686866 allocations: 6.34 GiB)
 16.401 s (163686468 allocations: 6.34 GiB)
 14.414 s (163685946 allocations: 6.34 GiB)
 16.001 s (163686462 allocations: 6.34 GiB)
 15.196 s (163686063 allocations: 6.34 GiB)
 16.388 s (163685017 allocations: 6.34 GiB)
 13.463 s (163685017 allocations: 6.34 GiB)
 13.596 s (163686586 allocations: 6.34 GiB)
 14.700 s (163685017 allocations: 6.34 GiB)

 13.386 s (163690226 allocations: 6.34 GiB)
 14.605 s (163690321 allocations: 6.34 GiB)
 12.528 s (138196722 allocations: 5.82 GiB)
 12.518 s (138196599 allocations: 5.82 GiB)
 12.068 s (138196474 allocations: 5.82 GiB)
 14.046 s (138196296 allocations: 5.82 GiB)
 12.447 s (138028744 allocations: 5.81 GiB)
 12.789 s (138070584 allocations: 5.81 GiB)
 12.861 s (138028791 allocations: 5.81 GiB)
 13.007 s (138028743 allocations: 5.81 GiB)
 12.314 s (138027697 allocations: 5.81 GiB)
 =#

 #=
using BenchmarkTools

function nodes2dofs(nodes, ndofs)

    n = length(nodes)
    gdl = zeros(Int64, n * ndofs)
    vec = Vector(1:ndofs)
    for i in 1:n
        gdl[(i-1)*ndofs.+vec] = (nodes[i] - 1) * ndofs .+ vec
    end
    return gdl
end

@btime nodes2dofs([1, 2], 6)


function nodes2dofs2(nodes, degreespernode)
    n = length(nodes)
    dofs = Vector{Int64}(undef, n * degreespernode)
    for i = 1:n
        for j = 1:degreespernode
            @inbounds dofs[(i-1)*degreespernode+j] = degreespernode * (nodes[i] - 1) + j
        end
    end
    return dofs
end


@btime nodes2dofs([1, 2], 6)

=#