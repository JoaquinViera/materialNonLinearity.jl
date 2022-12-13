
# Bending moment plot

function BendingMomentPlot(timesPlot, Mesh, PlotSettings, matM)

    ndofs = 3
    dofsM = [3, 6]
    nelems = size(Mesh.conecMat, 1)

    figsMat = Array{Plots.Plot{Plots.GRBackend},1}()

    label = "M(x)"
    for i in timesPlot
        title = "Bending Moment, step $i"
        fig = plot(minorgrid=PlotSettings.minorGridBool, legend=PlotSettings.legendPosition, title=title)

        for j = 1:nelems
            nodeselem = Mesh.conecMat[j, ndofs]
            xₑ = [Mesh.nodesMat[nodeselem[1], 1], Mesh.nodesMat[nodeselem[2], 1]]
            Mₑ = abs.(matM[j][i][dofsM])

            if j != nelems
                plot!(fig, xₑ, Mₑ, color=PlotSettings.color, lw=PlotSettings.lw, ms=PlotSettings.ms, label="")
            else
                plot!(fig, xₑ, Mₑ, color=PlotSettings.color, lw=PlotSettings.lw, ms=PlotSettings.ms, label=label)
            end
            xlabel!("x")
            ylabel!("M")
        end
        push!(figsMat, fig)
    end
    return figsMat
end

# Deformed shape ----- TO DO

# function DeformedShape(ndivs, nelems, matUk, timesPlot, Mesh, Section, MaterialModel, PlotSettings, matM)

#     ndofs = 3
#     rotXYXZ = Diagonal(ones(4, 4))
#     rotXYXZ[2, 2] = -1
#     rotXYXZ[4, 4] = -1
#     dofsbe = [2, 3, 5, 6]
#     dofsM = [3, 6]

#     title = "Deformed Shape"
#     label = "v(x)"
#     fig = plot(minorgrid=PlotSettings.minorGridBool, legend=PlotSettings.legendPosition, title=title)
#     for i in timesPlot
#         for j = 1:nelems
#             nodeselem = Mesh.conecMat[j, ndofs]
#             elemdofs = nodes2dofs(nodeselem[:], ndofs)
#             R, l = element_geometry(view(Mesh.nodesMat, nodeselem[1], :), view(Mesh.nodesMat, nodeselem[2], :), ndofs)
#             UₖₑL = R[dofsbe, dofsbe]' * matUk[i][elemdofs[dofsbe]]

# xₑ = collect(Mesh.nodesMat[nodeselem[1], 1]:l/ndivs:Mesh.nodesMat[nodeselem[2], 1])
# xplot = collect(0:l/ndivs:l)
# Mₑ = zeros(length(xₑ))
# Mₑ = abs.(matM[j][i][dofsM])

# for m in 1:length(xₑ)
#     Bₑ = intern_function(xplot[m], l) * rotXYXZ
#     κₑ = (Bₑ*UₖₑL)[1]
#     Mₑ[m] = -κₑ * MaterialModel.E * Section.Iy
# end

# if j != nelems
# plot!(fig, xₑ, Mₑ, color=PlotSettings.color, lw=PlotSettings.lw, ms=PlotSettings.ms, label="")
# else
# plot!(fig, xₑ, Mₑ, color=PlotSettings.color, lw=PlotSettings.lw, ms=PlotSettings.ms, label=label)
#             end
#             xlabel!("x")
#             ylabel!("v")
#         end
#     end
#     return fig
# end