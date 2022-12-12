
# Bending moment plot

function BendingMomentPlot(ndivs, nelems, matUk, timesPlot, Mesh, Section, MaterialModel, PlotSettings)

    ndofs = 3
    rotXYXZ = Diagonal(ones(4, 4))
    rotXYXZ[2, 2] = -1
    rotXYXZ[4, 4] = -1
    dofsbe = [2, 3, 5, 6]

    title = "Bending Moment"
    label = "M(x)"
    fig = plot(minorgrid=PlotSettings.minorGridBool, legend=PlotSettings.legendPosition, title=title)
    for i in timesPlot
        for j = 1:nelems
            nodeselem = Mesh.conecMat[j, ndofs]
            elemdofs = nodes2dofs(nodeselem[:], ndofs)
            R, l = element_geometry(view(Mesh.nodesMat, nodeselem[1], :), view(Mesh.nodesMat, nodeselem[2], :), ndofs)
            UₖₑL = R[dofsbe, dofsbe]' * matUk[i][elemdofs[dofsbe]]

            xₑ = collect(Mesh.nodesMat[nodeselem[1], 1]:l/ndivs:Mesh.nodesMat[nodeselem[2], 1])
            xplot = collect(0:l/ndivs:l)
            Mₑ = zeros(length(xₑ))

            for m in 1:length(xₑ)
                Bₑ = intern_function(xplot[m], l) * rotXYXZ
                κₑ = (Bₑ*UₖₑL)[1]
                Mₑ[m] = -κₑ * MaterialModel.E * Section.Iy
            end
            println(Mₑ[m])
            if j != nelems
                plot!(fig, xₑ, Mₑ, color=PlotSettings.color, lw=PlotSettings.lw, ms=PlotSettings.ms, label="")
            else
                plot!(fig, xₑ, Mₑ, color=PlotSettings.color, lw=PlotSettings.lw, ms=PlotSettings.ms, label=label)
            end
            xlabel!("x")
            ylabel!("M")
        end
    end
    return fig
end