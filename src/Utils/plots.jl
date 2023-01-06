
# Bending moment plot

function BendingMomentPlot(timesPlot, Mesh, PlotSettings, matFinte)

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
            Mₑ = abs.(matFinte[j][i][dofsM])

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

# Deformed shape plot

function DeformedShapePlot(timesPlot, Mesh, PlotSettings, matUk)

    ndivs = 10
    ndofs = 3
    rotXYXZ = Diagonal(ones(4, 4))
    rotXYXZ[2, 2] = -1
    rotXYXZ[4, 4] = -1
    dofsbe = [2, 3, 5, 6]
    # dofsM = [3, 6]
    nelems = size(Mesh.conecMat, 1)

    figsMat = Array{Plots.Plot{Plots.GRBackend},1}()

    label = "v(x)"
    for i in timesPlot
        title = "Deformed Shape"
        fig = plot(minorgrid=PlotSettings.minorGridBool, legend=PlotSettings.legendPosition, title=title)
        for j = 1:nelems
            nodeselem = Mesh.conecMat[j, ndofs]
            elemdofs = nodes2dofs(nodeselem[:], ndofs)
            R, l = element_geometry(view(Mesh.nodesMat, nodeselem[1], :), view(Mesh.nodesMat, nodeselem[2], :), ndofs)
            UₖₑL = R[dofsbe, dofsbe]' * matUk[i][elemdofs[dofsbe]]

            xₑ = collect(Mesh.nodesMat[nodeselem[1], 1]:l/ndivs:Mesh.nodesMat[nodeselem[2], 1])
            xplot = collect(0:l/ndivs:l)
            vₑ = zeros(length(xₑ))

            for m in 1:length(xₑ)
                Nₑ = intern_function(xplot[m], l, 0) * rotXYXZ
                vₑ[m] = (Nₑ*UₖₑL)[1]
            end

            if j != nelems
                plot!(fig, xₑ, vₑ, color=PlotSettings.color, lw=PlotSettings.lw, ms=PlotSettings.ms, label="")
            else
                plot!(fig, xₑ, vₑ, color=PlotSettings.color, lw=PlotSettings.lw, ms=PlotSettings.ms, label=label)
            end
            xlabel!("x")
            ylabel!("v")
        end
        push!(figsMat, fig)
    end
    return figsMat
end

# Constitutive model plot
function ConstitutiveModelPlot(model, ε::Vector, divs::Integer, factor_ε::Float64, factor_σ::Float64)
    ε_min, ε_max = (maximum(ε), minimum(ε))
    ε_Vector = LinRange(ε_min, ε_max, divs)

    title = "Constitutive Model σ-ε"
    fig = plot(minorgrid=1, title=title)
    σ_Vector = zeros(length(ε_Vector))

    for i in 1:length(ε_Vector)
        σ, nothing = constitutive_model(model, ε_Vector[i])
        σ_Vector[i] = σ
    end
    plot!(fig, ε_Vector * factor_ε, σ_Vector * factor_σ, color="blue", lw=3, ms=2, legend=false)
    plot!(fig, [0.0, 0.0], [-maximum(abs.(σ_Vector)), maximum(abs.(σ_Vector))] * factor_σ, lw=1.5, ms=1, label="", color=:"black")
    plot!(fig, [minimum(ε_Vector), maximum(ε_Vector)] * factor_ε, [0.0, 0.0], lw=1.5, ms=1, label="", color=:"black")
    xlabel!("ε")
    ylabel!("σ")
    return fig
end