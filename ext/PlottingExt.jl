module PlottingExt

using SunnyAnalysisTools, GLMakie, FFTW


# Draw a bounding box for the a set of parameters that would be given to the bin sampling functions.
# If recipvecs keyword given a value, then draws boundary in absolute units (not RLU).
function SunnyAnalysisTools.draw_boundary!(ax, q, directions, bounds; color=:black, linewidth=2.0, recipvecs=nothing, kwargs...)
    corners = GLMakie.Point3f.(SunnyAnalysisTools.corners_of_parallelepiped(directions, bounds; offset=q))
    corners = !isnothing(recipvecs) ? map(q -> recipvecs*q, corners) : corners
    boundaries = [
        [corners[1], corners[2], corners[4], corners[3], corners[1]],
        [corners[5], corners[6], corners[8], corners[7], corners[5]],
        [corners[7], corners[3]], [corners[8], corners[4]],
        [corners[5], corners[1]], [corners[6], corners[2]],
    ]
    for boundary in boundaries
        lines!(ax, boundary; color, linewidth, kwargs...)
    end
end


function SunnyAnalysisTools.visualize_binning(binning::SunnyAnalysisTools.UniformBinning; labframe=false)
    fig = Figure()
    (; crystal, directions, Δs, qcenters) = binning
    bounds = [(-Δ/2, Δ/2) for Δ in Δs[1:3]]

    ax1 = LScene(fig[1,1])

    for q in qcenters
        draw_boundary!(ax1, q, directions, bounds; recipvecs = labframe ? crystal.recipvecs : nothing)
    end

    return fig
end


function SunnyAnalysisTools.plot_binned_data!(panel, data, binning::SunnyAnalysisTools.AbstractBinning; title="", plotopts=NamedTuple(), axisopts=NamedTuple())
    (; Us, Vs, Ws, Es, labels) = binning

    # Determine which dimensions are relevant
    nenergies = size(data, 1)
    liveqidcs = findall(n -> n > 1, size(data)[2:end]) .+ 1 # Indices of spatial axes with more than one point
    deadidcs = filter(idx -> !(idx in liveqidcs), [2, 3, 4])
    if nenergies == 1
        push!(deadidcs, 1)
    end
    deadidcs = sort(deadidcs)
    nqdims = length(liveqidcs)

    if nqdims == 0
        # Single-Q energy-intensity plot
        ax = Axis(panel; xlabel = "Energy transfer", ylabel="Intensity", title, axisopts...)
        lines!(ax, binning.Es, data[:,1,1,1])

    elseif nqdims == 1
        if nenergies == 1
            # Constant-energy Q-intensity plot 
            error("This plotting functionality not yet implemented.")
        else
            # 2D path plot 
            xlabel, ylabel = labels[liveqidcs[1]], "ΔE"
            ax = Axis(panel; xlabel, ylabel, title, axisopts...)

            # Downselect the live indices
            plottingdata = data 
            for (i, deadidx) in enumerate(sort(deadidcs))
                plottingdata = selectdim(plottingdata, deadidx - i + 1, 1)
            end

            # Get the values along each axis
            axisvals = [Es, Us, Vs, Ws]
            xs = axisvals[liveqidcs[1]]
            ys = axisvals[1]

            heatmap!(ax, xs, ys, plottingdata'; plotopts...)

        end
    elseif nqdims == 2
        @assert nenergies == 1 "Three dimensional plotting is not yet supported."

        xlabel, ylabel = labels[liveqidcs[1]], labels[liveqidcs[2]]
        ax = Axis(panel; xlabel, ylabel, title, axisopts...)

        # Downselect the live indices
        plottingdata = data 
        for (i, deadidx) in enumerate(sort(deadidcs))
            plottingdata = selectdim(plottingdata, deadidx - i + 1, 1)
        end

        # Get the values along each axis
        axisvals = [Es, Us, Vs, Ws]
        xs = axisvals[liveqidcs[1]]
        ys = axisvals[liveqidcs[2]]

        heatmap!(ax, xs, ys, plottingdata; plotopts...)

    elseif nqdims == 3
        @assert nenergies != 1 "Four dimensional plotting is unsupportable."
        error("Three dimensional plotting is not yet supported.")
    end
end


function SunnyAnalysisTools.plot_binned_data!(panel, observation::SunnyAnalysisTools.Observation)
    (; ints, binning) = observation
    plot_binned_data!(panel, ints, binning)
end


function SunnyAnalysisTools.plot_binned_data!(panel, calc::SunnyAnalysisTools.ModelCalculation)
    (; data, binning) = calc
    plot_binned_data!(panel, data, binning)
end


end