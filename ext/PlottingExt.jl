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

function SunnyAnalysisTools.visualize_convolution(q_target, bi, spec)
    fig = Figure(size=(1800, 400))
    (; crystal) = bi

    ax1 = LScene(fig[1,1])
    ax2 = LScene(fig[1,2])
    ax3 = LScene(fig[1,3])
    ax4 = LScene(fig[1,4])

    pts_abs = map(q -> crystal.recipvecs*(q+q_target), spec.points)
    xs = [p[1] for p in pts_abs[:,1,1]]
    ys = [p[2] for p in pts_abs[1,:,1]]
    zs = [p[3] for p in pts_abs[1,1,:]]

    scatter!(ax1, pts_abs[:])
    draw_boundary!(ax1, q_target, bi.directions, bi.bounds; recipvecs=crystal.recipvecs)

    volume!(ax2, xs[1]..xs[end], ys[1]..ys[end], zs[1]..zs[end], real.(ifft(spec.kernel)))
    draw_boundary!(ax2, q_target, bi.directions, bi.bounds; recipvecs=crystal.recipvecs)

    scatter!(ax3, map(p -> p + q_target, spec.points[:]))
    draw_boundary!(ax3, q_target, bi.directions, bi.bounds)

    scatter!(ax4, map(p -> p + q_target, spec.points[spec.interior_idcs]))
    draw_boundary!(ax4, q_target, bi.directions, bi.bounds)

    return fig
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


end