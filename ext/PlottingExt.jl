module PlottingExt

using SunnyHelpersORNL, GLMakie


# Draw a bounding box for the a set of parameters that would be given to the bin sampling functions.
# If recipvecs keyword given a value, then draws boundary in absolute units (not RLU).
function SunnyHelpersORNL.draw_boundary!(ax, q, directions, bounds; color=:black, linewidth=2.0, recipvecs=nothing, kwargs...)
    corners = GLMakie.Point3f.(SunnyHelpersORNL.corners_of_parallelepiped(directions, bounds; offset=q))
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

end