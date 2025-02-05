"""
    BinInfo(crystal, directions, bounds)

Construct a `BinInfo` using a Sunny `Crystal`, the local reference frame
(expressed in RLU), and the integration bounds along each axis of this reference
frame. The latter two pieces of information are what are given to Shiver when
produces a 1-D slice.
"""
struct BinInfo
    crystal    :: Sunny.Crystal
    directions :: Array{Float64, 2}
    bounds     :: Vector{Tuple{Float64, Float64}}
end


function corners_of_parallelepiped(directions, bounds; offset=Sunny.Vec3(0, 0, 0))
    points = []
    b1, b2, b3 = bounds
    for k in 1:2, j in 1:2, i in 1:2 
        q_corner = offset + directions * Sunny.Vec3(b1[i], b2[j], b3[k])
        push!(points, q_corner)
    end
    return points
end
corners_of_parallelepiped(bi::BinInfo; offset=Sunny.Vec3(0, 0, 0)) = corners_of_parallelepiped(bi.directions, bi.bounds; offset)


"""
    sample_bin_uniform_abs(q, bininfo; linear_distance=1.0)

Sample points inside the bin uniformly in absolute units. `linear_distance` is
given in Å⁻¹.
"""
function sample_bin_uniform_abs(q, bininfo::BinInfo, linear_distance=1.0)
    make_odd_up(num::Int64) = num - num%2 + 1
    (; crystal, directions, bounds) = bininfo
    (; recipvecs) = crystal
    b1, b2, b3 = bounds

    # Make oversized set of points to fill integration region.
    corners_abs = map(point -> recipvecs*point, corners_of_parallelepiped(directions, bounds))
    extremas = [extrema([corner[i] for corner in corners_abs]) for i in 1:3]
    na, nb, nc = [floor(Int64, make_odd_up(round(Int64, abs(hi - lo)/linear_distance)) / 2) for (hi, lo) in extremas]
    grid_abs = [linear_distance * Sunny.Vec3(a, b, c) for a in -na:na, b in -nb:nb, c in -nc:nc] 

    # Convert into coordinates with respect to direction vectors and filter ou# 
# t
    # points outside of integration region.
    grid_basis = [inv(directions)*inv(recipvecs)*q for q in grid_abs]
    good_points_basis = filter(grid_basis) do q
        x, y, z = q
        if b1[1] <= x <= b1[2] && b2[1] <= y <= b2[2] && b3[1] <= z <= b3[2]
            return true
        end
        false
    end

    # Add back bin center and convert to RLU before returning.
    return [Sunny.Vec3(q) + directions*point for point in good_points_basis]
end


"""
    sample_bin_uniform_rlu(q, bininfo; linear_distance=1.0)

Sample points inside the bin uniformly in RLU. `linear_distance` is
given in RLU.
"""
function sample_bin_uniform_rlu(q, bininfo::BinInfo; linear_distance=1.0)
    make_odd_down(num::Int64) = num + num%2 - 1
    (; crystal, directions, bounds) = bininfo
    (; recipvecs) = crystal
    q = Sunny.Vec3(q)
    n1, n2, n3 = map(zip(bounds, eachcol(directions))) do (bounds, direction)
        distance = norm(recipvecs * (bounds[2]*direction - bounds[1]*direction))
        floor(Int64, make_odd_down(round(Int64, distance / linear_distance)) / 2)
    end
    unit_directions = hcat([inv(recipvecs) * linear_distance * normalize(recipvecs*direction) for direction in eachcol(directions)]...)
    return [q + unit_directions * Sunny.Vec3(a, b, c) for a in -n1:n1, b in -n2:n2, c in -n3:n3][:]
end


"""
    sample_bin_mc(q, bininfo; nsamples=10)

Sample points in bin at random uniformly along each axis. The number of samples
is set with `nsamples` keyword. 
"""
function sample_bin_mc(q, bininfo::BinInfo; nsamples=10)
    (; directions, bounds) = bininfo
    sample_coef((min, max)) = min + (max-min)*rand()
    return [Sunny.Vec3(q) + directions * [sample_coef(bounds) for bounds in bounds] for _ in 1:nsamples] 
end


"""
    uniform_grid_from_bin(q, bininfo::BinInfo, spacing)


"""
function uniform_grid_from_bin(q, crystal, directions, bounds, spacing)
    make_odd_up(num::Int64) = num - num%2 + 1
    (; recipvecs) = crystal

    # Make oversized set of points to fill integration region.
    corners_abs = map(point -> recipvecs*point, corners_of_parallelepiped(directions, bounds))
    extremas = [extrema([corner[i] for corner in corners_abs]) for i in 1:3]
    na, nb, nc = [floor(Int64, make_odd_up(round(Int64, abs(hi - lo)/spacing)) / 2) for (hi, lo) in extremas]
    grid_abs = [spacing * Sunny.Vec3(a, b, c) for a in -na:na, b in -nb:nb, c in -nc:nc] 

    # Convert to RLU and add back center 
    return [q + inv(recipvecs)*point for point in grid_abs]
end

"""
    find_points_in_bin(bincenter, bininfo, points)

Given a list of `points`, return the indices of points that lie within the
specified bin.
"""
function find_points_in_bin(bincenter, bininfo, points)
    make_odd_up(num::Int64) = num - num%2 + 1
    (; directions, bounds) = bininfo
    b1, b2, b3 = bounds
    to_local_frame = inv(directions)
    bincenter_abs = to_local_frame*bincenter 
    points_grid_basis = [to_local_frame*q for q in points] 
    return findall(points_grid_basis) do q
        x, y, z = q
        x_c, y_c, z_c = bincenter_abs
        if b1[1] <= x - x_c <= b1[2] && b2[1] <= y - y_c <= b2[2] && b3[1] <= z - z_c <= b3[2]
            return true
        end
        false
    end
end
