################################################################################
# Types
################################################################################

abstract type AbstractBinning end

mutable struct UniformBinning <: AbstractBinning
    const crystal     :: Sunny.Crystal

    const directions  :: Sunny.Mat3        # Definition of scattering coordinates (relative to RLU)
    labels            :: Vector{String} 

    const Us          :: Vector{Float64}
    const Vs          :: Vector{Float64}
    const Ws          :: Vector{Float64}
    const Es          :: Vector{Float64}

    const qcenters    :: Array{Vector{Float64}, 3}
    const qbase       :: Sunny.Vec3 

    const Δs          :: Vector{Float64}


    function UniformBinning(crystal, directions, Us, Vs, Ws, Es; labels=["U", "V", "W"])
        # Ensure that spacing is uniform 
        Δs = map([Us, Vs, Ws, Es]) do vals
            Δs = vals[2:end] .- vals[1:end-1]
            @assert all(Δ -> Δ ≈ Δs[1], Δs) "Step sizes must all be equal for a UniformBinning"
            Δs[1]
        end

        # If only bounds are given (as opposed to a list) determine center point.
        Us, Vs, Ws, Es = map([Us, Vs, Ws, Es]) do vals
            if length(vals) == 2
                [(vals[2]+vals[1])/2]
            else
                vals
            end
        end

        # Assemble bincenters in RLU
        qcenters = [directions*[U, V, W] for U in Us, V in Vs, W in Ws]
        qbase = qcenters[1,1,1] .- 0.5*(directions*Δs[1:3])

        new(crystal, directions, labels, Us, Vs, Ws, Es, qcenters, qbase, Δs)
    end
end 

function Base.show(io::IO, binning::UniformBinning)
    (; qcenters, Δs, crystal, directions) = binning
    println(io, "UniformBinning")
    nH, nK, nL = size(qcenters)

    println(io, crystal)
    println(io, "Directions: ", directions) # Make this nicer
    print(io, "H: ")

    # Make this nicer with select dim or similar
    if nH == 1
        println(io, "$(qcenters[1,1,1][1]), ΔH=$(round(Δs[1], digits=3))")
    else
        println(io, "$(qcenters[1,1,1][1])...$(qcenters[end,1,1][1]), ΔH=$(Δs[1])")
    end
    print(io, "K: ")
    if nK == 1
        println(io, "$(qcenters[1,1,1][2]), ΔK=$(Δs[2])")
    else
        println(io, "$(qcenters[1,1,1][2])...$(qcenters[1,end,1][2]), ΔH=$(Δs[2])")
    end
    print(io, "L: ")
    if nL == 1
        println(io, "$(qcenters[1,1,1][3]), ΔL=$(round(Δs[3], digits=3))")
    else
        println(io, "$(qcenters[1,1,1][3])...$(qcenters[end,1,1][3]), ΔH=$(round(Δs[3], digits=3))")
    end

end

################################################################################
# Methods 
################################################################################

function sample_binning(binning::UniformBinning; nperbin=1, nghosts=0, nperebin=1)
    (; crystal, directions, Δs, qcenters, qbase, Es) = binning
    (; recipvecs) = crystal
    nperbins = isa(nperbin, Number) ? nperbin * ones(3) : nperbin
    nghosts = isa(nghosts, Number) ? nghosts * ones(3) : nghosts

    # Sample q points
    directions_abs = recipvecs*directions
    increments = directions_abs * diagm(Δs[1:3] ./ nperbins)
    tops = [N+nghosts for (N, nghosts) in zip([size(qcenters)...], nghosts)] .* nperbins .- 1
    bottoms = [-nghosts for nghosts in nghosts] .* nperbins
    bounds = [b:t for (b, t) in zip(bottoms,tops)]
    offset = inv(recipvecs)*increments*[0.5, 0.5, 0.5] 
    qpoints = [qbase + offset + inv(recipvecs)*increments*[Na, Nb, Nc] for Na in bounds[1], Nb in bounds[2], Nc in bounds[3]]

    # Sample energies
    ΔE = Δs[4]
    increment = Δs[4] / nperebin
    ebase = Es[1] - ΔE/2
    offset = increment*0.5
    epoints = [ebase + offset + increment*N for N in 0:(length(Es)*nperebin-1)]

    return (; qpoints, epoints)
end


function find_points_in_bin(bincenter, directions, bounds, points)
    b1, b2, b3 = bounds

    # Convert all information into local bin coordinates. 
    to_local_frame = inv(directions)

    bincenter = to_local_frame*bincenter 
    points = [to_local_frame*q for q in points] 

    return findall(points) do q
        x, y, z = q
        x_c, y_c, z_c = bincenter
        if b1[1] <= x - x_c <= b1[2] && b2[1] <= y - y_c <= b2[2] && b3[1] <= z - z_c <= b3[2]
            return true
        end
        false
    end
end


function corners_of_parallelepiped(directions, bounds; offset=[0., 0, 0])
    points = []
    b1, b2, b3 = bounds
    for k in 1:2, j in 1:2, i in 1:2 
        q_corner = offset + directions * [b1[i], b2[j], b3[k]]
        push!(points, q_corner)
    end
    return points
end