struct TripleAxisPath
    crystal         :: Sunny.Crystal
    HKLs            :: Vector{Sunny.Vec3}
    Es              :: Vector{Float64}
    projection      :: Sunny.Mat3 
end


# Assumes Q list and E list are identical length (e.g., for
# a 1D). Think more about interface and data format.
function TripleAxisPath(crystal, HKLs, Es; projection=I(3))
    @assert length(HKLs) == length(Es) "Number of wave vectors must equal number of energies"
    TripleAxisPath(crystal, HKLs, Es, projection)
end

struct TripleAxis2DContour
    crystal         :: Sunny.Crystal
    HKLs            :: Vector{Sunny.Vec3}
    Es              :: Vector{Float64}
    projection      :: Sunny.Mat3 
end

# For a product
function TripleAxis2DContour(crystal, HKLs, Es; projection=I(3))
    TripleAxis2DContour(crystal, HKLs, Es, projection)
end