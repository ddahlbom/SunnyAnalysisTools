struct TripleAxisPath
    crystal         :: Sunny.Crystal
    HKLs            :: Vector{Sunny.Vec3}
    Es              :: Vector{Float64}
    projection      :: Sunny.Mat3 
end

function TripleAxisPath(crystal, HKLs, Es; projection=I(3))
    TripleAxisPath(crystal, HKLs, Es, projection)
end