module PythonToolsExt

using PythonCall, SunnyAnalysisTools


function SunnyAnalysisTools.cncs(; Ei, Δθ = 1.5)
    Instruments = pyimport("PyChop.Instruments")
    CNCS = Instruments.Instrument("CNCS")

    Δθ = Δθ * (π/180) # Convert to radians
    distances = CNCS.chopper_system.getDistances()
    x0, _, x1, x2, xm = map(val -> pyconvert(Float64, val), distances)  # xa not used
    moderator_width  = CNCS.moderator.getWidthSquared(Ei) |> x -> pyconvert(Float64, x) |> sqrt
    chopper_width = CNCS.chopper_system.getWidthSquared(Ei)[0] |> x -> pyconvert(Float64, x) |> sqrt

    L1, L2, L3 = x0 - xm, x1, x2
    Δtp = moderator_width * (1 - (xm/x0))
    Δtc = chopper_width 
    Δtd = 0.0

    return SunnyAnalysisTools.ChopperSpec("CNCS", Ei, L1, L2, L3, Δtp, Δtc, Δtd, Δθ)
end

function SunnyAnalysisTools.spins()
    params = Dict{String, Any}()
    return SunnyAnalysisTools.TripleAxisSpec("SPINS", params)
end

function SunnyAnalysisTools.hmi(; hires=false)
    params = Dict{String, Any}()
    if hires
        params["resolution"] = "high"
    else
        params["resolution"] = "low"
    end
    return SunnyAnalysisTools.TripleAxisSpec("HMI", params)
end

end