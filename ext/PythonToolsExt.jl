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

    SunnyAnalysisTools.ChopperSpec("CNCS", Ei, L1, L2, L3, Δtp, Δtc, Δtd, Δθ)
end

end