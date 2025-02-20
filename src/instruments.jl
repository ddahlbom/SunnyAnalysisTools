abstract type AbstractInstrumentSpec end

abstract type ChopperSpec <: AbstractInstrumentSpec end
abstract type TripleAxisSpec <: AbstractInstrumentSpec end

struct CNCSSpec <: ChopperSpec
    Ei   :: Float64
    L1   :: Float64
    L2   :: Float64
    L3   :: Float64
    Δtp  :: Float64
    Δtc  :: Float64
    Δtd  :: Float64
    Δθ   :: Float64

    function CNCSSpec(Ei; Δθ = 1.5)
        Instruments = pyimport("pychop.Instruments")
        CNCS = Instruments.Instrument("CNCS") 
        Δθ = Δθ * (π/180)

        distances = CNCS.chopper_system.getDistances()
        x0, _, x1, x2, xm = map(val -> pyconvert(Float64, val), distances)  # xa not used
        moderator_width  = CNCS.moderator.getWidthSquared(Ei) |> x -> pyconvert(Float64, x) |> sqrt
        chopper_width = CNCS.chopper_system.getWidthSquared(Ei)[0] |> x -> pyconvert(Float64, x) |> sqrt

        L1, L2, L3 = x0 - xm, x1, x2
        Δtp = moderator_width * (1 - (xm/x0))
        Δtc = chopper_width 
        Δtd = 0.0

        new(Ei, L1, L2, L3, Δtp, Δtc, Δtd, Δθ)
    end
end

