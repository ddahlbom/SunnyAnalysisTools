abstract type AbstractInstrumentSpec end

# abstract type TripleAxisSpec <: AbstractInstrumentSpec end

struct ChopperSpec <: AbstractInstrumentSpec 
    Ei   :: Float64
    L1   :: Float64
    L2   :: Float64
    L3   :: Float64
    Δtp  :: Float64
    Δtc  :: Float64
    Δtd  :: Float64
    Δθ   :: Float64
end

