abstract type AbstractInstrumentSpec end

# abstract type TripleAxisSpec <: AbstractInstrumentSpec end

struct ChopperSpec <: AbstractInstrumentSpec 
    name :: String
    Ei   :: Float64
    L1   :: Float64
    L2   :: Float64
    L3   :: Float64
    Δtp  :: Float64
    Δtc  :: Float64
    Δtd  :: Float64
    Δθ   :: Float64
end


function Base.show(io::IO, obs::ChopperSpec)
    (; name, Ei) = obs
    printstyled(io, "Chopper Instrument\n"; bold=true, color=:underline)
    println(io, "Instrument name: ", name)
    println(io, "Incident energy: $Ei meV")
end