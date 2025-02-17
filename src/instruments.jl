abstract type AbstractInstrumentSpec end

abstract type ChopperSpec <: AbstractInstrumentSpec end
abstract type TripleAxisSpec <: AbstractInstrumentSpec end

struct CNCSSpec <: ChopperSpec
    chopper    :: Symbol  # :hiflux, :hires, :intermediate
    Ei         :: Float64
    disk_freq  :: Float64
    fermi_freq :: Float64
end

function CNCSSpec(Ei; chopper=:hiflux, disk_freq=300., fermi_freq=60.)
    # Eventually use this space to calculate the actually physically relevant parameters.
    CNCSSpec(chopper, Ei, disk_freq, fermi_freq)
end
