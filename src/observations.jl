struct Observation
    instrument  :: Union{Nothing, AbstractInstrumentSpec}   # Data about the instrument -- stub for now
    binning     :: AbstractBinning                          # Binning used in observation
    ints        :: Array{Float64, 4}                        # Observed intensities
    errs        :: Array{Float64, 4}                        # Errors
    mask        :: Array{Float64, 4}                        # 1.0s and NaNs -- useful for plotting
    mask_idcs   :: Vector{CartesianIndex{4}}                # Indices of non-NaN values in ints and errs
    model       :: Union{Nothing, AbstractDataModel}        # Model of data -- stub for now.
    background  :: Union{Nothing, Function}                 # Model of background, given as a function f(q, E) (not sure worthy of a type yet)
end

function Observation(binning, ints, errs; instrument=nothing, model=nothing, background=nothing, filtersubzeros=false)

    # Check that intensities and errors are compatible with given binning.
    (; ecenters, qcenters) = binning
    composite_size = (length(ecenters), size(qcenters)...)
    @assert composite_size == size(ints) == size(errs) "Size or errors and/or intensities not compatible with given binning scheme"

    # Remove negative values (from background subtraction)
    if filtersubzeros
        for idx in eachindex(ints)
            if ints[idx] < 0
                ints[idx] = errs[idx] = NaN
            end
        end
    end

    # Determine the mask (assuming empty bins are NaNs).
    mask_idcs = findall(val -> !isnan(val), ints)
    mask = NaN * ones(composite_size)
    mask[mask_idcs] .= 1.0

    return Observation(instrument, binning, ints, errs, mask, mask_idcs, model, background)
end
