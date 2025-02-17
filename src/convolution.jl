################################################################################
# Energy broadening
################################################################################

# This should be rewritten based on more fundamental parameters so it works for
# additional values. See equation from Savici.
function energy_resolution_kernel(spec::CNCSSpec)
    (; Ei) = spec
    fwhm = if Ei ≈ 1.55
        E -> (2.4994e-4)*sqrt((Ei-E)^3 * ( (211.41149*(0.052+0.123*(Ei/(Ei-E))^1.5))^2 + (57.27700*(1.052+0.123*(Ei/(Ei-E))^1.5))^2) ) / 2.355 
    elseif Ei ≈ 2.5
        E -> (2.4994e-4)*sqrt((Ei-E)^3 * ( (183.25988*(0.052+0.123*(Ei/(Ei-E))^1.5))^2 + (57.27700*(1.052+0.123*(Ei/(Ei-E))^1.5))^2) ) / 2.355
    elseif Ei ≈ 6.59
        E -> (2.4994e-4)*sqrt((Ei-E)^3 * ( (124.84362*(0.052+0.123*(Ei/(Ei-E))^1.5))^2 + (57.27700*(1.052+0.123*(Ei/(Ei-E))^1.5))^2) ) / 2.355
    else
        error("No such incident energy setting known. Available eᵢ values are: 1.55, 2.5, 6.59")
        _ -> ()
    end
    return Sunny.NonstationaryBroadening((b, ω) -> exp(-(ω-b)^2/2fwhm(b)^2) / √(2π*fwhm(b)^2))
end


################################################################################
# Q-broadening only
################################################################################

# An AbstractConvolution should provide full specifications for how to treat
# intensities calculations from Sunny, including procedures for both energy
# and momentum convolution.
abstract type AbstractConvolution end

# Add IFFT plan
struct UniformQBroadening <: AbstractConvolution
    binning  :: UniformBinning

    qfwhm     :: Union{Float64, Array{Float64, 2}}   # Note needed after construction, but useful to be able to report
    qkernel   :: Array{ComplexF64, 3}                # Fourier trasnformed convolution kernel
    qpoints   :: Array{Sunny.Vec3, 3}                # Sampled points in momentum space
    epoints   :: Vector{Float64}                # Sampled points in momentum space
    qidcs    
    eidcs

    ekernel  :: Sunny.AbstractBroadening
end


function UniformQBroadening(binning::UniformBinning, bin_fwhm, ekernel; nperbin, nghosts)
    (; crystal, Δs, qcenters, ecenters, directions) = binning

    # Convert broadening parameters to a covariance matrix.
    σ = bin_fwhm / 2√(2log(2))
    Σ = isa(σ, Number) ? σ*I(3) : σ

    # Generate nperbin uniformly spaced samples for each bin, including in
    # padding bins (the number of which are determined by nghosts)
    (; qpoints, epoints) = sample_binning(binning; nperbin, nghosts)

    # Keep track of which points are in which q-bins. Note that the indices need to be FFT shifted.
    bounds = [(-Δ/2, Δ/2) for Δ in Δs[1:3]]
    points_shifted = fftshift(qpoints)
    qidcs = map(q0 -> find_points_in_bin(q0, directions, bounds, points_shifted), qcenters)

    # Binning indices for energy axis.
    ΔE = Δs[4]
    eidcs = map(ecenters) do E0
        findall(E -> abs(E0 - E) < ΔE/2, epoints)
    end

    # Generate the convlution kernel.
    binning_center = sum(qcenters)/length(qcenters)
    qkernel = gaussian_md(map(p -> crystal.recipvecs*(p-binning_center), qpoints), [0., 0, 0], Σ)
    qkernel_ft = fft(qkernel, (1, 2, 3))

    return UniformQBroadening(binning, bin_fwhm, qkernel_ft, qpoints, epoints, qidcs, eidcs, ekernel)
end


function intensities(swt, broadening_spec::UniformQBroadening; kwargs...)
    (; qpoints, epoints, qidcs, eidcs, qkernel, ekernel, binning) = broadening_spec
    (; qcenters, ecenters) = binning

    # Calculate intensities for all points in subsuming grid around bin.
    res = Sunny.intensities(swt, qpoints[:]; energies=epoints, kernel=ekernel, kwargs...)
    data = reshape(res.data, (length(epoints), size(qpoints)...))

    # Convolve along Q-axes only using an FFT. Unfortunately, energy is the fast
    # axis. 
    data_ft = fft(data, (2, 3, 4))
    for i in axes(data_ft, 1)
        data_ft[i,:,:,:] .*= qkernel
    end
    data_conv = real.(ifft(data_ft, (2, 3, 4))) ./ prod(size(qpoints))

    # Sum over samples that lie within each bin and normalize by number of
    # samples.
    res = zeros(length(ecenters), size(qcenters)...)
    for i in CartesianIndices(qcenters), j in eachindex(ecenters)
        for (ei, qi) in Iterators.product(eidcs[j], qidcs[i])
            res[j, i] += data_conv[ei,qi]
        end
        res[j, i] /= length(eidcs[j]) * length(qidcs[i])
    end

    return res
end