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
abstract type AbstractConvolutionSpec end

# Add IFFT plan
struct UniformQBroadening <: AbstractConvolutionSpec
    binning  :: UniformBinning

    # Q-broadening
    qfwhm     :: Union{Float64, Array{Float64, 2}}   # Kernel FWHM. Not needed after construction, but useful as a record. 
    qkernel   :: Array{ComplexF64, 3}                # Fourier trasnformed convolution kernel
    qpoints   :: Array{Sunny.Vec3, 3}                # Sampled points in momentum space
    qidcs                                            # Indices corresponding to interior of q bins. Same dimensions as binning.qcenters.

    # E-broadening
    ekernel   :: Sunny.AbstractBroadening            # Sunny broadening to be passed to `Sunny.intensities`
    epoints   :: Vector{Float64}                     # Sampled points in momentum space
    eidcs                                            # Indices corresponding to interior of energy bins. Same dimensions as binning.ecenters.
end


"""
    UniformQBroadening(binning::UniformBinning, qfwhm, ekernel; nperbin, nghosts, nperebin=1)

Generates an intensities calculation specification that performs convolution
over both energy and inverse Anstroms dimensions. The convolution kernel is
Gaussian and separable, i.e., the energy and spatial directions are
uncorrelated. 

- `qfwhm` specifies the spatial convolution kernel. It may be either a single
  number or a matrix. Units are inverse Angstrom.
- `ekernel` is a Sunny energy broadening kernel, e.g. `gaussian(fwhm)`.
- `nperbin` determines how many samples are included per linear dimension of a
  bin. This may be either a single number, or three numbers (one for each linear
  dimension).
- `nghosts` determines a number of padding bins to be included about the given
  binning scheme. These should be included if the size of the q-convolution
  kernel is large relative to the bin size.
- `nperebin` determines the number of samples to be included along the energy
  dimension of each bin. By default just the bin center is used.

"""
function UniformQBroadening(binning::UniformBinning, qfwhm, ekernel; nperbin, nghosts, nperebin=1)
    (; crystal, Δs, qcenters, ecenters, directions) = binning

    # Convert broadening parameters to a covariance matrix.
    σ = qfwhm / 2√(2log(2))
    Σ = isa(σ, Number) ? σ*I(3) : σ

    # Generate nperbin uniformly spaced samples for each bin, including in
    # padding bins.
    (; qpoints, epoints) = sample_binning(binning; nperbin, nghosts)

    # Keep track of which sample points are in which q-bins. *Note that the
    # indices need to be FFT shifted.*
    bounds = [(-Δ/2, Δ/2) for Δ in Δs[1:3]]
    points_shifted = fftshift(qpoints)
    qidcs = map(q0 -> find_points_in_bin(q0, directions, bounds, points_shifted), qcenters)

    # Calculate the convlution kernel.
    binning_center = sum(qcenters)/length(qcenters)
    qkernel = gaussian_md(map(p -> crystal.recipvecs*(p-binning_center), qpoints), [0., 0, 0], Σ)
    qkernel = fft(qkernel, (1, 2, 3))

    # Binning indices for energy axis.
    ΔE = Δs[4]
    eidcs = map(ecenters) do E0
        findall(E -> abs(E0 - E) < ΔE/2, epoints)
    end

    return UniformQBroadening(binning, qfwhm, qkernel, qpoints, qidcs, ekernel, epoints, eidcs)
end