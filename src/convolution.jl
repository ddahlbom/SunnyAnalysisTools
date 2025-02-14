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

    qfwhm    :: Union{Float64, Array{Float64, 2}}   # Note needed after construction, but useful to be able to report
    qkernel  :: Array{ComplexF64, 3}                # Fourier trasnformed convolution kernel
    points   :: Array{Sunny.Vec3, 3}                # Sampled points in momentum space
    interior_idcs                                   # Indices
    interior_idcs_ft

    ekernel  :: Sunny.AbstractBroadening
end


function UniformQBroadening(binning::UniformBinning, bin_fwhm, ekernel; nperbin, nghosts)
    (; crystal, Δs, bincenters, directions) = binning

    # Convert broadening parameters to a covariance matrix.
    σ = bin_fwhm / 2√(2log(2))
    Σ = isa(σ, Number) ? σ*I(3) : σ

    # Generate nperbin uniformly spaced samples for each bin, including in
    # padding bins (the number of which are determined by nghosts)
    points = sample_binning(binning; nperbin, nghosts)
    points_shifted = fftshift(points)

    # Keep track of which points are in which bins.
    bounds = [(-Δ/2, Δ/2) for Δ in Δs[1:3]]
    interior_idcs = map(q0 -> find_points_in_bin(q0, directions, bounds, points), bincenters)
    interior_idcs_ft = map(q0 -> find_points_in_bin(q0, directions, bounds, points_shifted), bincenters)

    # Generate the convlution kernel.
    binning_center = sum(bincenters)/length(bincenters)
    qkernel = gaussian_md(map(p -> crystal.recipvecs*(p-binning_center), points), [0., 0, 0], Σ)
    qkernel_ft = fft(qkernel, (1, 2, 3))

    return UniformQBroadening(binning, bin_fwhm, qkernel_ft, points, interior_idcs, interior_idcs_ft, ekernel)
end


function intensities(swt, broadening_spec::UniformQBroadening; kwargs...)
    (; points, qkernel, interior_idcs_ft, ekernel, binning) = broadening_spec
    (; bincenters, Ecenters, crystal) = binning

    energies = Ecenters

    # Calculate intensities for all points in subsuming grid around bin.
    res = Sunny.intensities(swt, points[:]; energies, kernel=ekernel, kwargs...)
    data = reshape(res.data, (length(energies), size(points)...))

    # Convolve along Q-axes only using an FFT. Unfortunately, energy is the fast
    # axis. Can tweak later.
    data_ft = fft(data, (2, 3, 4))
    for i in axes(data_ft, 1)
        data_ft[i,:,:,:] .*= qkernel
    end
    data_conv = real.(ifft(data_ft, (2, 3, 4))) ./ prod(size(points))

    # Integrate over those slices that lie within the bin and normalize to the
    # number of samples.
    res = zeros(length(energies), size(bincenters)...)
    for i in CartesianIndices(size(bincenters))  
        res[:,i] = sum(data_conv[:, interior_idcs_ft[i]], dims=(2,3,4)) ./ length(interior_idcs_ft[i])
    end

    return res
end