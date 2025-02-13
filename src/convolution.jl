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
struct SeparableUniformQBroadening <: AbstractConvolution
    crystal  :: Sunny.Crystal                       # Keep around for recipvecs

    Qfwhm    :: Union{Float64, Array{Float64, 2}}   # Note needed after construction, but useful to be able to report
    Qkernel  :: Array{ComplexF64, 3}                # Fourier trasnformed convolution kernel
    points   :: Array{Sunny.Vec3, 3}                # Sampled points in momentum space
    interior_idcs                                   # Indices
    interior_idcs_ft

    Ekernel 
end


function SeparableUniformQBroadening(bi::BinSpec, bin_fwhm, Ekernel; nsigmas=3, sampdensity=1)
    (; crystal) = bi

    # Set up Gaussian kernel parameters
    σ = bin_fwhm / 2√(2log(2))
    Σ = if isa(σ, Number)
        σ*I(3)
    else
        σ
    end

    # Add heuristics for determining extension of box -- probably an analysis of
    # the extrema as a function of boundary values in local coordinates. For
    # now, just make the volume 27 times larger.
    q0 = [0., 0, 0]
    points = sample_bin_and_surrounds(q0, bi; sigma=Σ, nsigmas, sampdensity)
    interior_idcs = find_points_in_bin(q0, bi, points)
    interior_idcs_ft = find_points_in_bin(q0, bi, fftshift(points))
    Qkernel = gaussian_md(map(p -> crystal.recipvecs*p, points), q0, Σ)
    Qkernel_ft = fft(Qkernel, (1, 2, 3))

    return SeparableUniformQBroadening(crystal, bin_fwhm, Qkernel_ft, points, interior_idcs, interior_idcs_ft, Ekernel)
end

function SeparableUniformQBroadening(binning::UniformBinning, bin_fwhm, Ekernel; nsigmas=3, sampdensity=1)
    (; binspec, bincenters) = binning
    (; crystal, directions, bounds) = binspec

    # Set up Gaussian kernel parameters
    σ = bin_fwhm / 2√(2log(2))
    Σ = if isa(σ, Number)
        σ*I(3)
    else
        σ
    end

    # Add heuristics for determining extension of box -- probably an analysis of
    # the extrema as a function of boundary values in local coordinates. For
    # now, just make the volume 27 times larger.
    q0 = [0., 0, 0]
    points = sample_binning_and_surrounds(binning, Σ; nsigmas, sampdensity)
    interior_idcs = find_points_in_bin(q0, bi, points)
    interior_idcs_ft = find_points_in_bin(q0, bi, fftshift(points))
    Qkernel = gaussian_md(map(p -> crystal.recipvecs*p, points), q0, Σ)
    Qkernel_ft = fft(Qkernel, (1, 2, 3))

    return SeparableUniformQBroadening(crystal, bin_fwhm, Qkernel_ft, points, interior_idcs, interior_idcs_ft, Ekernel)
end


function Base.show(io::IO, broadening_spec::SeparableUniformQBroadening)
    println(io, "SeparableUniformQBroadening. FWHM: ")
end



function intensities_instrument(swt, q; energies, broadening_spec, kwargs...)
    (; points, Qkernel, interior_idcs_ft, crystal, Ekernel) = broadening_spec

    # Calculate intensities for all points in subsuming grid around bin.
    res = Sunny.intensities(swt, map(p -> p + q, points[:]); energies, kernel=Ekernel, kwargs...)
    data = reshape(res.data, (length(energies), size(points)...))

    # Convolve along Q-axes only using an FFT. Unfortunately, energy is the fast
    # axis. Can tweak later.
    data_ft = fft(data, (2, 3, 4))
    for i in axes(data_ft, 1)
        data_ft[i,:,:,:] .*= Qkernel
    end
    data_conv = real.(ifft(data_ft, (2, 3, 4))) ./ prod(size(points))

    # Integrate over those slices that lie within the bin and normalize to the
    # number of samples.
    slice = sum(data_conv[:, interior_idcs_ft], dims=(2,3,4)) ./ length(interior_idcs_ft)

    return Sunny.Intensities(
        crystal,
        Sunny.QPoints([Sunny.Vec3(q)]),
        collect(energies),
        slice,
    )
end