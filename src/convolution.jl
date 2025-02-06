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
abstract type AbstractQBroadening end

# Add FFT plan
struct UniformQBroadening <: AbstractQBroadening
    fwhm    :: Union{Float64, Array{Float64, 2}}
    kernel  :: Array{ComplexF64, 3}
    points  :: Array{Sunny.Vec3, 3}
    crystal :: Sunny.Crystal
    interior_idcs   # Ultimately don't need to keep these -- keep around for debugging.
    interior_idcs_ft
end

function widen_bounds_by_factor(bounds, factor)
    Δs = [(b[2] - b[1]) for b in bounds]
    centers = [(b[1] + b[2])/2 for b in bounds]
    return [(c - factor*Δ/2, c + factor*Δ/2) for (c, Δ) in zip(centers, Δs)]
end


function UniformQBroadening(bi::BinInfo, fwhm; nsigmas=3, sampdensity=1)
    (; crystal) = bi

    # Set up Gaussian kernel parameters
    σ = fwhm / 2√(2log(2))
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
    kernel = gaussian_md(map(p -> crystal.recipvecs*p, points), q0, Σ)
    kernel_ft = fft(kernel, (1, 2, 3))

    return UniformQBroadening(fwhm, kernel_ft, points, crystal, interior_idcs, interior_idcs_ft)
end


function intensities_instrument(swt, q, qbroadeningspec; energies, energy_kernel, kwargs...)
    (; points, kernel, interior_idcs_ft, crystal) = qbroadeningspec

    # Calculate intensities for all points in subsuming grid around bin.
    res = Sunny.intensities(swt, map(p -> p + q, points[:]); energies, kernel=energy_kernel, kwargs...)
    data = reshape(res.data, (length(energies), size(points)...))

    # Convolve along Q-axes only using an FFT. Unfortunately, energy is the fast
    # axis. Can tweak later.
    data_ft = fft(data, (2, 3, 4))
    for i in axes(data_ft, 1)
        data_ft[i,:,:,:] .*= kernel
    end
    data_conv = real.(ifft(data_ft, (2, 3, 4))) ./ prod(size(points))

    # Integrate over those slices that lie within the bin.
    slice = sum(data_conv[:, interior_idcs_ft], dims=(2,3,4)) ./ length(interior_idcs_ft)

    return Sunny.Intensities(
        crystal,
        Sunny.QPoints([Sunny.Vec3(q)]),
        collect(energies),
        slice,
    )
end