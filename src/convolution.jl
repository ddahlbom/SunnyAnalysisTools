################################################################################
# Energy broadening
################################################################################

# This should be rewritten based on more fundamental parameters so it works for
# additional values. See equation from Savici.
function E_resolution_kernel(spec::CNCSSpec)
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
    fwhm    :: Float64
    spacing :: Float64
    kernel  :: Array{ComplexF64, 3}
    points  :: Array{Sunny.Vec3, 3}
    interior_idcs
end

function widen_bounds_by_factor(bounds, factor)
    Δs = [(b[2] - b[1]) for b in bounds]
    centers = [(b[1] + b[2])/2 for b in bounds]
    return [(c - factor*Δ/2, c + factor*Δ/2) for (c, Δ) in zip(centers, Δs)]
end


function UniformQBroadening(bincenter, bi::BinInfo, fwhm, spacing=0.1)
    (; crystal) = bi

    # Set up Gaussian kernel parameters
    σ = fwhm / 2√(2log(2))
    K = if length(σ) == 1
        Matrix(1/σ * I(3))
    else 
        diagm(1 ./ σ)
    end

    # Add heuristics for determining extension of box -- probably an analysis of
    # the extrema as a function of boundary values in local coordinates.
    (; crystal, directions, bounds) = bi 
    bounds_new = widen_bounds_by_factor(bounds, 2.0)
    bi_convolution = BinInfo(crystal, directions, bounds_new)
    points = uniform_grid_from_bin(bincenter, bi_convolution, spacing)
    interior_idcs = find_points_in_bin(bincenter, bi, points)
    kernel = gaussian_md(points, crystal.recipvecs*bincenter, K)
    kernel_ft = fft(kernel)

    return UniformQBroadening(fwhm, spacing, kernel_ft, points, interior_idcs)
end


function instrument_intensities(swt, qbroadeningspec; energies, measure, E_kernel)
    (; points, kernel_ft, interior_idcs) = qbroadeningspec

    # Calculate intensities for all points in subsuming grid around bin
    res = intensities(swt, points[:]; energies, measure, kernel=E_kernel)
    data = reshape(res.data, (length(energies), size(points)...))

    # Convolve along Q-axes only using an FFT
    data_ft = fft(data, (2, 3, 4))
    for i in axes(data_ft, 1)
        data_ft[i,:,:,:] .*= kernel_ft
    end
    data_conv = real.(ifft(data_ft, (2, 3, 4))) ./ prod(size(points))^2

    # Integrate over those slices that lie within the bin
    slice = sum(data_conv[:, interior_idcs], dims=(2,3,4)) ./ length(interior_idcs)
    return dropdims(slice; dims=(2,))
end