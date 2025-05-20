function calculate_intensities(swtmodel::SWTModel, broadening_spec::StationaryQConvolution; kwargs...)
    (; swt) = swtmodel
    (; qpoints, epoints, qidcs, eidcs, qkernel, ekernel, binning) = broadening_spec
    (; qcenters, Es, binvol, crystvol) = binning

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
    res = zeros(length(Es), size(qcenters)...)
    for i in CartesianIndices(qcenters), j in eachindex(Es)
        for (ei, qi) in Iterators.product(eidcs[j], qidcs[i])
            res[j, i] += data_conv[ei,qi]
        end
        res[j, i] /= length(eidcs[j]) * length(qidcs[i])
    end
    res .*= binvol

    # spec and params
    ModelCalculation(res, binning, broadening_spec, swtmodel.params)
end


function calculate_intensities(swtmodel::SWTModel, broadening_spec::UniformSampling; kwargs...)
    (; swt) = swtmodel
    (; qpoints, epoints, qidcs, eidcs, ekernel, binning) = broadening_spec
    (; qcenters, Es, binvol) = binning

    # Calculate intensities for all points in subsuming grid around bin.
    res = Sunny.intensities(swt, qpoints[:]; energies=epoints, kernel=ekernel, kwargs...)
    data = reshape(res.data, (length(epoints), size(qpoints)...))

    # Sum over samples that lie within each bin and normalize by number of
    # samples.
    res = zeros(length(Es), size(qcenters)...)
    for i in CartesianIndices(qcenters), j in eachindex(Es)
        for (ei, qi) in Iterators.product(eidcs[j], qidcs[i])
            res[j, i] += data[ei,qi]
        end
        res[j, i] /= length(eidcs[j]) * length(qidcs[i])
    end
    res .*= binvol

    # spec and params
    ModelCalculation(res, binning, broadening_spec, swtmodel.params)
end

gaussian_func(v, μ, K) = exp(-((v-μ)' * K * (v-μ))/2)
sample_q(μ, Σ, N) = rand(MvNormal(μ, (Σ + Σ')/2), N)
sample_q(rng, μ, Σ, N) = rand(rng, MvNormal(μ, (Σ + Σ')/2), N)

function tax_convolved_intensity_mc(intfunc, qe0, K, nsamps)
    Σ = inv(K)
    qes = sample_q(qe0, Σ, nsamps)
    hkls = [Sunny.Vec3(qe[1], qe[2], qe[3]) for qe in eachcol(qes)]
    hkls = reshape(hkls, length(hkls))

    # (; data, disp) = intensities_bands(swt, hkls)
    (; data, disp) = intfunc(hkls)

    cumval = 0.0
    for (iq, q) in enumerate(hkls), iband in axes(disp, 1)
        qe = SVector{4, Float64}(q..., disp[iband, iq])
        cumval += data[iband, iq] * gaussian_func(qe, qe0, K) # Shouldn't have to multiply by gaussian_func...
        # cumval += data[iband, iq] 
    end

    return cumval/nsamps 
end


function calculate_intensities(swtmodel::SWTModel, taxspec::TripleAxisMC{1}; kwargs...)
    (; path, N, Ks) = taxspec 
    (; HKLs, Es, projection) = path

    buf = zero(path.Es)
    intfunc(hkls) = Sunny.intensities_bands(swtmodel.swt, hkls; kwargs...)
    for (n, (HKL, K, E)) in enumerate(zip(HKLs, Ks, Es))
        q = projection*HKL
        buf[n] = tax_convolved_intensity_mc(intfunc, SVector{4, Float64}(q..., E), K, N)
    end
    return buf
end

function calculate_intensities(swtmodel::SWTModel, taxspec::TripleAxisMC{2}; kwargs...)
    (; path, N, Ks) = taxspec 
    (; HKLs, Es, projection) = path

    buf = zeros(length(HKLs), length(path.Es))
    intfunc(hkls) = Sunny.intensities_bands(swtmodel.swt, hkls; kwargs...)
    #for (n, (HKL, K, E)) in enumerate(zip(HKLs, Ks, Es))
    for (n, ((j, HKL), (k, E))) in enumerate(Iterators.product(enumerate(HKLs), enumerate(Es)))
        K = Ks[n]
        q = projection*HKL
        buf[j,k] = tax_convolved_intensity_mc(intfunc, SVector{4, Float64}(q..., E), K, N)
    end
    return buf
end

function calculate_intensities(func::Function, taxspec::TripleAxisMC{2}; kwargs...)
    (; path, N, Ks) = taxspec 
    (; HKLs, Es, projection) = path

    buf = zeros(length(HKLs), length(path.Es))
    function intfunc(hkls) 
        N = length(hkls)
        data = ones(1, N)
        disp = zeros(1, N) 

        for (n, hkl) in enumerate(hkls)
            disp
        end

        return (; data, disp)
    end
    #for (n, (HKL, K, E)) in enumerate(zip(HKLs, Ks, Es))
    for (n, ((j, HKL), (k, E))) in enumerate(Iterators.product(enumerate(HKLs), enumerate(Es)))
        K = Ks[n]
        q = projection*HKL
        buf[j,k] = tax_convolved_intensity_mc(intfunc, SVector{4, Float64}(q..., E), K, N)
    end
    return buf
end
