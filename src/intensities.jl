function calculate_intensities(swtmodel::SWTModel, broadening_spec::StationaryQConvolution; kwargs...)
    (; swt) = swtmodel
    (; qpoints, epoints, qidcs, eidcs, qkernel, ekernel, binning) = broadening_spec
    (; qcenters, Es) = binning

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

    # spec and params
    ModelCalculation(res, binning, broadening_spec, swtmodel.params)
end


function calculate_intensities(swtmodel::SWTModel, broadening_spec::UniformSampling; kwargs...)
    (; swt) = swtmodel
    (; qpoints, epoints, qidcs, eidcs, ekernel, binning) = broadening_spec
    (; qcenters, Es) = binning

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

    # spec and params
    ModelCalculation(res, binning, broadening_spec, swtmodel.params)
end