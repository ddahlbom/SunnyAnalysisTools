function chi_square(obs::Observation, calc::ModelCalculation; scale=1.0, intensity_normalization=false)
    (; ints, errs, mask_idcs, background) = obs
    (; data) = calc
    @assert size(data) == size(ints)

    calc_scaled = scale .* data

    if !isnothing(background)
        # TODO: Add addition of background to calculation here
    end

    χ2 = 0.0
    dof = 0
    total_intensity = 0.0
    for idx in mask_idcs
        if !iszero(errs[idx])
            χ2 += (ints[idx] - calc_scaled[idx])^2 / errs[idx]^2
            dof += 1
            total_intensity += ints[idx]
        end
    end

    return intensity_normalization ? (χ2 / total_intensity) / dof : χ2 / dof
end

function similarity_measure1(obs::Observation, calc::ModelCalculation; scale=1.0)
    (; ints, errs, mask_idcs, background) = obs
    (; data) = calc
    @assert size(data) == size(ints)

    calc_scaled = scale .* data

    if !isnothing(background)
        # TODO: Add addition of background to calculation here
    end

    err = 0.0
    for idx in mask_idcs
        if !iszero(errs[idx])
            err += (ints[idx] - calc_scaled[idx])^2
        end
    end

    return err 
end