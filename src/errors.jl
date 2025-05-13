"""
    chi_square(obs::TimeOfFlightObservation, calc::ModelCalculation; scale=1.0, ndof=1)

Returns 

```math
    χ² = ∑_{j} (I_{exp}_j - I_{calc}_j)² / (ν σ_j²)
```

where `ν` is the number of degrees of freedom in the model, set with the keyword
`ndof`.
"""
function chi_square(obs::TimeOfFlightObservation, calc::ModelCalculation; scale=1.0, ndof=1)
    (; ints, errs, mask_idcs, background) = obs
    (; data) = calc
    @assert size(data) == size(ints)

    calc_scaled = scale .* data

    if !isnothing(background)
        warn("Background corrections not implemented")
    end

    χ2 = 0.0
    for idx in mask_idcs
        if !iszero(errs[idx])
            χ2 += (ints[idx] - calc_scaled[idx])^2 / errs[idx]^2
        end
    end

    return χ2 / ndof 
end

"""
    weighted_residual(obs::TimeOfFlightObservation, calc::ModelCalculation; scale=1.0)

Returns 

```math
    R = ∑_{j} (I_{exp}_j - I_{calc}_j)² / (I_{calc}_j)^2
```

where `ν` is the number of degrees of freedom in the model, set with the keyword
`ndof`.
"""
function weighted_residual(obs::TimeOfFlightObservation, calc::ModelCalculation; scale=1.0)
    (; ints, errs, mask_idcs, background) = obs
    (; data) = calc
    @assert size(data) == size(ints)

    calc_scaled = scale .* data

    if !isnothing(background)
        warn("Background corrections not implemented")
    end

    residual = 0.0
    for idx in mask_idcs
        if !iszero(errs[idx])
            residual += (ints[idx] - calc_scaled[idx])^2 / (calc_scaled[idx]^2)
        end
    end

    return residual 
end

"""
    similarity_measure(obs::TimeOfFlightObservation, calc::ModelCalculation; scale=1.0)

Returns 

```math
    R_{wp} = √[∑_{j} ((I_{exp}_j - I_{calc}_j) / σ_{exp})² / (I_{exp}_j / σ_{exp})²]
```

where `ν` is the number of degrees of freedom in the model, set with the keyword
`ndof`.
"""
function similarity_measure(obs::TimeOfFlightObservation, calc::ModelCalculation; scale=1.0)
    (; ints, errs, mask_idcs, background) = obs
    (; data) = calc
    @assert size(data) == size(ints)

    calc_scaled = scale .* data

    if !isnothing(background)
        warn("Background corrections not implemented")
    end

    Rwp = 0.0
    for idx in mask_idcs
        if !iszero(errs[idx])
            Rwp += ((ints[idx] - calc_scaled[idx])/errs[idx])^2 / (ints[idx]/errs[idx])^2
        end
    end

    return √Rwp 
end


"""
    squared_difference(obs::TimeOfFlightObservation, calc::ModelCalculation; scale=1.0)

Returns 

```math
    χ² = ∑_{j} (I_{exp}_j - I_{calc}_j)²
```

where `ν` is the number of degrees of freedom in the model, set with the keyword
`ndof`.
"""
function squared_difference(obs::TimeOfFlightObservation, calc::ModelCalculation; scale=1.0)
    (; ints, errs, mask_idcs, background) = obs
    (; data) = calc
    @assert size(data) == size(ints)

    calc_scaled = scale .* data

    if !isnothing(background)
        warn("Background corrections not implemented")
    end

    err = 0.0
    for idx in mask_idcs
        if !iszero(errs[idx])
            err += (ints[idx] - calc_scaled[idx])^2
        end
    end

    return err
end
