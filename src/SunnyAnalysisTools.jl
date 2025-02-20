module SunnyAnalysisTools

import Sunny
using LinearAlgebra, Random, LsqFit, FFTW

ENV["JULIA_CONDAPKG_BACKEND"] = "Null"
ENV["JULIA_PYTHONCALL_EXE"] = "/home/uud/miniforge3/envs/neutron/bin/python" # How to automate this?

using PythonCall

################################################################################
# Files and exports 
################################################################################
include("constants.jl")

include("util.jl")

include("instruments.jl")
export CNCSSpec

include("data_modeling.jl")
export gaussian_mixture_model

include("binning.jl")
export UniformBinning

include("calculation_spec.jl")
export Calculation, StationaryQConvolution, UniformSampling

include("resolution.jl")
export energy_resolution_kernel

include("observations.jl")
export Observation

include("models.jl")
export SWTModel

include("intensities.jl")
export CalculationResult, calculate_intensities

include("parsing.jl")
# export read_shiver_ascii

include("errors.jl")
export chi_square


################################################################################
# Only load plotting functions if GLMakie imported. See Sunny.jl as reference
################################################################################
function is_pkg_loaded(pkg::Symbol)
    return any(k -> Symbol(k.name) == pkg, keys(Base.loaded_modules))
end

extension_fns = [
    :GLMakie => [:draw_boundary!, :visualize_binning],
]

for (_pkg, fns) in extension_fns
    for fn in fns
        @eval function $fn end
        @eval export $fn
    end
end

function __init__()
    # Notify user if extension function requires package import
    if isdefined(Base.Experimental, :register_error_hint)
        Base.Experimental.register_error_hint(MethodError) do io, exc, _argtypes, _kwargs
            fn = Symbol(exc.f)
            for (pkg, fns) in extension_fns
                if in(fn, fns) && !is_pkg_loaded(pkg)
                    pkgstr = (pkg == :Makie) ? "a variant of Makie" : "package $pkg"
                    printstyled(io, "\nImport $pkgstr to enable `$fn`.\n"; bold=true)
                end
            end
        end
    end
end

# Access to PlottingExt module for developer convenience
PlottingExt() = Base.get_extension(@__MODULE__, :PlottingExt)

end