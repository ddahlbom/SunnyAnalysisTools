module SunnyHelpersORNL

import Sunny
using LinearAlgebra, Random, LsqFit, FFTW

################################################################################
# Files and exports 
################################################################################
include("util.jl")

include("instrument_configs.jl")
export CNCSSpec

include("data_modeling.jl")
export gaussian_mixture_model

include("binning.jl")
export BinSpec, sample_bin_uniform_abs, sample_bin_uniform_rlu, sample_bin_mc, find_points_in_bin, uniform_grid_from_bin

include("convolution.jl")
export energy_resolution_kernel, UniformQBroadening, intensities


################################################################################
# Only load plotting functions if GLMakie imported. See Sunny.jl as reference
################################################################################
function is_pkg_loaded(pkg::Symbol)
    return any(k -> Symbol(k.name) == pkg, keys(Base.loaded_modules))
end

extension_fns = [
    :GLMakie => [:draw_boundary!, :visualize_bin, :visualize_binning],
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