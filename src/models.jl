struct SWTModel
    swt
    params :: Union{Nothing, Dict{Any, Any}}
end

function Base.show(io::IO, ::SWTModel)
    printstyled(io, "Spin Wave Model\n"; bold=true, color=:underline)
end