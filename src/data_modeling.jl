abstract type AbstractDataModel end

function gaussian_mixture_model(x, p, N=4)
    g(x, μ, σ) = (1/(σ *√(2π)))*exp.(-0.5((x - μ)/σ)^2)
    offset = p[end]
    out = zero(x) 
    for i in 1:N
        a, μ, σ = p[(i-1)*3+1:i*3] 
        out .+= a * g.(x, μ, σ)
    end
    return out .+= offset
end