# Un-normalized Gaussian function where K = Σ⁻¹.
gaussian_func(v, μ, K) = exp(-((v-μ)' * K * (v-μ))/2)

# Normalized Gaussian function. 
function gaussian_md(vs, μ, Σ)
    k = length(vs[1])
    coef =  1/sqrt(det(Σ)*(2π)^k)
    map(v -> coef*exp(-((v-μ)' * inv(Σ) * (v-μ))/2), vs)
end

# Multidimensional Gaussian samplers. Note that given Σs from hdf5 file are only
# approximately symmetric, so an explicit symmetrization is performed.
sample_q(μ, Σ, N) = rand(MvNormal(μ, (Σ + Σ')/2), N)
sample_q(rng, μ, Σ, N) = rand(rng, MvNormal(μ, (Σ + Σ')/2), N)