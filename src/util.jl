# Un-normalized Gaussian function where K = Σ⁻¹.
gaussian_func(v, μ, K) = exp(-((v-μ)' * K * (v-μ))/2)

# Normalized Gaussian function. K = Σ⁻¹
function gaussian_md(vs, μ, K)
    k = length(vs[1])
    coef =  sqrt(det(K)/((2π)^k))
    map(v -> coef*exp(-0.5((v-μ)' * K * (v-μ))), vs)
end

# Multidimensional Gaussian samplers. Note that given Σs from hdf5 file are only
# approximately symmetric, so an explicit symmetrization is performed.
sample_q(μ, Σ, N) = rand(MvNormal(μ, (Σ + Σ')/2), N)
sample_q(rng, μ, Σ, N) = rand(rng, MvNormal(μ, (Σ + Σ')/2), N)