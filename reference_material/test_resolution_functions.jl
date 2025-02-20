using Sunny, LinearAlgebra, GLMakie, SunnyAnalysisTools, LsqFit
import SunnyAnalysisTools: energy_resolution, dQx, dQy, dQ

ENV["JULIA_CONDAPKG_BACKEND"] = "Null"
ENV["JULIA_PYTHONCALL_EXE"] = "/home/uud/miniforge3/envs/neutron/bin/python" # How to automate this?
using PythonCall


mantid = pyimport("mantid")
Instruments = pyimport("pychop.Instruments")

CNCS = Instruments.Instrument("CNCS")
Ei = 6.59

# Instrument Parameters
distances = CNCS.chopper_system.getDistances()
x0, xa, x1, x2, xm = map(val -> pyconvert(Float64, val), distances) 
moderator_width  = CNCS.moderator.getWidthSquared(Ei) |> x -> pyconvert(Float64, x) |> sqrt
chopper_width = CNCS.chopper_system.getWidthSquared(Ei)[0] |> x -> pyconvert(Float64, x) |> sqrt

L1 = x0 - xm
L2 = x1
L3 = x2
Δtp = moderator_width * (1 - (xm/x0))
Δtc = chopper_width 
Δθ = 1.5 * (pi / 180.0)*0.5
Δtd = 0.0

################################################################################

using Statistics

begin
Es = range(-Ei, Ei, 100)
Qs = range(0.8, 1.2, 10)

dE_data, dQ_data = [], []
for E in Es
    dE_list, dQ_list = [], []
    for Q in Qs
        dE = SunnyAnalysisTools.energy_resolution_full(Ei, E, L1, L2, L3, Δtp, Δtd, Δtc)
        dQval = dQ(Ei, E, Q, L1, L2, L3, Δtp, Δtd, Δtc, Δθ)
        push!(dE_list, dE)
        push!(dQ_list, dQval)
    end
    push!(dE_data, [E, mean(dE_list), std(dE_list)])
    push!(dQ_data, [E, mean(dQ_list), std(dQ_list)])
end
dE_data = hcat(dE_data...)'
dQ_data = hcat(dQ_data...)'

@. poly3(x, p) = p[1]*x^3 + p[2]*x^2 + p[3]*x + p[4]

fit_dE = curve_fit(poly3, dE_data[:,1], dE_data[:,2], ones(4))
fit_dQ = curve_fit(poly3, dQ_data[:,1], dQ_data[:,2], ones(4))

# Plot energy resolution
fig = let
    fig = Figure()
    ax = Axis(fig[1,1]; title="Energy Resolution", xticks=-6:1:6, xlabel="Energy Transfer (meV)", ylabel="dE (meV)")
    lines!(ax, dE_data[:,1], poly3(dE_data[:,1], fit_dE.param); color=:orange)
    scatter!(ax, dE_data[:,1], dE_data[:,2])
    fig
end
fig
end

# Plot momentum resolution vs energy transfer
fig = let
    fig = Figure()
    ax = Axis(fig[1,1]; xticks = -06:1:6, xlabel="Energy Transfer", ylabel="dQ (1/Å)")
    # ylims!(ax, 0.005, 0.08)
    lines!(ax, dQ_data[:,1], poly3(dQ_data[:,1], fit_dQ.param); color=:orange, label="fit", linestyle=:dash)
    scatter!(ax, dQ_data[:,1], dQ_data[:,2]; markersize=12, label="dA data")
    errorbars!(ax, dQ_data[:,1], dQ_data[:,2], dQ_data[:,3])
    axislegend(ax)
    fig
end

a = b = 5.711
c = 21.16411
latvecs = lattice_vectors(a, b, c, 90, 90, 120)
crystal = Crystal(latvecs, [[0, 0, 0]])
recipvecs = crystal.recipvecs
a_star, b_star, c_star = SunnyAnalysisTools.reciprocal_lattice(a, b, c, 90, 90, 120)


dQa_rlu = dQ_data[:,2] / a_star
dQb_rlu = dQ_data[:,2] / b_star
dQc_rlu = dQ_data[:,2] / c_star

fig = let
    fig = Figure()
    ax = Axis(fig[1,1]; xlabel="Energy transfer (meV)", ylabel="dQ (RLU)")
    markersize = 12
    scatter!(ax, dQ_data[:, 1], dQa_rlu; markersize, label="dQ along a*")
    errorbars!(ax, dQ_data[:, 1], dQa_rlu, dQ_data[:,3] / a_star)
    scatter!(ax, dQ_data[:, 1], dQb_rlu; markersize, label="dQ along b*")
    errorbars!(ax, dQ_data[:, 1], dQb_rlu, dQ_data[:,3] / b_star)
    scatter!(ax, dQ_data[:, 1], dQc_rlu; markersize, label="dQ along c*")
    errorbars!(ax, dQ_data[:, 1], dQc_rlu, dQ_data[:,3] / c_star)
    axislegend(ax)
    fig
end