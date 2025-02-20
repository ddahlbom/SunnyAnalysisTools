import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from mantid.simpleapi import *
from pychop.Instruments import Instrument


# Constants
m_n = 1.675e-27  # Neutron mass in kg
hbar = 1.0545718e-34  # Reduced Planck constant in J·s
h = 6.626e-34  # Planck constant in J·s
meV_to_J = 1.60218e-22  # Conversion factor from eV to Joules
angstrom_to_meter = 1e10  # Conversion factor for Q from Å^-1 to m^-1

# Functions
def energy_to_velocity(E_meV):
    return np.sqrt(2 * E_meV * meV_to_J / m_n)

def velocity_to_wavelength(v):
    return h / (m_n * v) * 1e10  # Wavelength in Å

def calculate_two_theta(Ei, E, Q):
    Ef = Ei - E
    v_i = energy_to_velocity(Ei)
    v_f = energy_to_velocity(Ef)
    Q = Q * angstrom_to_meter

    cos_2theta = (v_i**2 + v_f**2 - (Q**2 * hbar**2 / m_n**2)) / (2 * v_i * v_f)
    cos_2theta = np.clip(cos_2theta, -1, 1)  # Ensure within valid range
    return np.arccos(cos_2theta)

def energy_resolution_explicit(E_i, E, L1, L2, L3, delta_tp, delta_tc, delta_td):
    E_i_J = E_i * meV_to_J
    E_f_J = (E_i - E) * meV_to_J
    vi_cubed = (2 * E_i_J / m_n)**(3/2)
    vf_cubed = (2 * E_f_J / m_n)**(3/2)
    
    term_tp = (vi_cubed / L1 + vf_cubed * L2 / (L1 * L3))**2 * delta_tp**2
    term_tc = (vi_cubed / L1 + vf_cubed * (L1 + L2) / (L1 * L3))**2 * delta_tc**2
    term_td = (vf_cubed / L3)**2 * delta_td**2

    delta_E_J = m_n * np.sqrt(term_tp + term_tc + term_td)
    return delta_E_J * 6.24151e21  # Convert back to meV

def dQ(Ei, E, Q, L1, L2, L3, delta_tp, delta_td, delta_tc, delta_theta):
    dQx_val = dQx(Ei, E, Q, L1, L2, L3, delta_tp, delta_td, delta_tc, delta_theta)
    dQy_val = dQy(Ei, E, Q, L1, L2, L3, delta_tp, delta_td, delta_tc, delta_theta)

    two_theta = calculate_two_theta(Ei, E, Q)
    Qx = (m_n / hbar) * (energy_to_velocity(Ei) - energy_to_velocity(Ei - E) * np.cos(two_theta)) / angstrom_to_meter
    Qy = (m_n / hbar) * (-energy_to_velocity(Ei - E) * np.sin(two_theta)) / angstrom_to_meter
    
    return (1 / Q) * np.sqrt((Qx**2 * dQx_val**2) + (Qy**2 * dQy_val**2))

def dQx(Ei, E, Q, L1, L2, L3, delta_tp, delta_td, delta_tc, delta_theta):
    Ef = Ei - E
    v_i = energy_to_velocity(Ei)
    v_f = energy_to_velocity(Ef)
    two_theta = calculate_two_theta(Ei, E, Q)

    cos_2theta = np.cos(two_theta)
    sin_2theta = np.sin(two_theta)

    term1 = (v_i**2 + v_f**2 * (L2 / L3) * cos_2theta)**2 * delta_tp**2 / L1**2
    term2 = (v_i**2 + v_f**2 * ((L1 + L2) / L3) * cos_2theta)**2 * delta_tc**2 / L1**2
    term3 = (v_f**2 / L3 * cos_2theta)**2 * delta_td**2
    term4 = (v_f**2 * sin_2theta**2) * delta_theta**2

    return (m_n / hbar) * np.sqrt(term1 + term2 + term3 + term4) / angstrom_to_meter

def dQy(Ei, E, Q, L1, L2, L3, delta_tp, delta_td, delta_tc, delta_theta):
    Ef = Ei - E
    v_i = energy_to_velocity(Ei)
    v_f = energy_to_velocity(Ef)
    two_theta = calculate_two_theta(Ei, E, Q)

    cos_2theta = np.cos(two_theta)
    sin_2theta = np.sin(two_theta)

    term1 = (v_f**2 * (L2 / (L1 * L3)) * sin_2theta)**2 * delta_tp**2
    term2 = (v_f**2 * ((L1 + L2) / (L1 * L3)) * sin_2theta)**2 * delta_tc**2
    term3 = (v_f**2 / L3 * sin_2theta)**2 * delta_td**2
    term4 = (v_f**2 * cos_2theta**2) * delta_theta**2

    return (m_n / hbar) * np.sqrt(term1 + term2 + term3 + term4) / angstrom_to_meter





# Parameters
E_i = 6.59
E_values = np.linspace(-E_i, E_i, 100)
Q_values = np.linspace(0.8, 1.2, 10)
CNCS = Instrument('CNCS')
x0, xa, x1, x2, xm = CNCS.chopper_system.getDistances()


#E = np.arange(-E_i,E_i, 0.01*E_i)  # meV (energy transfer)
#E = 0
L1 = x0 - xm # m
L2 = x1  # m
L3 = x2  # m
delta_tp = np.sqrt(CNCS.moderator.getWidthSquared(E_i)) * (1 - (xm / x0))
delta_tc = np.sqrt(CNCS.chopper_system.getWidthSquared(E_i))[0]
delta_theta = 1.5 * np.pi / 180.0*0.5  # radians
delta_td = 0

# Compute dE and dQ
dE_data, dQ_data = [], []
for E in E_values:
    dE_list, dQ_list = [], []
    for Q in Q_values:
        dE = energy_resolution_explicit(E_i, E, L1, L2, L3, delta_tp, delta_tc, delta_td)
        dQ_val = dQ(E_i, E, Q, L1, L2, L3, delta_tp, delta_td, delta_tc, delta_theta)
        if np.isfinite(dE) and np.isfinite(dQ_val):
            dE_list.append(dE)
            dQ_list.append(dQ_val)
    if dE_list and dQ_list:
        dE_data.append((E, np.mean(dE_list), np.std(dE_list)))
        dQ_data.append((E, np.mean(dQ_list), np.std(dQ_list)))

dE_data, dQ_data = np.array(dE_data), np.array(dQ_data)

# Fit to third-order polynomial
def poly3(x, a, b, c, d):
    return a * x**3 + b * x**2 + c * x + d

popt_dE, _ = curve_fit(poly3, dE_data[:, 0], dE_data[:, 1])
popt_dQ, _ = curve_fit(poly3, dQ_data[:, 0], dQ_data[:, 1])

# Print polynomial coefficients
print("dE fit coefficients:", popt_dE)
print("dQ fit coefficients:", popt_dQ)

# Plot Energy Resolution vs. Energy Transfer
plt.figure()
plt.errorbar(dE_data[:, 0], dE_data[:, 1], yerr=dE_data[:, 2], fmt='o', label="dE data")
plt.plot(dE_data[:, 0], poly3(dE_data[:, 0], *popt_dE), label="Fit: dE", linestyle="--")
plt.xlabel("Energy Transfer (meV)")
plt.ylabel("dE (meV)")
plt.legend()
plt.title("Energy Resolution vs. Energy Transfer")
plt.grid()
plt.show()

# Plot Momentum Resolution vs. Energy Transfer
plt.figure()
plt.errorbar(dQ_data[:, 0], dQ_data[:, 1], yerr=dQ_data[:, 2], fmt='o', label="dQ data")
plt.plot(dQ_data[:, 0], poly3(dQ_data[:, 0], *popt_dQ), label="Fit: dQ", linestyle="--")
plt.xlabel("Energy Transfer (meV)")
plt.ylabel("dQ (1/Å)")
plt.legend()
plt.title("Momentum Resolution vs. Energy Transfer")
plt.grid()
plt.show()

# Function to compute reciprocal lattice basis vectors from direct lattice parameters
def reciprocal_lattice(a, b, c, alpha, beta, gamma):
    alpha, beta, gamma = np.radians([alpha, beta, gamma])  # Convert to radians

    # Volume of the direct unit cell
    V = a * b * c * np.sqrt(
        1 - np.cos(alpha)**2 - np.cos(beta)**2 - np.cos(gamma)**2
        + 2 * np.cos(alpha) * np.cos(beta) * np.cos(gamma)
    )

    # Reciprocal lattice parameters
    a_star = (b * c * np.sin(alpha)) / V
    b_star = (a * c * np.sin(beta)) / V
    c_star = (a * b * np.sin(gamma)) / V

    return 2 * np.pi * a_star, 2 * np.pi * b_star, 2 * np.pi * c_star  # Include 2π factor

# User-defined lattice parameters
a, b, c = 5.711, 5.711, 21.16411
alpha, beta, gamma = 90, 90, 120

# Compute reciprocal lattice unit vectors with the correct 2π factor
a_star, b_star, c_star = reciprocal_lattice(a, b, c, alpha, beta, gamma)

# Convert dQ from Å⁻¹ to RLU using the corrected conversion
dQ_a_star = dQ_data[:, 1] / a_star
dQ_b_star = dQ_data[:, 1] / b_star
dQ_c_star = dQ_data[:, 1] / c_star

# Plot dQ in reciprocal lattice units
plt.figure()
plt.errorbar(dQ_data[:, 0], dQ_a_star, yerr=dQ_data[:, 2] / a_star, fmt='o', label="dQ along a*")
plt.errorbar(dQ_data[:, 0], dQ_b_star, yerr=dQ_data[:, 2] / b_star, fmt='o', label="dQ along b*")
plt.errorbar(dQ_data[:, 0], dQ_c_star, yerr=dQ_data[:, 2] / c_star, fmt='o', label="dQ along c*")
plt.xlabel("Energy Transfer (meV)")
plt.ylabel("dQ (RLU)")
plt.legend()
plt.title("Momentum Resolution in Reciprocal Lattice Units")
plt.grid()
plt.show()