import os
import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP
from scipy.interpolate import interp1d
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.optimize import fsolve
from unit_conversions import psia_to_pa, lbm_per_hr_to_kg_per_s, fahrenheit_to_celsius, inches_to_meters, square_feet_to_square_meters


# 1) VOLUMETRIC HEAT GENERATION RATE
def volumetric_heat_generation(z, r, P_nom, n_rods, D, H_active, R_eq, F_q):
    """
    Volumetric heat generation profile along z (chopped cosine).
    D = pellet diameter (m). r, R_eq are only used by the radial Bessel term
    (commented out: the calculation is for the central sub-channel, r=0 -> J0(0)=1).
    Returns: qv_profile (W/m^3), H_e (m), q_avg (W/m^3), q_v_max (W/m^3).
    """
    Tot_power = P_nom * 0.974  # W (power generated inside the fuel)
    q_avg = Tot_power / (n_rods * np.pi * (D/2)**2 * H_active)  # W/m^3
    q_v_max = q_avg * F_q  # W/m^3

    # All quantities here are in metres (the PDF gives them in cm with the same
    # numbers): delta = (D_C/D_R) * L_R is dimensionless * [m] = [m] (cancellation).
    lambda_tr = 0.0029     # m (transport mean free path, 0.29 cm)
    D_c = lambda_tr / 3    # m (diffusion coefficient in the core)
    D_r = 0.16             # m (diffusion coefficient in the reflector)
    L_r = 2.85             # m (diffusion length in the reflector)
    delta = D_c / D_r * L_r  # m (reflector saving ~ 1.72 cm)
    H_e = H_active + 1.42 * lambda_tr + 2 * delta  # m (extrapolated height)

    # qv_profile = q_v_max * np.cos(np.pi * z / H_e) * j0(np.pi * 2.4048*r / R_eq)  # W/m^3
    qv_profile = q_v_max * np.cos(np.pi * z / H_e)
    return qv_profile, H_e, q_avg, q_v_max


# 2) AVERAGE MASS VELOCITY
def average_mass_velocity(m_flow_eff, A_flow_eff):
    """Average mass velocity G = effective mass flow / flow area (kg/m^2 s)."""
    G_avg = m_flow_eff / A_flow_eff
    print(f"Average mass velocity G_avg: {G_avg:.2f} kg/(m^2 s)")
    return G_avg


# 3) COOLANT SPECIFIC ENTHALPY
def coolant_specific_enthalpy(z, G_avg, A_c, q_v_max, H_e, D_pellet, p_sys, H_active, T_in):
    """
    Coolant specific enthalpy profile along z (analytical integral of the
    chopped cosine distribution).
    Returns: h_profile (J/kg), W_hc (kg/s, sub-channel mass flow rate).
    """
    W_hc = G_avg * A_c  # kg/s (sub-channel mass flow rate)
    A_fuel = np.pi / 4 * D_pellet**2  # m^2 (fuel cross-section)

    h_in = CP.PropsSI('H', 'T', T_in + 273.15, 'P', p_sys, 'Water')  # J/kg

    # 1.0267 = 1/0.974: q_v_max is based on 97.4% of the power (energy in the
    # fuel), but the coolant receives the total power (energy transferred to
    # the channel), so the full 100% is "restored" here.
    h_profile = h_in + 1.0267 * (q_v_max * A_fuel * H_e) / (W_hc * np.pi) * (np.sin(np.pi * z / H_e) + np.sin(np.pi * H_active / 2 / H_e))

    return h_profile, W_hc


# 4) TEMPERATURE PROFILE
def temperature_profile(h_profile, p_sys):
    """
    Coolant temperature profile (°C) from h and p, clamped at T_sat.
    Returns: T_profile (°C), first_sat_idx (index of the first saturated node).
    """
    T_profile = CP.PropsSI('T', 'H', h_profile, 'P', p_sys, 'Water') - 273.15  # °C
    # The temperature cannot exceed the saturation value
    T_sat = CP.PropsSI('T', 'P', p_sys, 'Q', 0, 'Water') - 273.15  # °C
    T_profile = np.minimum(T_profile, T_sat)

    sat_indices = np.where(T_profile >= T_sat)[0]
    first_sat_idx = sat_indices[0] if len(sat_indices) > 0 else None

    return T_profile, first_sat_idx


# 5) EQUILIBRIUM QUALITY PROFILE
def equilibrium_quality_profile(h_profile, p_sys):
    """
    Equilibrium quality profile along z: x_eq = (h - h_ls) / (h_vs - h_ls).
    Returns the clamped profile (>= 0), the full profile and H_fg.
    """
    h_ls = CP.PropsSI('H', 'P', p_sys, 'Q', 0, 'Water')  # saturated liquid enthalpy
    h_vs = CP.PropsSI('H', 'P', p_sys, 'Q', 1, 'Water')  # saturated vapour enthalpy
    H_fg = h_vs - h_ls  # latent heat of vaporisation
    x_eq_completo = (h_profile - h_ls) / (h_vs - h_ls)
    x_eq_profile = np.maximum(0, x_eq_completo)  # force negative values to 0
    return x_eq_profile, x_eq_completo, H_fg


# 6) CALCULATION OF THE OUTER CLADDING TEMPERATURE
def safe_props(prop, T_K, P, T_sat_K):
    """Robust CoolProp property: near/above saturation use the saturated liquid."""
    try:
        if T_K >= T_sat_K - 0.01:
            return CP.PropsSI(prop, 'P', P, 'Q', 0, 'Water')
        return CP.PropsSI(prop, 'T', T_K, 'P', P, 'Water')
    except ValueError:
        return CP.PropsSI(prop, 'P', P, 'Q', 0, 'Water')


def T_outer_cladding_profile(T_sat_K, G_avg, D_eq, T_profile, p_sys, C, q_flux):
    """
    Outer cladding temperature: single-phase Dittus-Boelter vs Jens-Lottes
    (boiling), the minimum is taken. Also returns the ONB index.
    """
    # Properties evaluated node by node (T converted to Kelvin)
    mu_profile = np.array([safe_props('V', T + 273.15, p_sys, T_sat_K) for T in T_profile])
    k_profile = np.array([safe_props('L', T + 273.15, p_sys, T_sat_K) for T in T_profile])
    cp_profile = np.array([safe_props('C', T + 273.15, p_sys, T_sat_K) for T in T_profile])

    # Dimensionless numbers
    Re_profile = (G_avg * D_eq) / mu_profile
    Pr_profile = (cp_profile * mu_profile) / k_profile

    # Dittus-Boelter correlation (fluid heating, Pr exponent = 0.4)
    Nu_profile = C * (Re_profile**0.8) * (Pr_profile**0.4)
    h_single_phase = (Nu_profile * k_profile) / D_eq

    # Boiling with Jens-Lottes (q in MW/m^2, p in bar, T in °C)
    q_flux_MW = q_flux * 1e-6
    p_sys_bar = p_sys / 1e5
    T_co_JL = (T_sat_K - 273.15) + 25*(q_flux_MW)**0.25 * np.exp(-p_sys_bar/62)  # °C

    # Single-phase wall temperature
    T_co_SP = T_profile + q_flux / h_single_phase  # °C

    # Combine the two profiles taking the minimum (i.e. phase change occurs)
    T_co = np.minimum(T_co_JL, T_co_SP)

    # ONB index: where Jens-Lottes drops below the single-phase curve
    onb_indices = np.where(T_co_JL < T_co_SP)[0]
    first_onb_idx = onb_indices[0] if len(onb_indices) > 0 else None

    return h_single_phase, T_co, T_co_JL, T_co_SP, first_onb_idx


def detachment(h_l, q_D, T_sat_K, T_cool, zz):
    """
    Bubble detachment point (Griffith model): occurs where T_cool > T_l,
    with T_l = T_sat - q_D/(5*h_l).
    Returns (T_det, z_det, idx) or (None, None, None) if it does not occur.
    """
    T_sat_C = T_sat_K - 273.15
    T_l = T_sat_C - q_D / (5 * h_l)
    detachment_indices = np.where(T_cool > T_l)[0]

    if len(detachment_indices) > 0:
        first_detachment_idx = detachment_indices[0]
        T_det = T_cool[first_detachment_idx]
        z_det = zz[first_detachment_idx]
        return T_det, z_det, first_detachment_idx
    return None, None, None


def flow_quality_two_phase(h_l, T_sat_K, T_cool, index_det, H_fg, p_sys, q_flux, zz, z_det, p_H, W):
    """
    Flow quality x(z) after detachment (Bowring-Rouhani model):
    x = integral of p_H*(q - q_SP) / (H_fg*W*(1+eps)) from z_det to z.
    Computed with Rouhani eps (eps_R) and Bowring eps (eps_B = 1.6).
    """
    # Slice the vectors from the detachment point to the end
    h_l = h_l[index_det:]
    T_cool = T_cool[index_det:]

    # Single-phase heat flux contribution
    T_sat_C = T_sat_K - 273.15  # °C
    q_SP = h_l * (T_sat_C - T_cool)  # W/m^2

    rho_g_sat = CP.PropsSI('D', 'T', T_sat_K, 'Q', 1, 'Water')   # saturated vapour density
    i_l_sat = CP.PropsSI('H', 'T', T_sat_K, 'Q', 0, 'Water')     # J/kg, saturated liquid enthalpy
    rho_l = np.array([safe_props('D', T + 273.15, p_sys, T_sat_K) for T in T_cool])  # liquid density at T_cool
    i_l = np.array([safe_props('H', T + 273.15, p_sys, T_sat_K) for T in T_cool])    # J/kg, liquid enthalpy at T_cool
    eps_R = rho_l / (rho_g_sat*H_fg) * (i_l_sat - i_l)  # ROUHANI
    eps_B = 1.6                                         # BOWRING

    def calc_x_z(eps_arr):
        q_flux_cut = q_flux[index_det:]
        integr = p_H * (q_flux_cut - q_SP) / (H_fg * W * (1 + eps_arr))
        z_int = zz[index_det:]
        x_z_partial = cumulative_trapezoid(integr, z_int, initial=0)

        x_z_total = np.zeros_like(zz)
        x_z_total[index_det:] = x_z_partial
        return x_z_total

    x_flow_R = calc_x_z(eps_R)
    x_flow_B = calc_x_z(eps_B)

    return x_flow_R, x_flow_B


# VOID FRACTION
def void_fraction(p_sys, D_eq, zz, i_ONB, i_Det, G_avg, X_flow):
    """
    Void fraction profile in 3 regions:
      - single phase (z < ONB): alpha = 0
      - ONB -> detachment: linear from 0 to alpha_D (Maurer + Bowring/Rouhani)
      - after detachment: Zuber-Findlay (Collier-Thome) and Fauske/Moody slip ratio
    Returns the 3 profiles: Zuber-Findlay, Fauske, Moody.
    """
    void_fraction_profile = np.zeros_like(zz)
    # Void fraction at detachment from Maurer (alpha = 4*delta/D_h)
    p_bar = p_sys / 1e5  # bar
    R_d = 2.37e-3 / p_bar**0.237              # m (bubble radius, ROUHANI)
    bubble_layer_thickness = 0.0666 * R_d     # m (bubble layer thickness, BOWRING)
    void_fraction_D = 4*bubble_layer_thickness / D_eq  # void fraction at detachment (MAURER)
    # Linear segment between ONB and detachment
    void_fraction_profile[i_ONB:i_Det] = np.linspace(0, void_fraction_D, i_Det - i_ONB)

    # Post-detachment, model 1: ZUBER-FINDLAY (Collier-Thome, valid in this p range)
    C_0 = 1.13
    C_1 = 0.5 * (1.18 + 1.4)
    rho_g_sat = CP.PropsSI('D', 'P', p_sys, 'Q', 1, 'Water')
    rho_l_sat = CP.PropsSI('D', 'P', p_sys, 'Q', 0, 'Water')
    sigma = CP.PropsSI('surface_tension', 'P', p_sys, 'Q', 0, 'Water')
    g = 9.81

    # Vapour quality from the detachment point onwards
    x = X_flow[i_Det:]

    drift_term = (sigma * g * (rho_l_sat - rho_g_sat) / rho_l_sat**2)**0.25
    alpha_post_det = (x / rho_g_sat) / (C_0 * (x / rho_g_sat + (1 - x) / rho_l_sat) + C_1 * (1 / G_avg) * drift_term)

    # NB: at z = z_det the quality x = 0, so Z-F would give alpha = 0, creating
    # a discontinuity with the linear segment ending at alpha_D. To ensure
    # continuity alpha_D is added ("wall bubble layer" + flow void).
    # Conservative approximation: slightly overestimates alpha in the bulk.
    alpha_post_det = void_fraction_D + alpha_post_det

    void_fraction_profile_ZF = void_fraction_profile.copy()
    void_fraction_profile_ZF[i_Det:] = alpha_post_det

    # Post-detachment, model 2: SLIP RATIO
    S_F = (rho_l_sat / rho_g_sat)**0.5     # FAUSKE
    S_M = (rho_l_sat / rho_g_sat)**(1/3)   # MOODY

    # Matched by adding void_fraction_D (as for Z-F)
    alpha_F = void_fraction_D + x/(x + (1-x)*rho_g_sat/rho_l_sat * S_F)
    alpha_M = void_fraction_D + x/(x + (1-x)*rho_g_sat/rho_l_sat * S_M)

    void_fraction_profile_F = void_fraction_profile.copy()
    void_fraction_profile_M = void_fraction_profile.copy()
    void_fraction_profile_F[i_Det:] = alpha_F
    void_fraction_profile_M[i_Det:] = alpha_M

    return void_fraction_profile_ZF, void_fraction_profile_F, void_fraction_profile_M


# 7) T INNER CLADDING TEMPERATURE PROFILE
def T_inner_cladding_profile(T_co, q_vol, A_f, D_ci, D_co):
    """
    Inner cladding temperature from hollow-cylinder conduction with k_Zr = A + B*T:
    integrating, A*(T_ci-T_co) + B/2*(T_ci^2-T_co^2) = q_v*A_f/(2pi)*ln(D_co/D_ci).
    Quadratic equation in T_ci, the positive physical root is taken.
    """
    T_inner = np.zeros_like(T_co)
    A = 11.45      # k_Zr = A + B*T  (W/m°C), WCAP 3269-4-1
    B = 1.425e-2

    for i in range(len(T_co)):
        T_co_loc = T_co[i]
        rhs = q_vol[i] * A_f / (2 * np.pi) * np.log(D_co / D_ci)

        # (B/2)*T_ci^2 + A*T_ci + c = 0
        a = B / 2
        b = A
        c = -rhs - (A * T_co_loc + B/2 * T_co_loc**2)
        delta = b**2 - 4*a*c
        if delta < 0:
            raise ValueError("Negative discriminant, no real solution for T_inner")
        T_inner[i] = (-b + np.sqrt(delta)) / (2*a)

    return T_inner


# 8) T PELLET SURFACE TEMPERATURE PROFILE
def pellet_thermal_conductivity(T):
    """UO2 pellet thermal conductivity (Westinghouse). T in °C -> W/(m·K)."""
    k_UO2 = 1/(11.8+0.0238*T) + 8.775e-13*T**3  # W/(cm·K)
    return k_UO2 * 100  # W/(m·K)


# --- Gap conductance: thermal expansion, elastic deformation, h_g_T ---
def thermal_expansion(D_Ta, T_mean, T_amb, component):
    """Diameter change due to thermal expansion (fuel or cladding)."""
    def alpha(T):  # T in °C, returns 1/°C
        if component == 'fuel':
            return 7.87e-6 + 3.9e-9*T
        elif component == 'cladding':
            return 5.62e-6 + 3.162e-9*T
        else:
            raise ValueError("Unknown component")
    return D_Ta*alpha(T_mean)*(T_mean - T_amb)


def elastic_deformation(r_ci, r_co, T_c, p_i, p_e):
    """Elastic radial deformation of the cladding due to the pressure difference."""
    gamma = r_co/r_ci
    nu = 0.43
    def young_modulus(T):
        return 1.148e11 - 5.99e7 * T  # Pa (Young's modulus vs temperature)
    E = young_modulus(T_c + 273.15)
    delta_r_ci = (r_ci / (E * (gamma**2 - 1))) * (p_i * ((1 - nu) + (1 + nu) * gamma**2) - 2 * gamma**2 * p_e)
    return delta_r_ci


def calculate_gap_conductance(delta0, D_pellet, D_in_clad, D_out_clad, T_amb, T_c_avg, p_i, p_sys, T_fuel_avg, T_f_S, T_ci):
    """Total gap conductance accounting for expansion, deformation and radiation."""
    # 1. Expansion and deformation (they change the gap thickness)
    delta_f = thermal_expansion(D_pellet, T_fuel_avg, T_amb, 'fuel')
    delta_cl = thermal_expansion(D_in_clad, T_c_avg, T_amb, 'cladding')
    delta_def = elastic_deformation(D_in_clad/2, D_out_clad/2, T_c_avg, p_i, p_sys)

    # Final gap: delta_f shrinks it, delta_cl widens it, delta_def shrinks/widens it
    delta = delta0 - delta_f/2 + delta_cl/2 + delta_def

    if np.any(delta <= 0):
        print("Warning: pellet-cladding contact detected!")

    # 2. Gas (He) conductivity at the mean of pellet surface and inner cladding
    T_gas_avg_K = (T_f_S + T_ci) / 2 + 273.15
    k_gas = 0.1763e-2 * T_gas_avg_K**0.77163  # W/mK

    # Gas conductive coefficient with jump distance
    h_gap_gas = k_gas / (delta + 2.54e-5)

    # 3. Radiative coefficient (emissivity = 1 -> upper bound for h_rad)
    sigma = 5.67e-8
    T_f_S_K = T_f_S + 273.15
    T_ci_K = T_ci + 273.15
    h_rad_val = sigma * (T_f_S_K**2 + T_ci_K**2) * (T_f_S_K + T_ci_K)

    return h_gap_gas + h_rad_val, delta


def calculate_T_pellet_surface_iterative(T_ci, T_co, q_vol_profile, A_fuel, D_pellet, D_in_clad, D_out_clad, T_amb, p_sys, delta0):
    """
    Pellet surface temperature: T_f_S = T_ci + q''_f / h_g_T.
    Iterative because h_g_T (gap conductance) and q''_f depend on the pellet
    temperatures, which are initially unknown.
    Returns: T_f_S, h_tot, delta (gap), T_center (estimate), n. of iterations.
    """
    T_f_S = T_ci.copy() + 100.0  # initial guess
    T_fuel_avg = T_f_S + 200.0   # initial guess for the pellet mean temperature
    T_c_avg = (T_ci + T_co) / 2

    p_i = 7e6  # Pa (fill-gas internal pressure)
    tol = 1e-3
    max_iter = 1000

    for iteration in range(max_iter):
        T_f_S_old = T_f_S.copy()

        # 1. Pellet properties from the previous-step temperatures
        k_pellet = pellet_thermal_conductivity(T_fuel_avg)
        T_center = T_f_S_old + (q_vol_profile * A_fuel) / (4 * np.pi * k_pellet)
        T_fuel_avg = T_f_S_old + (q_vol_profile * A_fuel) / (8 * np.pi * k_pellet)

        # 2. Gap conductance and effective gap
        h_tot, delta = calculate_gap_conductance(
            delta0, D_pellet, D_in_clad, D_out_clad,
            T_amb, T_c_avg, p_i, p_sys, T_fuel_avg, T_f_S_old, T_ci
        )

        # 3. New surface temperature imposed by the flux: q'' = q_vol * A_fuel / (pi * D_pellet)
        heat_flux_pellet = (q_vol_profile * A_fuel) / (np.pi * D_pellet)
        T_f_S_new = T_ci + heat_flux_pellet / h_tot

        # 4. Convergence check
        error = np.linalg.norm(T_f_S_new - T_f_S_old) / np.linalg.norm(T_f_S_old)
        T_f_S = T_f_S_new
        if error < tol:
            print(f"Convergence reached after {iteration} iterations.")
            break
    else:
        print("WARNING: max_iter reached, the pellet surface temperature is not fully converged!")

    if np.any(delta <= 0):
        print("WARNING: pellet-cladding contact detected at convergence!")

    return T_f_S, h_tot, delta, T_center, iteration


# 9) FUEL CENTER LINE TEMPERATURE
def calculate_fuel_centerline_temperature(T_f_S_profile, q_vol_profile, A_fuel):
    """
    Fuel centerline temperature, solved node by node from the conductivity
    integral equation: int_{T_S}^{T_CL} k_UO2 dT = q' * f_robertson / (4 pi),
    where f_robertson is the centerline flux depression factor.
    """
    # Analytical integral of k_UO2 (Westinghouse): k = 100*[1/(11.8+0.0238T) + 8.775e-13 T^3]
    def integral_k(T):
        return (100 / 0.0238) * np.log(11.8 + 0.0238 * T) + 100 * (8.775e-13 / 4) * (T**4)

    f_robertson = 0.965
    T_centerline = np.zeros_like(T_f_S_profile)

    for i in range(len(T_f_S_profile)):
        T_surface_loc = T_f_S_profile[i]
        RHS = (q_vol_profile[i] * A_fuel * f_robertson) / (4 * np.pi)

        # Root of: integral_k(T_CL) - integral_k(T_surface) - RHS
        def objective_function(T_CL_guess):
            return integral_k(T_CL_guess) - integral_k(T_surface_loc) - RHS

        T_centerline[i] = fsolve(objective_function, T_surface_loc + 500.0)[0]

    return T_centerline


def solve_pellet_radial_temperature(z, q_vol_profile, T_surface, R_pellet, pellet_thermal_conductivity, Nr=100, tol=1e-6, max_iter=500, alpha=0.5):
    """
    Solves radial conduction in the pellet at each axial location z.

    Equation:
        1/r d/dr ( r k(T) dT/dr ) + q'''(z) = 0
    BC:
        dT/dr = 0          at r = 0
        T = T_surface(z)   at r = R_pellet
    """
    f_rob = 0.965  # Robertson factor
    q_vol_profile = q_vol_profile * f_rob
    r = np.linspace(0, R_pellet, Nr)
    dr = r[1] - r[0]

    Nz = len(z)
    T = np.zeros((Nr, Nz))
    T_center = np.zeros(Nz)

    for iz in range(Nz):
        q = q_vol_profile[iz]
        T_s = T_surface[iz]
        T_old = np.ones(Nr) * T_s  # initial guess

        for iteration in range(max_iter):
            k = pellet_thermal_conductivity(T_old)

            A = np.zeros((Nr, Nr))
            b = np.zeros(Nr)

            # Centre BC: dT/dr = 0  -> T[0] = T[1]
            A[0, 0] = 1.0
            A[0, 1] = -1.0
            b[0] = 0.0

            # Internal nodes
            for i in range(1, Nr - 1):
                rp = r[i] + dr / 2
                rm = r[i] - dr / 2
                kp = 0.5 * (k[i] + k[i + 1])
                km = 0.5 * (k[i] + k[i - 1])

                A[i, i - 1] = rm * km / dr**2
                A[i, i] = -(rp * kp + rm * km) / dr**2
                A[i, i + 1] = rp * kp / dr**2
                b[i] = -q * r[i]

            # Surface BC: T(R) = T_s
            A[-1, -1] = 1.0
            b[-1] = T_s

            T_calc = np.linalg.solve(A, b)

            # Under-relaxation
            T_new = alpha * T_calc + (1 - alpha) * T_old
            error = np.max(np.abs(T_new - T_old))
            if error < tol:
                break
            T_old = T_new.copy()

        T[:, iz] = T_new
        T_center[iz] = T_new[0]

    return r, T, T_center


# 10) CRITICAL FLUX BY W3 AND GRID FACTOR
def critical_flux_uniform(G_avg, D_eq, p_sys, h_in, h_lsat, x_c, L_active_in):
    """
    Uniform critical heat flux: W-3 correlation (15x15) with the 0.88 correction
    for the 17x17 bundle, multiplied by the grid spacer factor Fs.
    Returns: qc_w3_17 (corrected W-3), qc_eu (W-3 + grid factor), in kW/m^2.
    """
    # --- W-3 (units: p in MPa, G in kg/m^2 s, i in kJ/kg, D_h in m) ---
    p_mpa = p_sys / 1e6
    t1 = (2.022 - 0.06238*p_mpa)
    t2 = (0.1722 - 0.01427*p_mpa) * np.exp((18.177 - 0.5987*p_mpa)*x_c)
    t3 = (0.1484 - 1.596*x_c + 0.1729*x_c*np.abs(x_c))*(2.326*G_avg) + 3271
    t4 = (1.157 - 0.869*x_c)
    t5 = (0.2664 + 0.8357*np.exp(-124.1*D_eq))
    t6 = (0.8285 + 0.0003413*(h_lsat - h_in))
    qc_w3_15 = (t1+t2)*t3*t4*t5*t6      # kW/m^2
    qc_w3_17 = qc_w3_15 * 0.88          # 15x15 -> 17x17 correction

    # --- Grid spacer factor Fs (units: p in psi, L in ft, G in lb/ft^2 hr) ---
    p_psi = p_sys * 0.000145038
    G_lb = G_avg * 737.338
    L_ft = L_active_in / 12
    alpha = 0.038

    # Ks from the table as a function of L_s = active length / n. of grids.
    # L_s = 16.8 in falls below the table range [20, 32] in -> linear extrapolation.
    N_grids = 10
    L_s = L_active_in / N_grids  # inches
    L_s_vect = np.array([20, 26, 32])      # inches (increasing)
    Ks = np.array([0.066, 0.046, 0.027])
    Ks_interp = interp1d(L_s_vect, Ks, kind='linear', fill_value='extrapolate')(L_s)

    Fs = (p_psi/225.896)**0.5*(1.445-0.0371*L_ft)*(np.exp((x_c+0.2)**2)-0.73) + Ks_interp*G_lb/1e6*(alpha/0.019)**0.35
    qc_eu = qc_w3_17 * Fs
    return qc_w3_17, qc_eu


def critical_flux_non_uniform(qc_eu, G_avg, z_array, q_flux_array, z_ONB_idx, x_c):
    """
    Non-uniform critical heat flux: Tong F-factor (Lin 1991 form).
    F(z) = C * integral[ONB->z] q(z') exp(-C(z-z')) dz' / (q(z)*(1-exp(-C(z-z_ONB)))).
    Upstream of ONB F = 1. Returns: F_profile, qc_NU = qc_eu / F.
    """
    F_profile = np.ones_like(z_array)
    G_avg_Mlb = G_avg * 737.338e-6  # kg/(m^2 s) -> Mlb/(ft^2 hr)
    z_ONB = z_array[z_ONB_idx]

    # F is evaluated only at the nodes downstream of ONB
    for i in range(z_ONB_idx + 1, len(z_array)):
        z_current = z_array[i]
        q_current = q_flux_array[i]

        # Tong-Lin C: C in [in^-1] with G in Mlb/hr/ft^2 -> converted to [m^-1]
        C_in = 0.15 * (1 - x_c[i])**4.31 / G_avg_Mlb**0.478
        C = C_in / 0.0254

        # Dummy vectors z' from ONB to the current node; trapezoidal integral
        z_prime = z_array[z_ONB_idx : i+1]
        q_prime = q_flux_array[z_ONB_idx : i+1]
        integrand = q_prime * np.exp(-C * (z_current - z_prime))
        integral_value = trapezoid(integrand, z_prime)

        denominator = q_current * (1.0 - np.exp(-C * (z_current - z_ONB)))
        F_profile[i] = (C / denominator) * integral_value

    qc_NU = qc_eu / F_profile  # kW/m^2
    return F_profile, qc_NU


# 11) DNBR and MINIMUM DNBR
def DNBR_calculation(q_flux, qc_nu):
    """
    DNBR(z) = q_c,NU(z) / q(z).
    IMPORTANT: q_flux and qc_nu must have the SAME unit (both kW/m^2 or both
    W/m^2). In the main q_flux/1e3 is used because qc_nu is in kW/m^2 while
    q_flux elsewhere is in W/m^2.
    """
    DNBR = qc_nu / q_flux  # dimensionless
    MDNBR = np.min(DNBR)   # Minimum Departure from Nucleate Boiling Ratio
    return DNBR, MDNBR


# ============================================================================
# MAIN
# ============================================================================
if __name__ == "__main__":

    # ------------------- INPUT: operating conditions -------------------
    P_nom = 3400e6  # W
    p_sys = psia_to_pa(2250)  # Pa
    n_rods = 157 * 264
    m_flow_tot = lbm_per_hr_to_kg_per_s(113.5e6)  # kg/s
    F_q = 2.6
    Bypass_fraction = 0.059
    m_flow_eff = m_flow_tot * (1 - Bypass_fraction)
    A_flow_eff = square_feet_to_square_meters(41.8)  # m^2
    T_sat = CP.PropsSI('T', 'P', p_sys, 'Q', 0, 'Water')  # K
    T_in = fahrenheit_to_celsius(535)  # °C

    # ------------------- INPUT: geometry (square array) ----------------
    D_out_clad = inches_to_meters(0.374)
    H_active_in = 168
    H_active = inches_to_meters(H_active_in)
    w = inches_to_meters(0.496)  # pitch
    s_clad = inches_to_meters(0.0225)
    D_in_clad = D_out_clad - 2 * s_clad
    D_pellet = inches_to_meters(0.3225)
    R_eq = inches_to_meters(119.7)

    # ------------------- Derived geometry ------------------------------
    P_wet = np.pi * D_out_clad
    A_c = w**2 - np.pi/4 * D_out_clad**2
    D_eq = 4 * A_c / P_wet
    A_fuel = np.pi/4 * D_pellet**2
    R_pellet = D_pellet / 2
    C = 0.042*w/D_out_clad - 0.024  # Dittus-Boelter coefficient

    # ------------------- Axial and radial grid -------------------------
    dzz = 0.0254/2  # half an inch
    z = np.arange(-H_active/2, H_active/2 + dzz, dzz)
    r = np.arange(0, R_eq + dzz, dzz)

    # ===================================================================
    # 1-2) Volumetric heat generation and heat flux on the cladding
    # ===================================================================
    qv_profile, H_e, q_avg, q_v_max = volumetric_heat_generation(
        z, r, P_nom, n_rods, D_pellet, H_active, R_eq, F_q)
    q_flux = qv_profile * A_fuel / P_wet  # W/m^2
    
    # ===================================================================
    # 3) Average mass velocity
    # ===================================================================
    G_avg = average_mass_velocity(m_flow_eff, A_flow_eff)

    # ===================================================================
    # 4) Coolant enthalpy and temperature
    # ===================================================================
    h_profile, W_hc = coolant_specific_enthalpy(
        z, G_avg, A_c, q_v_max, H_e, D_pellet, p_sys, H_active, T_in)
    T_profile, first_sat_idx = temperature_profile(h_profile, p_sys)
    z_sat = z[first_sat_idx] if first_sat_idx is not None else None

    # ===================================================================
    # 5) Equilibrium quality
    # ===================================================================
    x_eq_profile, x_eq_completo, H_fg = equilibrium_quality_profile(h_profile, p_sys)

    # ===================================================================
    # 6) Outer cladding, flow quality, void fraction
    # ===================================================================
    h_single_phase, T_co, T_co_JL, T_co_SP, first_onb_idx = T_outer_cladding_profile(
        T_sat, G_avg, D_eq, T_profile, p_sys, C, q_flux)
    z_NB = z[first_onb_idx] if first_onb_idx is not None else None        # ONB position
    T_co_NB = T_co[first_onb_idx] if first_onb_idx is not None else None  # cladding T at ONB

    T_det, z_det, first_detachment_idx = detachment(
        h_single_phase, q_flux, T_sat, T_profile, z)

    x_flow_R, x_flow_B = flow_quality_two_phase(
        h_single_phase, T_sat, T_profile, first_detachment_idx,
        H_fg, p_sys, q_flux, z, z_det, P_wet, W_hc)

    alpha_ZF_R, alpha_F_R, alpha_M_R = void_fraction(
        p_sys, D_eq, z, first_onb_idx, first_detachment_idx, G_avg, x_flow_R)
    alpha_ZF_B, alpha_F_B, alpha_M_B = void_fraction(
        p_sys, D_eq, z, first_onb_idx, first_detachment_idx, G_avg, x_flow_B)

    # ===================================================================
    # 7-8-9) Temperature profiles from inner cladding to pellet centerline
    # ===================================================================
    T_ci = T_inner_cladding_profile(T_co, qv_profile, A_fuel, D_in_clad, D_out_clad)

    # ------------------- Export per Assignment 3 -----------------------
    # Salva i profili di temperatura della guaina (parete interna/esterna) in un
    # CSV, cosi' l'Assignment 3 puo' rileggerli senza ricalcolare la catena TH.
    out_csv = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'output_cladding_T.csv')
    np.savetxt(out_csv, np.column_stack([z, T_ci, T_co]),
               delimiter=',', header='z_m,T_ci_C,T_co_C', comments='')
    print(f"Cladding temperatures exported to: {out_csv}")

    T_amb = 25
    delta0 = 0.5*(D_in_clad - D_pellet)
    T_f_S, h_tot, delta_out, T_center, iteration = calculate_T_pellet_surface_iterative(
        T_ci, T_co, qv_profile, A_fuel, D_pellet, D_in_clad, D_out_clad,
        T_amb, p_sys, delta0)

    r, T_radial, T_center = solve_pellet_radial_temperature(
        z=z, q_vol_profile=qv_profile, T_surface=T_f_S, R_pellet=R_pellet,
        pellet_thermal_conductivity=pellet_thermal_conductivity,
        Nr=100, tol=1e-6, max_iter=500, alpha=0.5)

    T_centerline = calculate_fuel_centerline_temperature(T_f_S, qv_profile, A_fuel)

    # ===================================================================
    # 10) Critical flux W3, grid factor, Tong F-factor, MDNBR
    # ===================================================================
    h_in_kJ = h_profile[0] / 1e3
    h_lsat_kJ = CP.PropsSI('H', 'P', p_sys, 'Q', 0, 'Water') / 1e3
    qc_w3, qc_eu = critical_flux_uniform(
        G_avg, D_eq, p_sys, h_in_kJ, h_lsat_kJ, x_eq_completo, H_active_in)
    F_profile, qc_NU = critical_flux_non_uniform(
        qc_eu, G_avg, z, q_flux, first_onb_idx, x_eq_completo)
    DNBR, MDNBR = DNBR_calculation(q_flux/1e3, qc_NU)  # both in kW/m^2
    MDNBR_idx = np.argmin(DNBR)
    z_MDNBR = z[MDNBR_idx]


    # ===================================================================
    # PLOT
    # ===================================================================
    plot_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'plots')
    os.makedirs(plot_dir, exist_ok=True)

    def plot_axial(filename, xlabel, title, curves, *, ylim=None, markers=(),
                   hlines=(), figsize=(10, 6)):
        """Save an x-vs-z plot with one or more curves. Each curve is a dict
        with key 'x' (data) + optional plt.plot kwargs (label, linestyle, ...)."""
        plt.figure(figsize=figsize)
        show_legend = False
        for c in curves:
            kwargs = {k: v for k, v in c.items() if k != 'x'}
            plt.plot(c['x'], z, **kwargs)
            show_legend = show_legend or 'label' in kwargs
        for mk in markers:
            plt.plot(mk['x'], mk['z'], mk.get('fmt', 'ro'), label=mk.get('label'))
            show_legend = True
        for hl in hlines:
            plt.axhline(hl['z'], color=hl.get('color', 'red'),
                        linestyle=hl.get('linestyle', ':'), alpha=hl.get('alpha', 0.4))
        plt.title(title)
        plt.ylabel('z (m)')
        plt.xlabel(xlabel)
        if ylim is not None:
            plt.ylim(ylim)
        if show_legend:
            plt.legend()
        plt.grid()
        plt.savefig(os.path.join(plot_dir, filename))
        plt.close()

    # 1) Volumetric heat generation
    plot_axial('1_volumetric_heat_generation.png', 'q_v (W/m^3)',
               'Volumetric Heat Generation Rate along the z-axis',
               [{'x': qv_profile}])

    # 4.1) Coolant enthalpy
    plot_axial('4.1_coolant_specific_enthalpy.png', 'h (J/kg)',
               'Coolant Specific Enthalpy Profile along the z-axis',
               [{'x': h_profile}])

    # 4.2) Coolant temperature with saturation point
    T_markers = []
    if z_sat is not None:
        T_markers.append({'x': T_sat - 273.15, 'z': z_sat, 'fmt': 'go',
                          'label': f'Saturated Point (T={T_sat-273.15:.1f} °C, z={z_sat:.2f} m)'})
    plot_axial('4.2_coolant_temperature_profile.png', 'T (°C)',
               'Coolant Temperature Profile along the z-axis',
               [{'x': T_profile, 'label': 'T_coolant (°C)'}], markers=T_markers)

    # 5) Equilibrium quality
    plot_axial('5_equilibrium_quality.png', 'x (kg/kg)',
               'Equilibrium Quality Profile along the z-axis',
               [{'x': x_eq_completo}])

    # 6.1) Outer cladding temperature, with the ONB point
    onb_markers = []
    if T_co_NB is not None:
        onb_markers.append({'x': T_co_NB, 'z': z_NB, 'fmt': 'bo',
                            'label': f'ONB Point (T={T_co_NB:.1f} °C, z={z_NB:.2f} m)'})
    plot_axial('6.1_outer_cladding_temperature.png', 'Temperature (°C)',
               'Outer Cladding Temperature Profile along the z-axis',
               [{'x': T_co, 'label': 'T_co (Actual)', 'color': 'black', 'linewidth': 2},
                {'x': T_co_JL, 'label': 'T_co_JL (Jens-Lottes)', 'linestyle': '--'},
                {'x': T_co_SP, 'label': 'T_co_SP (Single Phase)', 'linestyle': '-.'}],
               markers=onb_markers)

    # 6.2) Flow quality
    plot_axial('6.2_flow_quality.png', 'x (kg/kg)',
               'Flow Quality Profile along the z-axis',
               [{'x': x_flow_R, 'label': 'Rouhani eps'},
                {'x': x_flow_B, 'label': 'Bowring eps=1.6', 'linestyle': '--'}])

    # 6.3) Void fraction
    plot_axial('6.3_void_fraction.png', 'Void Fraction',
               'Void Fraction Profile along the z-axis',
               [{'x': alpha_ZF_R, 'label': 'Void Fraction (Zuber-Findlay + Rouhani)', 'color': 'blue'},
                {'x': alpha_F_R,  'label': 'Void Fraction (Fauske + Rouhani)', 'color': 'blue', 'linestyle': '--'},
                {'x': alpha_M_R,  'label': 'Void Fraction (Moody + Rouhani)',  'color': 'blue', 'linestyle': '-.'},
                {'x': alpha_ZF_B, 'label': 'Void Fraction (Zuber-Findlay + Bowring)', 'color': 'green'},
                {'x': alpha_F_B,  'label': 'Void Fraction (Fauske + Bowring)', 'color': 'green', 'linestyle': '--'},
                {'x': alpha_M_B,  'label': 'Void Fraction (Moody + Bowring)',  'color': 'green', 'linestyle': '-.'}],
               ylim=(-0.25, max(z)))

    # 7) Inner cladding temperature
    plot_axial('7_inner_cladding_temperature.png', 'T_inner (°C)',
               'Inner Cladding Temperature Profile along the z-axis',
               [{'x': T_ci, 'label': 'T_inner_cladding (°C)', 'color': 'red'},
                {'x': T_co, 'label': 'T_outer_cladding (°C)', 'color': 'black', 'linestyle': '--'}])

    # 8) Fuel surface temperature
    plot_axial('8_fuel_surface_temperature.png', 'T_surface_fuel (°C)',
               'Pellet Surface Temperature Profile along the z-axis',
               [{'x': T_f_S, 'label': 'T_surface_fuel (°C)', 'color': 'red'}])

    # 8bis) Gap thickness with its minimum
    min_idx = np.argmin(delta_out)
    plot_axial('8_gap_thickness.png', 'Gap Thickness (m)',
               'Gap Thickness Profile along the z-axis',
               [{'x': delta_out, 'label': 'Gap Thickness (m)', 'color': 'red'}],
               markers=[{'x': delta_out[min_idx], 'z': z[min_idx], 'fmt': 'bo',
                         'label': f'Min gap: {delta_out[min_idx]:.2e} m (z = {z[min_idx]:.2f} m)'}])

    # 9) Fuel centerline (analytical)
    plot_axial('9_fuel_centerline_temperature.png', 'T_centerline_fuel (°C)',
               'Fuel Center Line Temperature Profile along the z-axis',
               [{'x': T_centerline, 'label': 'T_centerline_fuel (°C)', 'color': 'red'}])

    # 9bis) Fuel center from the radial solution
    plot_axial('9bis_center_fuel_temp.png', 'T center [K]',
               'Pellet Center Temperature along the z-axis',
               [{'x': T_center, 'label': 'T_centerline_fuel (°C)', 'color': 'red'}])

    # 9 all) All temperatures
    T_fuel_max_idx = int(np.argmax(T_centerline))
    plot_axial('9_all_temperatures.png', 'Temperature (°C)',
               'All Temperatures Profile along the z-axis',
               [{'x': T_profile, 'label': 'Coolant (T_profile)', 'color': 'blue'},
                {'x': T_co, 'label': 'Cladding Outside (T_co)', 'color': 'black'},
                {'x': T_ci, 'label': 'Cladding Inside (T_ci)', 'color': 'green'},
                {'x': T_f_S, 'label': 'Fuel Surface (T_f_S)', 'color': 'orange'},
                {'x': T_centerline, 'label': 'Fuel Centerline', 'color': 'red'}],
               markers=[{'x': T_centerline[T_fuel_max_idx], 'z': z[T_fuel_max_idx], 'fmt': 'ro',
                         'label': f'T_fuel max = {T_centerline[T_fuel_max_idx]:.1f} °C (z = {z[T_fuel_max_idx]:.2f} m)'}],
               figsize=(10, 8))

    # 11) Critical heat flux + MDNBR
    plot_axial('11_critical_flux_uniform.png', 'q (kW/m^2)',
               'Critical HeatFlux Profile along the z-axis',
               [{'x': qc_eu, 'label': 'q_c uniform (W3 + grid factor)'},
                {'x': qc_NU, 'label': 'q_c non-uniform (Tong F-factor)', 'linestyle': '--'},
                {'x': q_flux/1e3, 'label': 'q_flux (actual)', 'linestyle': ':'}],
               markers=[{'x': q_flux[MDNBR_idx]/1e3, 'z': z_MDNBR, 'fmt': 'ro',
                         'label': f'MDNBR = {MDNBR:.2f} (z = {z_MDNBR:.2f} m)'}],
               hlines=[{'z': z_MDNBR}])

    # ===================================================================
    # 12) Design thermal limits verification
    # ===================================================================
    # Typical PWR limits (assignment slides + Westinghouse manuals):
    #   MDNBR >= 1.85 at nominal power
    #   MDNBR >= 1.30 at maximum overpower (typically 115% of P_nom)
    #   T_fuel,CL < 2800 degC (safety limit, UO2 melting ~2840 degC)
    #   T_clad,out < 350 degC at nominal operation (design limit)
    MDNBR_LIMIT_NOMINAL   = 1.85
    MDNBR_LIMIT_OVERPOWER = 1.30
    OVERPOWER_FACTOR      = 1.15
    T_FUEL_LIMIT_C        = 2800.0
    T_CLAD_LIMIT_C        = 350.0

    MDNBR_at_overpower = MDNBR / OVERPOWER_FACTOR  # estimate: q rises, q_c ~constant
    T_fuel_max = float(np.max(T_centerline))
    T_clad_max = float(np.max(T_co))

    def _status(value, limit, higher_is_safer):
        ok = (value >= limit) if higher_is_safer else (value <= limit)
        return "OK" if ok else "FAIL"

    print("\n" + "=" * 70)
    print(" DESIGN THERMAL LIMITS VERIFICATION")
    print("=" * 70)
    print(f" MDNBR nominal           = {MDNBR:6.3f}   "
          f"(limit >= {MDNBR_LIMIT_NOMINAL:.2f})   "
          f"[{_status(MDNBR, MDNBR_LIMIT_NOMINAL, higher_is_safer=True)}]")
    print(f"   z(MDNBR) location     = {z_MDNBR:+.3f} m")
    print(f" MDNBR @ {OVERPOWER_FACTOR*100:.0f}% overpower  = {MDNBR_at_overpower:6.3f}   "
          f"(limit >= {MDNBR_LIMIT_OVERPOWER:.2f})   "
          f"[{_status(MDNBR_at_overpower, MDNBR_LIMIT_OVERPOWER, higher_is_safer=True)}]")
    print(f" T fuel centerline max   = {T_fuel_max:6.1f} degC  "
          f"(limit <= {T_FUEL_LIMIT_C:.0f} degC) "
          f"[{_status(T_fuel_max, T_FUEL_LIMIT_C, higher_is_safer=False)}]")
    print(f" T cladding outer max    = {T_clad_max:6.1f} degC  "
          f"(limit <= {T_CLAD_LIMIT_C:.0f} degC)  "
          f"[{_status(T_clad_max, T_CLAD_LIMIT_C, higher_is_safer=False)}]")
    print("=" * 70)

