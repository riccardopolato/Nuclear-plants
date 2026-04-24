import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP
import numpy as np
import scipy
from scipy.integrate import cumulative_trapezoid


# 1) VOLUMETRIC HEAT GENERATION RATE
def volumetric_heat_generation(z, P_nom, n_rods, D_in_clad, H_active, F_q):
    """
    Calcola il profilo di generazione di calore volumetrico lungo l'asse z.
    
    Parameters:
    - z: array delle posizioni assiali (m)
    - P_nom: potenza nominale (W)
    - n_rods: numero di barre combustibili
    - D_in_clad: diametro interno della guaina (m)
    - H_active: altezza attiva del core (m)
    - F_q: fattore di picco del flusso termico
    
    Returns:
    - qv_profile: profilo di generazione di calore volumetrico (W/m^3)
    """
    Tot_power = P_nom * 0.974  # W (potenza totale generata)
    q_avg = Tot_power / (n_rods * np.pi * (D_in_clad/2)**2 * H_active)  # W/m^3
    q_v_max = q_avg * F_q  # W/m^3 (tasso massimo)
    
    lambda_tr = 0.0029  # m (lunghezza di trasporto)
    D_c = lambda_tr / 3  # diffusion coeff. in the core
    D_r = 0.16  # diffusion coeff. in the reflector
    L_r = 2.85  # m (lunghezza del riflettore)
    delta = D_c / D_r * L_r  # transport length
    H_e = H_active + 1.42 * lambda_tr + 2 * delta  # m altezza estrapolata
    
    qv_profile = q_v_max * np.cos(np.pi * z / H_e)  # W/m^3
    
    return qv_profile


# 2) AVERAGE MASS VELOCITY
def average_mass_velocity(m_flow_eff, A_flow_eff):
    """
    Calcola la velocità di massa media.
    
    Parameters:
    - m_flow_eff: portata effettiva (kg/s)
    - A_flow_eff: area di flusso effettiva (m^2)
    
    Returns:
    - G_avg: velocità di massa media (kg/(m^2 s))
    """
    G_avg = m_flow_eff / A_flow_eff  # kg/(m^2 s)
    print(f"Average mass velocity G_avg: {G_avg:.2f} kg/(m^2 s)")
    
    return G_avg


# 3) COOLANT SPECIFIC ENTHALPY
def coolant_specific_enthalpy(z, G_avg, A_c, w, Dout_clad, D_pellet, p_sys, 
                               P_nom, n_rods, D_in_clad, H_active, F_q):
    """
    Calcola il profilo di entalpia specifica del refrigerante lungo l'asse z.
    
    Parameters:
    - z: array delle posizioni assiali (m)
    - G_avg: velocità di massa media (kg/(m^2 s))
    - A_c: area di flusso per canale (m^2)
    - w: pitch tra le barre (m)
    - Dout_clad: diametro esterno della guaina (m)
    - D_pellet: diametro del pellet (m)
    - p_sys: pressione del sistema (Pa)
    - P_nom: potenza nominale (W)
    - n_rods: numero di barre combustibili
    - D_in_clad: diametro interno della guaina (m)
    - H_active: altezza attiva del core (m)
    - F_q: fattore di picco del flusso termico
    
    Returns:
    - h_profile: profilo di entalpia specifica (J/kg)
    """
    W_hc = G_avg * A_c  # kg/s (portata per canale)
    A_fuel = np.pi / 4 * D_pellet**2  # m^2 (area del combustibile)
    
    T_in = 279.44 + 273.15  # K (temperatura di ingresso)
    h_in = CP.PropsSI('H', 'T', T_in, 'P', p_sys, 'Water')  # J/kg
    
    # Calcolo dei parametri necessari per il profilo di generazione di calore
    Tot_power = P_nom * 0.974
    q_avg = Tot_power / (n_rods * np.pi * (D_in_clad/2)**2 * H_active)
    q_v_max = q_avg * F_q
    
    lambda_tr = 0.0029
    D_c = lambda_tr / 3
    D_r = 0.16
    L_r = 2.85
    delta = D_c / D_r * L_r
    H_e = H_active + 1.42 * lambda_tr + 2 * delta
    
    h_profile = h_in + 1.0267 * (q_v_max * A_fuel * H_e) / (W_hc * np.pi) * \
                (np.sin(np.pi * z / H_e) + np.sin(np.pi * H_active / 2 / H_e))  # J/kg
    
    return h_profile, W_hc


# 4) TEMPERATURE PROFILE
def temperature_profile(h_profile, p_sys):
    """
    Calcola il profilo di temperatura del refrigerante lungo l'asse z.
    
    Parameters:
    - h_profile: profilo di entalpia specifica (J/kg)
    - p_sys: pressione del sistema (Pa)
    
    Returns:
    - T_profile: profilo di temperatura (°C)
    """
    T_profile = CP.PropsSI('T', 'H', h_profile, 'P', p_sys, 'Water') - 273.15  # °C
        # La temperatura non può superare quella di saturazione nel calcolo dell'entalpia
    T_sat = CP.PropsSI('T', 'P', p_sys, 'Q', 0, 'Water') - 273.15 # °C
    T_profile = np.minimum(T_profile, T_sat)
    return T_profile

# 5) EQUILIBRIUM QUALITY PROFILE
def equilibrium_quality_profile(h_profile, p_sys):
    """
    Calcola il profilo di titolo di equilibrio (equilibrium quality) lungo l'asse z.
    x_eq = (h - h_ls) / (h_vs - h_ls)
    """
    h_ls = CP.PropsSI('H', 'P', p_sys, 'Q', 0, 'Water')  # Entalpia liquido saturo
    h_vs = CP.PropsSI('H', 'P', p_sys, 'Q', 1, 'Water')  # Entalpia vapore saturo
    H_fg = h_vs - h_ls  # Calore latente di vaporizzazione
    x_eq_profile = (h_profile - h_ls) / (h_vs - h_ls)
    return x_eq_profile, H_fg

# 6) CALCULATION OF THE OUTER CLADDING TEMPERATURE

def safe_props(prop, T_K, P, T_sat_K):
        try:
            # Riduciamo leggermente T per evitare conflitti con la curva di saturazione
            if T_K >= T_sat_K - 0.01:
                return CP.PropsSI(prop, 'P', P, 'Q', 0, 'Water')
            return CP.PropsSI(prop, 'T', T_K, 'P', P, 'Water')
        except ValueError:
            return CP.PropsSI(prop, 'P', P, 'Q', 0, 'Water')
        
### CONTROLLARE PROFILO T, TROPPO VERSO IL BASSO
def T_outer_cladding_profile(T_sat_K, G_avg, D_eq, T_profile, p_sys, C, q_flux):
    """
    Calcola il coefficiente di scambio termico (h) usando la legge di Dittus-Boelter.
    
    """

    # Vettorizzazione iterando sui valori di T_profile (con conversione in Kelvin)
    mu_profile = np.array([safe_props('V', T + 273.15, p_sys, T_sat_K) for T in T_profile])
    k_profile = np.array([safe_props('L', T + 273.15, p_sys, T_sat_K) for T in T_profile])
    cp_profile = np.array([safe_props('C', T + 273.15, p_sys, T_sat_K) for T in T_profile])
    
    # Calcolo dei numeri adimensionali
    Re_profile = (G_avg * D_eq) / mu_profile
    Pr_profile = (cp_profile * mu_profile) / k_profile
    
    # Correlazione di Dittus-Boelter (riscaldamento del fluido, esponente Pr = 0.4)
    Nu_profile = C * (Re_profile**0.8) * (Pr_profile**0.4)
    
    # Calcolo di h
    h_single_phase = (Nu_profile * k_profile) / D_eq
    
    # Consideriamo boiling con Jens-Lottes
    q_flux_MW = q_flux*1e-6  # MW/m^2
    p_sys_bar = p_sys / 1e5  # bar
    T_co_JL = (T_sat_K - 273.15) + 25*(q_flux_MW)**0.25 * np.exp(-p_sys_bar/62) # °C (temperatura di ebollizione con superheat)

    # confronto con la temperatura della single phase
    T_co_SP = T_profile + q_flux / h_single_phase  # °C (temperatura di ebollizione con superheat)

    # combino i due profili prendendo il minimo tra i due (ovvero ho cambio di fase)
    T_co = np.minimum(T_co_JL, T_co_SP)

    # trovo l'indice in cui la curva Jens-Lottes diventa minore della Single Phase (punto di ONB - Onset of Nucleate Boiling)
    onb_indices = np.where(T_co_JL < T_co_SP)[0]
    first_onb_idx = onb_indices[0] if len(onb_indices) > 0 else None
    
    return h_single_phase, T_co, T_co_JL, T_co_SP, first_onb_idx

# bubble detachment
def detachment(h_l, q_D, T_sat_K, T_cool, zz):
    # converto T_sat in Celsius se T_cool e' in Celsius
    T_sat_C = T_sat_K - 273.15
    
    # T detachment: T_l è un array se q_D e h_l lo sono
    T_l = T_sat_C - q_D / (5 * h_l)

    # cerco gli indici in cui la condizione è soddisfatta
    # T_cool > T_l
    detachment_indices = np.where(T_cool > T_l)[0]
    
    if len(detachment_indices) > 0:
        # Prendo il primo punto in cui si verifica il distacco
        first_detachment_idx = detachment_indices[0]
        T_det = T_cool[first_detachment_idx]
        z_det = zz[first_detachment_idx]
        return T_det, z_det, first_detachment_idx
    else:
        # Nessun distacco lungo il canale
        return None, None, None

def flow_quality_two_phase(h_l, T_sat_K, T_cool, index_det, H_fg, p_sys, q_flux, zz, z_det, p_H, W):
    # calcoliamo ogni contributo (evaporation, convection e heat transfer)
    
    # taglio i vettori dal punto di detachment alla fine
    h_l = h_l[index_det:]
    T_cool = T_cool[index_det:]

    # bubbles heat transfer
    T_sat_C = T_sat_K - 273.15  # °C
    q_SP = h_l * (T_sat_C - T_cool)  # W/m^2 (heat flux per single phase)

    # Densità vapore saturo (Q = 1) - Costante a pressione di saturazione
    rho_g_sat = CP.PropsSI('D', 'T', T_sat_K, 'Q', 1, 'Water')
    i_l_sat = CP.PropsSI('H', 'T', T_sat_K, 'Q', 0, 'Water')  # J/kg (entalpia liquido saturo)
    # Densità liquido (in funzione della temperatura locale T_cool)
    rho_l = np.array([safe_props('D', T + 273.15, p_sys, T_sat_K) for T in T_cool])
    i_l = np.array([safe_props('H', T + 273.15, p_sys, T_sat_K) for T in T_cool])  # J/kg (entalpia liquido a T_cool)

    eps = rho_l /(rho_g_sat*H_fg)*(i_l_sat - i_l)  # quality evaporativa
    
    # calcolo la quality                                                                                                                 
    integr = p_H * (q_flux - q_SP) / (H_fg * W * (1 + eps))
    z_int = z[index_det:]
    integranda_int = integr[index_det:]
    x_z_parziale = cumulative_trapezoid(integranda_int, z_int, initial=0)
    x_z_totale = np.zeros_like(z)
    x_z_totale[index_det:] = x_z_parziale
    return x_z_totale

# ============================================================================
# MAIN
# ============================================================================

if __name__ == "__main__":
    
    # DATI DA TABELLA
    P_nom = 3400e6  # W
    p_sys = 15.51e6  # Pa
    n_rods = 157 * 264
    m_flow_tot = 1.4559e4  # kg/s
    F_q = 2.6  # heat flux hot channel factor
    m_flow_eff = 13456  # kg/s (portata effettiva considerando 5,9% di bypass flow)
    A_flow_eff = 3.883  # m^2 (area di flusso effettiva)
    T_sat = CP.PropsSI('T', 'P', p_sys, 'Q', 0, 'Water') #[K]

    # DATI GEOMETRICI (square array)
    Dout_clad = 9.5e-3  # m
    H_active = 4.2672  # m
    w = 12.6e-3  # m (pitch, passo tra le barre)
    s_guaina = 0.57e-3  # m (spessore guaina)
    D_in_clad = Dout_clad - 2 * s_guaina  # m (diametro interno)
    D_pellet = 8.22e-3  # m (diametro pellet)
    P_wet = np.pi * Dout_clad # m (perimetro bagnato)
    C = 0.042*w/Dout_clad - 0.024  # coefficiente per la correlazione di Dittus-Boelter

    # DISCRETIZZAZIONE DI Z (ASSIALE)
    z = np.linspace(-H_active/2, H_active/2, 100)  # m (asse z, da -H/2 a H/2)
    
    # CALCOLO AREA DI FLUSSO PER CANALE
    A_c = w**2 - np.pi/4 * Dout_clad**2  # m^2
    D_eq = 4 * A_c / P_wet  # m (diametro equivalente)

    # CHIAMATA ALLE FUNZIONI
    qv_profile = volumetric_heat_generation(z, P_nom, n_rods, D_in_clad, H_active, F_q)
    
    # calcolo il q_flux
    q_flux = qv_profile * A_c/P_wet  # W/m^2 (heat flux sulla guaina)

    G_avg = average_mass_velocity(m_flow_eff, A_flow_eff)
    
    h_profile, W_hc = coolant_specific_enthalpy(z, G_avg, A_c, w, Dout_clad, D_pellet, 
                                          p_sys, P_nom, n_rods, D_in_clad, H_active, F_q)
    
    T_profile = temperature_profile(h_profile, p_sys)

    x_eq_profile, H_fg = equilibrium_quality_profile(h_profile, p_sys)
    
    h_single_phase, T_co, T_co_JL, T_co_SP, first_onb_idx = T_outer_cladding_profile(T_sat, G_avg, D_eq, T_profile, p_sys, C, q_flux)
    
    T_det, z_det, first_detachment_idx= detachment(h_single_phase, q_flux, T_sat, T_profile, z)

    x_flow = flow_quality_two_phase(h_single_phase, T_sat, T_profile, first_detachment_idx, H_fg, p_sys, q_flux, z, z_det, p_sys, W_hc)

    # PLOT
    plt.figure(figsize=(10, 6))
    plt.plot(qv_profile, z)
    plt.title('Volumetric Heat Generation Rate along the z-axis')
    plt.ylabel('z (m)')
    plt.xlabel('q_v (W/m^3)')
    plt.grid()
    
    plt.figure(figsize=(10, 6))
    plt.plot(h_profile, z)
    plt.title('Coolant Specific Enthalpy Profile along the z-axis')
    plt.ylabel('z (m)')
    plt.xlabel('h (J/kg)')
    plt.grid()
    
    plt.figure(figsize=(10, 6))
    plt.plot(T_profile, z, label='T_coolant (°C)')
    if T_det is not None:
        plt.plot(T_det, z_det, 'ro', label=f'Detachment Point (T={T_det:.1f} °C, z={z_det:.2f} m)')
    plt.title('Coolant Temperature Profile along the z-axis')
    plt.ylabel('z (m)')
    plt.xlabel('T (°C)')
    plt.legend()
    plt.grid()
    
    plt.figure(figsize=(10, 6))
    plt.plot(T_co, z, label='T_co (Actual)', color='black', linewidth=2)
    plt.plot(T_co_JL, z, label='T_co_JL (Jens-Lottes)', linestyle='--')
    plt.plot(T_co_SP, z, label='T_co_SP (Single Phase)', linestyle='-.')
    plt.title('Outer Cladding Temperature Profile along the z-axis')
    plt.ylabel('z (m)')
    plt.xlabel('Temperature (°C)')
    plt.legend()
    plt.grid()
    
    plt.figure(figsize=(10, 6))
    plt.plot(x_flow, z)
    plt.title('Flow Quality Profile along the z-axis')
    plt.ylabel('z (m)')
    plt.xlabel('x (kg/kg)')
    plt.grid()

    plt.show()   