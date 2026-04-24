import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP


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
    
    return h_profile


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
    
    x_eq_profile = (h_profile - h_ls) / (h_vs - h_ls)
    return x_eq_profile


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
    
    # DATI GEOMETRICI (square array)
    Dout_clad = 9.5e-3  # m
    H_active = 4.2672  # m
    w = 12.6e-3  # m (pitch, passo tra le barre)
    s_guaina = 0.57e-3  # m (spessore guaina)
    D_in_clad = Dout_clad - 2 * s_guaina  # m (diametro interno)
    D_pellet = 8.22e-3  # m (diametro pellet)
    
    # DISCRETIZZAZIONE DI Z (ASSIALE)
    z = np.linspace(-H_active/2, H_active/2, 100)  # m (asse z, da -H/2 a H/2)
    
    # CALCOLO AREA DI FLUSSO PER CANALE
    A_c = w**2 - np.pi/4 * Dout_clad**2  # m^2
    
    # CHIAMATA ALLE FUNZIONI
    qv_profile = volumetric_heat_generation(z, P_nom, n_rods, D_in_clad, H_active, F_q)
    
    G_avg = average_mass_velocity(m_flow_eff, A_flow_eff)
    
    h_profile = coolant_specific_enthalpy(z, G_avg, A_c, w, Dout_clad, D_pellet, 
                                          p_sys, P_nom, n_rods, D_in_clad, H_active, F_q)
    
    T_profile = temperature_profile(h_profile, p_sys)

    x_eq_profile = equilibrium_quality_profile(h_profile, p_sys)
    
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
    plt.plot(T_profile, z)
    plt.title('Coolant Temperature Profile along the z-axis')
    plt.ylabel('z (m)')
    plt.xlabel('T (°C)')
    plt.grid()
    plt.show()