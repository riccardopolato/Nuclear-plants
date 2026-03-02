import CoolProp.CoolProp as CP
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

def get_pipe_geometry(outside_diameter_inch, thickness_inch):
    """Calcola diametro interno [m] e area [m2] da pollici."""
    # uso la tabella scegliendo sempre da STD
    d_inner_m = (outside_diameter_inch - 2 * thickness_inch) * 0.0254
    area_m2 = np.pi * (d_inner_m / 2)**2
    return d_inner_m, area_m2

def calculate_massflow_and_velocity(Q_total, p, D, A):
    """Calcola portata massica [kg/s] e velocità delle fasi [m/s]."""
    h_liq = CP.PropsSI('H', 'P', p, 'Q', 0, 'Water')
    h_vap = CP.PropsSI('H', 'P', p, 'Q', 1, 'Water')
    # dato che sono a saturazione a p e T costanti, Q_thermico = m_dot * (h_vap - h_liq) 
    m_dot = Q_total / (h_vap - h_liq)
    
    rho_l = CP.PropsSI('D', 'P', p, 'Q', 0, 'Water')
    rho_v = CP.PropsSI('D', 'P', p, 'Q', 1, 'Water')
    
    v_liq = m_dot / (rho_l * A)
    v_vap = m_dot / (rho_v * A)
    
    return m_dot, v_liq, v_vap

def calculate_fluid_properties_and_friction(p, v_liq, v_vap, D):
    """Calcola densità, Reynolds e fattori d'attrito (per ora solo caso smooth)."""
    rho_liq = CP.PropsSI('D', 'P', p, 'Q', 0, 'Water')
    rho_vap = CP.PropsSI('D', 'P', p, 'Q', 1, 'Water')
    mu_liq = CP.PropsSI('V', 'P', p, 'Q', 0, 'Water')
    mu_vap = CP.PropsSI('V', 'P', p, 'Q', 1, 'Water')
    
    Re_liq = (rho_liq * v_liq * D) / mu_liq
    Re_vap = (rho_vap * v_vap * D) / mu_vap
    
    # calcolo il reynolds distinguendo tr laminare e turbolento
    def get_f(Re):
        if Re <= 0: return 0
        return 64 / Re if Re < 3000 else 0.316 / (Re**0.25)
    
    return mu_liq, mu_vap, rho_liq, rho_vap, Re_liq, Re_vap, get_f(Re_liq), get_f(Re_vap)

def calculate_height(rho_l, rho_v, v_l, v_v, f_l, f_v, k_l, k_v, L_l, L_v, D, g=9.81):
    """Risolve il bilancio di pressione per trovare l'altezza h [m]."""
    term_liquid = (f_l * (L_l / D) + k_l) * (0.5 * rho_l * v_l**2)
    term_vapor = (f_v * (L_v / D) + k_v) * (0.5 * rho_v * v_v**2)
    
    h = (term_liquid + term_vapor) / ((rho_l - rho_v) * g)
    return h

# --- SEZIONE DATI  ---

# Costanti fisiche
g = 9.81              # m/s^2

# Parametri termoidraulici
p = 70 * 10**5        # Pa (70 bar)
Q_thermal = 34.8e6    # W (Potenza termica)

# Geometria del circuito
L_liq = 10.0          # m (Lunghezza ramo liquido)
L_vap = 10.0          # m (Lunghezza ramo vapore - ipotizzata uguale a L1)

# Geometria del tubo (da pollici a metri)
OD_inch = 16         # inches
th_inch = 0.375      # inches

# Coefficienti di perdita localizzata
k_h = 40              # heater
k_c = 20              # cooler

# --- ESECUZIONE CALCOLI ---

# 1. Geometria
D, A = get_pipe_geometry(OD_inch, th_inch)

# 2. Portata e Velocità
m_dot, v_l, v_v = calculate_massflow_and_velocity(Q_thermal, p, D, A)

# 3. Proprietà e Attrito
mu_l, mu_v, rho_l, rho_v, Re_l, Re_v, f_l, f_v = calculate_fluid_properties_and_friction(p, v_l, v_v, D)

# 4. Calcolo Altezza h
h_finale = calculate_height(rho_l, rho_v, v_l, v_v, f_l, f_v, k_c, k_h, L_liq, L_vap, D, g)

# --- OUTPUT RISULTATI ---

print(f"Altezza h necessaria: {h_finale:.4f} m")
