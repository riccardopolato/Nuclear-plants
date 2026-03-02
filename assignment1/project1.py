import CoolProp.CoolProp as CP
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

# --- FUNZIONI DI CALCOLO ---

def get_pipe_geometry(outside_diameter_inch, thickness_inch):
    """Calcola diametro interno [m] e area [m2] da pollici."""
    d_inner_m = (outside_diameter_inch - 2 * thickness_inch) * 0.0254
    area_m2 = np.pi * (d_inner_m / 2)**2
    return d_inner_m, area_m2

def calculate_massflow_and_velocity(Q_total, p, D, A):
    """Calcolo termoidraulico delle fasi a saturazione."""
    h_liq = CP.PropsSI('H', 'P', p, 'Q', 0, 'Water')
    h_vap = CP.PropsSI('H', 'P', p, 'Q', 1, 'Water')
    m_dot = Q_total / (h_vap - h_liq)
    
    rho_l = CP.PropsSI('D', 'P', p, 'Q', 0, 'Water')
    rho_v = CP.PropsSI('D', 'P', p, 'Q', 1, 'Water')
    
    v_liq = m_dot / (rho_l * A)
    v_vap = m_dot / (rho_v * A)
    return m_dot, v_liq, v_vap, rho_l, rho_v

def calculate_friction_factors(p, v_liq, v_vap, D, rho_l, rho_v):
    """Calcola Reynolds e fattori d'attrito."""
    mu_liq = CP.PropsSI('V', 'P', p, 'Q', 0, 'Water')
    mu_vap = CP.PropsSI('V', 'P', p, 'Q', 1, 'Water')
    
    Re_liq = (rho_l * v_liq * D) / mu_liq
    Re_vap = (rho_v * v_vap * D) / mu_vap
    
    def get_f(Re):
        if Re <= 0: return 0
        return 64 / Re if Re < 3000 else 0.316 / (Re**0.25)
    
    return get_f(Re_liq), get_f(Re_vap)

def calculate_height(rho_l, rho_v, v_l, v_v, f_l, f_v, k_l, k_v, L_l, L_v, D, g=9.81):
    """Bilancio di pressione per altezza h."""
    term_liquid = (f_l * (L_l / D) + k_l) * (0.5 * rho_l * v_l**2)
    term_vapor = (f_v * (L_v / D) + k_v) * (0.5 * rho_v * v_v**2)
    h = (term_liquid + term_vapor) / ((rho_l - rho_v) * g)
    return h

# --- SEZIONE DATI ---

g = 9.81
p = 70 * 10**5
Q_thermal = 34.8e6
L_liq = 10.0
L_vap = 10.0
L_max = 10.0  # Limite fisico h <= L
k_h = 40      # Heater (hot)
k_c = 20      # Cooler (cold)

# --- 1. CARICAMENTO DATI DA FILE ESTERNO ---
# Carichiamo le due colonne: Outside Diameter e Thickness
data_table = np.loadtxt('assignment1/diameter_table.txt')

# --- 2. CALCOLO SINGOLO (Esempio per OD = 16) ---
OD_test, th_test = 16.0, 0.375
D_t, A_t = get_pipe_geometry(OD_test, th_test)
m_t, vl_t, vv_t, rl_t, rv_t = calculate_massflow_and_velocity(Q_thermal, p, D_t, A_t)
fl_t, fv_t = calculate_friction_factors(p, vl_t, vv_t, D_t, rl_t, rv_t)
h_test = calculate_height(rl_t, rv_t, vl_t, vv_t, fl_t, fv_t, k_c, k_h, L_liq, L_vap, D_t, g)

print(f"--- TEST SINGOLO (OD {OD_test}\") ---")
print(f"h calcolato: {h_test:.4f} m\n")

# --- 3. CICLO DI OTTIMIZZAZIONE ---
plot_od = []
plot_h = []

print(f"{'OD [in]':<10} | {'h [m]':<10} | {'Stato'}")
print("-" * 35)

for row in data_table:
    od_i = row[0]
    th_i = row[1]
    
    D_i, A_i = get_pipe_geometry(od_i, th_i)
    m_i, vl_i, vv_i, rl_i, rv_i = calculate_massflow_and_velocity(Q_thermal, p, D_i, A_i)
    fl_i, fv_i = calculate_friction_factors(p, vl_i, vv_i, D_i, rl_i, rv_i)
    h_i = calculate_height(rl_i, rv_i, vl_i, vv_i, fl_i, fv_i, k_c, k_h, L_liq, L_vap, D_i, g)
    
    if 0 < h_i <= L_max:
        plot_od.append(od_i)
        plot_h.append(h_i)
        print(f"{od_i:<10.3f} | {h_i:<10.4f} | ✅ OK")
    else:
        print(f"{od_i:<10.3f} | {h_i:<10.4f} | ❌ SCARTATO")

# --- 4. GRAFICO FINALE ---
if plot_od:
    plt.figure(figsize=(10, 6))
    plt.plot(plot_od, plot_h, color='navy', marker='o', linestyle='-', linewidth=2, markersize=6)
    plt.axhline(y=L_max, color='red', linestyle='--', label='Limite h = L')
    
    plt.yscale('log') # Scala logaritmica utile per range ampi di h
    plt.xlabel('Outside Diameter [inch]')
    plt.ylabel('Altezza h [m] (Scala Log)')
    plt.title('Dimensionamento Ottimale: h vs Diametro Esterno')
    plt.grid(True, which='both', alpha=0.3)
    plt.legend()
    plt.show()