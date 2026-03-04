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

def calculate_friction_factors(p, v_liq, v_vap, D, rho_l, rho_v, eps):
    """Calcola Reynolds e fattori d'attrito."""
    mu_liq = CP.PropsSI('V', 'P', p, 'Q', 0, 'Water')
    mu_vap = CP.PropsSI('V', 'P', p, 'Q', 1, 'Water')
    
    Re_liq = (rho_l * v_liq * D) / mu_liq
    Re_vap = (rho_v * v_vap * D) / mu_vap
    
    def get_f(Re):
        if eps > 0:
             term = ((eps/D)/3.7)**1.11 + 6.9/Re
             return (-1.8 * np.log10(term))**-2
        if eps == 0:
            if Re <= 0: return 0
        return 64 / Re if Re < 3000 else 0.316 / (Re**0.25)
    
    return get_f(Re_liq), get_f(Re_vap)

def calculate_height(rho_l, rho_v, v_l, v_v, f_l, f_v, k_l, k_v, L_l, L_v, D, g=9.81):
    """Bilancio di pressione per altezza h."""
    term_liquid = (f_l * (L_l / D) + k_l) * (0.5 * rho_l * v_l**2)
    term_vapor = (f_v * (L_v / D) + k_v) * (0.5 * rho_v * v_v**2)
    h = (term_liquid + term_vapor) / ((rho_l - rho_v) * g)
    return h

def calculate_L_for_h_equal_L(rho_l, rho_v, v_l, v_v, f_l, f_v, k_l, k_v, D, g=9.81):
    """Risolve analiticamente per L imponendo h = L (Punto 2)."""
    dp_conc = (k_l * 0.5 * rho_l * v_l**2) + (k_v * 0.5 * rho_v * v_v**2)
    fric_unit = (f_l / D * 0.5 * rho_l * v_l**2) + (f_v / D * 0.5 * rho_v * v_v**2)
    buoyancy_unit = (rho_l - rho_v) * g
    
    denom = buoyancy_unit - fric_unit
    if denom <= 0:
        return None 
    return dp_conc / denom

def load_diameter_data(filepath):
    """Carica i dati dal file TXT."""
    try:
        return np.loadtxt(filepath)
    except FileNotFoundError:
        print(f"Errore: Il file {filepath} non esiste.")
        return None

def solve_for_specific_diameter(od_inch, th_inch, Q_thermal, p, L_max, k_c, k_h, g, eps, find_L=False):
    """Calcola h per un singolo diametro specifico inserito dall'utente."""
    # Sfrutta le funzioni esistenti
    D, A = get_pipe_geometry(od_inch, th_inch)
    m_dot, vl, vv, rl, rv = calculate_massflow_and_velocity(Q_thermal, p, D, A)
    fl, fv = calculate_friction_factors(p, vl, vv, D, rl, rv, eps)
    
    if not find_L:
        h = calculate_height(rl, rv, vl, vv, fl, fv, k_c, k_h, L_max, L_max, D, g)
        return h
    else:
        L = calculate_L_for_h_equal_L(rl, rv, vl, vv, fl, fv, k_c, k_h, D, g)
        return L

def run_optimization_cycle(data_table, Q_thermal, p, L_max, L_liq, L_vap, k_c, k_h, g, eps, find_L=False):
    plot_x = []
    plot_y = []

    for row in data_table:
        od_i, th_i = row[0], row[1]
        D_i, A_i = get_pipe_geometry(od_i, th_i)
        m_i, vl_i, vv_i, rl_i, rv_i = calculate_massflow_and_velocity(Q_thermal, p, D_i, A_i)
        fl_i, fv_i = calculate_friction_factors(p, vl_i, vv_i, D_i, rl_i, rv_i, eps)
        
        if not find_L:
            h_i = calculate_height(rl_i, rv_i, vl_i, vv_i, fl_i, fv_i, k_c, k_h, L_liq, L_vap, D_i, g)
            if 0 < h_i <= L_max:
                plot_x.append(od_i)
                plot_y.append(h_i)
        else:
            L_req = calculate_L_for_h_equal_L(rl_i, rv_i, vl_i, vv_i, fl_i, fv_i, k_c, k_h, D_i, g)
            # Protezione contro i NoneType
            if L_req is not None and 1 <= L_req <= 100:
                plot_x.append(od_i)
                plot_y.append(L_req)
                        
    return plot_x, plot_y

def plot_results(od_list, h_list, eps, find_L=False):
    """Genera il grafico dei diametri validi."""
    if not od_list:
        print("Nessun dato valido da plottare.")
        return
    
    # Gestione etichetta rugosità
    condition_label = "Smooth" if eps == 0 else f"Roughness = {eps:.2e} m)"

    plt.figure(figsize=(10, 6))
    plt.plot(od_list, h_list, color='navy', marker='o', label=condition_label)
    plt.yscale('log')
    plt.xlabel('Outside Diameter [inch]')
    plt.ylabel('Height h [m]' if not find_L else 'Required Length L [m]')
    plt.title('Design Optimization: Height vs Diameter' if not find_L else 'Design Optimization: Length = Height')
    plt.grid(True, which='both', alpha=0.3)
    plt.legend()
    
# --- MAIN SCRIPT ---

G_CONST = 9.81 
P_SYS = 70e5 
Q_THERM = 34.8e6 
K_HOT, K_COLD = 40, 20

# Dati di input da tabella
OD_TEST = 14 
TH_TEST = 0.375 
table = load_diameter_data('assignment1/diameter_table.txt')

# --- PUNTO 1 ---
L_LIM = 10.0 
EPSILON = 0  # Tubi lisci

h_1 = solve_for_specific_diameter(OD_TEST, TH_TEST, Q_THERM, P_SYS, L_LIM, K_COLD, K_HOT, G_CONST, EPSILON, find_L=False)
print(f"For smooth pipe OD={OD_TEST} inch, h = {h_1:.2f} m")

v_od, v_h = run_optimization_cycle(table, Q_THERM, P_SYS, L_LIM, L_LIM, L_LIM, K_COLD, K_HOT, G_CONST, EPSILON, find_L=False)
plot_results(v_od, v_h, EPSILON, find_L=False)

# --- PUNTO 2 ---
h2 = solve_for_specific_diameter(OD_TEST, TH_TEST, Q_THERM, P_SYS, L_LIM, K_COLD, K_HOT, G_CONST, EPSILON, find_L=True)
print(f"For smooth pipe OD={OD_TEST} inch, L richiesto per h=L è: {h2:.2f} m")

v_od2, v_L2 = run_optimization_cycle(table, Q_THERM, P_SYS, L_LIM, L_LIM, L_LIM, K_COLD, K_HOT, G_CONST, EPSILON, find_L=True)
plot_results(v_od2, v_L2, EPSILON, find_L=True)

# --- PUNTO 3 ---
EPSILON_3 = 5.0e-5 # Rugosità fornita  
h_1 = solve_for_specific_diameter(OD_TEST, TH_TEST, Q_THERM, P_SYS, L_LIM, K_COLD, K_HOT, G_CONST, EPSILON_3, find_L=True)
print(f"For rough pipe OD={OD_TEST} inch, h = {h_1:.2f} m")

v_od, v_h = run_optimization_cycle(table, Q_THERM, P_SYS, L_LIM, L_LIM, L_LIM, K_COLD, K_HOT, G_CONST, EPSILON_3, find_L=True)
plot_results(v_od, v_h, EPSILON_3, find_L=True)
plt.show()