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

def load_diameter_data(filepath):
    """Carica i dati dal file TXT."""
    try:
        return np.loadtxt(filepath)
    except FileNotFoundError:
        print(f"Errore: Il file {filepath} non esiste.")
        return None


def run_optimization_cycle(data_table, Q_thermal, p, L_max, L_liq, L_vap, k_c, k_h, g):
    """Esegue il ciclo su tutti i diametri e filtra i risultati."""
    plot_od = []
    plot_h = []

    for row in data_table:
        od_i, th_i = row[0], row[1]
        
        D_i, A_i = get_pipe_geometry(od_i, th_i)
        m_i, vl_i, vv_i, rl_i, rv_i = calculate_massflow_and_velocity(Q_thermal, p, D_i, A_i)
        fl_i, fv_i = calculate_friction_factors(p, vl_i, vv_i, D_i, rl_i, rv_i)
        h_i = calculate_height(rl_i, rv_i, vl_i, vv_i, fl_i, fv_i, k_c, k_h, L_liq, L_vap, D_i, g)
        
        if 0 < h_i <= L_max:
            plot_od.append(od_i)
            plot_h.append(h_i)
                        
    return plot_od, plot_h

def plot_results(od_list, h_list, L_max):
    """Genera il grafico dei diametri validi."""
    if not od_list:
        print("Nessun dato valido da plottare.")
        return

    plt.figure(figsize=(10, 6))
    plt.plot(od_list, h_list, color='navy', marker='o', label='Soluzioni Valide')
    plt.axhline(y=L_max, color='red', linestyle='--', label=f'Limite L = {L_max}m')
    
    plt.yscale('log')
    plt.xlabel('Outside Diameter [inch]')
    plt.ylabel('Altezza h [m]')
    plt.title('Dimensionamento: h vs OD')
    plt.grid(True, which='both', alpha=0.3)
    plt.legend()
    plt.show()

# --- MAIN SCRIPT ---

# Parametri fissi
G_CONST = 9.81
P_SYS = 70e5
Q_THERM = 34.8e6
L_LIMIT = 10.0
K_HOT, K_COLD = 40, 20

# Risolvo per un diametro specifico 
od_test = 16 # inch
th_test = 0.375 # inch
D_test, A_test = get_pipe_geometry(od_test, th_test)
m_test, vl_test, vv_test, rl_test, rv_test = calculate_massflow_and_velocity(Q_THERM, P_SYS, D_test, A_test)
fl_test, fv_test = calculate_friction_factors(P_SYS, vl_test, vv_test, D_test, rl_test, rv_test)
h_test = calculate_height(rl_test, rv_test, vl_test, vv_test, fl_test, fv_test, K_COLD, K_HOT, L_LIMIT, L_LIMIT, D_test, G_CONST)
print(f"Per OD={od_test} inch, h = {h_test:.2f} m")

# 1. Caricamento
table = load_diameter_data('assignment1/diameter_table.txt')

if table is not None:
    # 2. Ottimizzazione
    valid_od, valid_h = run_optimization_cycle(table, Q_THERM, P_SYS, L_LIMIT, L_LIMIT, L_LIMIT, K_COLD, K_HOT, G_CONST)
    
    # 3. Plot
    plot_results(valid_od, valid_h, L_LIMIT)



