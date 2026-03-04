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
    """Calcola Reynolds e fattori d'attrito: smooth e rough."""
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

def calculate_friction_factor_colebrook(Re, eps, D, tolerance=1e-6, max_iterations=100):
    """Risolve l'equazione di Colebrook iterativamente partendo dal tubo liscio."""
    if Re <= 0:
        return 0
    
    # Stima iniziale 
    f = 64 / Re if Re < 3000 else 0.316 / (Re**0.25)
    
    # Iterazione 
    for _ in range(max_iterations):
        f_old = f
        sqrt_f = np.sqrt(f)
        term = (eps / D) / 3.7 + 2.51 / (Re * sqrt_f)
        f = 1 / (-2 * np.log10(term)) ** 2
        if abs(f - f_old) < tolerance:
            return f
    
    return f

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
    valid_od = None

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
                # Salva il primo diametro valido incontrato
                if valid_od is None:
                    valid_od = od_i
        else:
            L_req = calculate_L_for_h_equal_L(rl_i, rv_i, vl_i, vv_i, fl_i, fv_i, k_c, k_h, D_i, g)
            # Protezione contro i NoneType
            if L_req is not None and 0 <= L_req <= 100:
                plot_x.append(od_i)
                plot_y.append(L_req)
                # Salva il primo diametro valido incontrato
                if valid_od is None:
                    valid_od = od_i
                            
    return plot_x, plot_y, valid_od

def run_optimization_cycle_colebrook(data_table, Q_thermal, p, L_max, L_liq, L_vap, k_c, k_h, g, eps, find_L=False):
    """Stesso ciclo ma con Colebrook."""
    plot_x = []
    plot_y = []
    valid_od = None

    for row in data_table:
        od_i, th_i = row[0], row[1]
        D_i, A_i = get_pipe_geometry(od_i, th_i)
        m_i, vl_i, vv_i, rl_i, rv_i = calculate_massflow_and_velocity(Q_thermal, p, D_i, A_i)
        
        # Calcola Reynolds
        mu_liq = CP.PropsSI('V', 'P', p, 'Q', 0, 'Water')
        mu_vap = CP.PropsSI('V', 'P', p, 'Q', 1, 'Water')
        Re_liq = (rl_i * vl_i * D_i) / mu_liq
        Re_vap = (rv_i * vv_i * D_i) / mu_vap
        
        # Colebrook
        fl_i = calculate_friction_factor_colebrook(Re_liq, eps, D_i)
        fv_i = calculate_friction_factor_colebrook(Re_vap, eps, D_i)
        
        L_req = calculate_L_for_h_equal_L(rl_i, rv_i, vl_i, vv_i, fl_i, fv_i, k_c, k_h, D_i, g)
        if L_req is not None and 0 <= L_req <= 100:
            plot_x.append(od_i)
            plot_y.append(L_req)
            # Salva il primo diametro valido incontrato
            if valid_od is None:
                valid_od = od_i
    return plot_x, plot_y, valid_od

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

def plot_results_comparison(od_smooth, h_smooth, od_list_haaland, h_list_haaland, od_list_colebrook, h_list_colebrook, eps):
    """Grafico di confronto tra Haaland e Colebrook."""
    if not od_list_haaland or not od_list_colebrook:
        print("Nessun dato valido da plottare.")
        return
    
    ylabel = 'Required Length L [m]'
    title = 'Design Optimization with roughness : Length = Height'
    plt.figure(figsize=(10, 6))
    plt.plot(od_smooth, h_smooth, color='gray', linestyle='-', marker='^', label='Smooth', linewidth=2)
    plt.plot(od_list_haaland, h_list_haaland, color='green', linestyle='--', marker='o', label='Haaland', linewidth=2)
    plt.plot(od_list_colebrook, h_list_colebrook, color='orange', linestyle=':', marker='s', label='Colebrook', linewidth=2)
    plt.yscale('log')
    plt.xlabel('Outside Diameter [inch]')
    plt.ylabel(ylabel)
    plt.title(f'{title} (ε={eps:.2e} m)')
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

v_od, v_h, min_od_1 = run_optimization_cycle(table, Q_THERM, P_SYS, L_LIM, L_LIM, L_LIM, K_COLD, K_HOT, G_CONST, EPSILON, find_L=False)
plot_results(v_od, v_h, EPSILON, find_L=False)

# --- PUNTO 2 ---
h2 = solve_for_specific_diameter(OD_TEST, TH_TEST, Q_THERM, P_SYS, L_LIM, K_COLD, K_HOT, G_CONST, EPSILON, find_L=True)

v_od2, v_h2 ,min_od_2= run_optimization_cycle(table, Q_THERM, P_SYS, L_LIM, L_LIM, L_LIM, K_COLD, K_HOT, G_CONST, EPSILON, find_L=True)
plot_results(v_od2, v_h2, EPSILON, find_L=True)

# --- PUNTO 3 ---
EPSILON_3 = 5.0e-5 # Rugosità fornita  
h_3 = solve_for_specific_diameter(OD_TEST, TH_TEST, Q_THERM, P_SYS, L_LIM, K_COLD, K_HOT, G_CONST, EPSILON_3, find_L=True)

# Haaland
v_od3_h, v_h3_h, min_od_3_h = run_optimization_cycle(table, Q_THERM, P_SYS, L_LIM, L_LIM, L_LIM, K_COLD, K_HOT, G_CONST, EPSILON_3, find_L=True)

# Colebrook
v_od3_c, v_h3_c, min_od_3_c = run_optimization_cycle_colebrook(table, Q_THERM, P_SYS, L_LIM, L_LIM, L_LIM, K_COLD, K_HOT, G_CONST, EPSILON_3, find_L=True)

# Confronto
plot_results_comparison(v_od2, v_h2, v_od3_h, v_h3_h, v_od3_c, v_h3_c, EPSILON_3)

# --- GENERAZIONE TABELLA RIASSUNTIVA ---

print("\n" + "="*80)
print(f"{'TABELLA RIASSUNTIVA DEL PROGETTO':^80}")
print("="*80)

# Intestazione Colonne
header = f"{'Caso / Punto':<35} | {'Min OD Valido [in]':<20} | {f'h o L per OD={OD_TEST}\" [m]':<20}"
print(header)
print("-" * 80)

# Righe della tabella (3 cifre decimali per i risultati sull'OD specifico)
# Punto 1: Smooth - h calculation
print(f"{'Punto 1: Smooth (L = 10)':<35} | {min_od_1 if min_od_1 else 'N/A':<20} | {h_1:<20.3f}")

# Punto 2: Smooth - L critical
print(f"{'Punto 2: Smooth (L = h)':<35} | {min_od_2 if min_od_2 else 'N/A':<20} | {h2:<20.3f}")

# Punto 3: Rough - Haaland
print(f"{'Punto 3: Rough (L = h)':<35} | {min_od_3_h if min_od_3_h else 'N/A':<20} | {h_3:<20.3f}")

# Punto 3: Rough - Colebrook
# Assumiamo che h_3_c sia il risultato della funzione solve_for_specific_diameter usando il metodo Colebrook
# Se non l'hai calcolato a parte, puoi chiamare la funzione qui:
#print(f"{'Punto 3: Rough (Colebrook)':<35} | {min_od_3_c if min_od_3_c else 'N/A':<20} | {h_3_c:<20.3f}")

print("="*80)
54
plt.show()