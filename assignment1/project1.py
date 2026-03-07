import CoolProp.CoolProp as CP
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

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

    
    plt.subplot(2, 1, 1 if not find_L else 2)
    plt.plot(od_list, h_list, color='black' if not find_L else 'navy', marker='o', label=condition_label)
    if not find_L:
        plt.axhline(y=10, color='red', linestyle='--', linewidth=1.5, label='Limit L = 10m', zorder=2)
    plt.yscale('log')
    plt.xlabel('Outside Diameter [inch]')
    plt.ylabel('Height h [m]' if not find_L else 'Required Length L [m]')
    plt.title('Design Optimization: Height vs Diameter' if not find_L else 'Design Optimization: Length = Height')
    plt.grid(True, which='both', alpha=0.3)
    plt.legend(loc='upper right')

def plot_results_comparison(od_smooth, h_smooth, od_list_haaland, h_list_haaland, od_list_colebrook, h_list_colebrook, eps):
    """Grafico di confronto tra Haaland e Colebrook con Inset Plot (Zoom)."""
    if not od_list_haaland or not od_list_colebrook:
        print("Nessun dato valido da plottare.")
        return
    
    fig, ax = plt.subplots(figsize=(12, 7))
    
    # --- PLOT PRINCIPALE ---
    ax.plot(od_smooth, h_smooth, color='gray', linestyle='-', marker='^', label='Smooth', linewidth=1.5, alpha=0.6)
    ax.plot(od_list_haaland, h_list_haaland, color='navy', linestyle='--', marker='o', label='Haaland', linewidth=2)
    ax.plot(od_list_colebrook, h_list_colebrook, color='darkorange', linestyle=':', marker='s', label='Colebrook', linewidth=2)
    
    ax.set_yscale('log')
    ax.set_xlabel('Outside Diameter [inch]', fontweight='bold')
    ax.set_ylabel('Required Length L [m]', fontweight='bold')
    ax.set_title(f'Design Optimization: Length = Height (ε={eps:.2e} m)', fontsize=14)
    ax.grid(True, which='both', alpha=0.3)
    ax.legend(loc='lower left', fontsize=10)

    # --- CREAZIONE INSET (ZOOM) ---
    # Posizioniamo il riquadro in una zona vuota (es. in basso a sinistra 'lower left')
    # width e height sono percentuali del grafico principale
    ax_ins = inset_axes(ax, width="45%", height="45%", loc='upper right', borderpad=3)
    
    # Plottiamo gli stessi dati nell'inset
    ax_ins.plot(od_smooth, h_smooth, color='gray', linestyle='-', marker='^', alpha=0.6)
    ax_ins.plot(od_list_haaland, h_list_haaland, color='navy', linestyle='--', marker='o')
    ax_ins.plot(od_list_colebrook, h_list_colebrook, color='darkorange', linestyle=':', marker='s')
    
    # Definiamo i limiti dello zoom (8" - 14")
    x1, x2 = 8, 14
    # Troviamo i valori di y corrispondenti per centrare lo zoom (usiamo i valori di Colebrook come riferimento)
    # Filtriamo i valori di h nell'intervallo x1-x2
    y_vals = [h for x, h in zip(od_list_colebrook, h_list_colebrook) if x1 <= x <= x2]
    y1, y2 = min(y_vals)*0.9, max(y_vals)*1.1
    
    ax_ins.set_xlim(x1, x2)
    ax_ins.set_ylim(y1, y2)
    ax_ins.set_yscale('log') # Fondamentale mantenere la scala logaritmica anche qui
    ax_ins.grid(True, alpha=0.2)
    
    # Riduciamo la dimensione dei font per l'inset
    ax_ins.tick_params(axis='both', which='major', labelsize=8)
    
    # Aggiungiamo i connettori (le linee che uniscono il grafico principale allo zoom)
    mark_inset(ax, ax_ins, loc1=1, loc2=3, fc="none", ec="0.5", linestyle='--')

    plt.draw()

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
v_od, v_h, min_od_1 = run_optimization_cycle(table, Q_THERM, P_SYS, L_LIM, L_LIM, L_LIM, K_COLD, K_HOT, G_CONST, EPSILON, find_L=False)

# --- PUNTO 2 ---
v_od2, v_h2 ,min_od_2= run_optimization_cycle(table, Q_THERM, P_SYS, L_LIM, L_LIM, L_LIM, K_COLD, K_HOT, G_CONST, EPSILON, find_L=True)

# --- PUNTO 3 ---
EPSILON_3 = 5.0e-5 # Rugosità fornita  
# Haaland
v_od3_h, v_h3_h, min_od_3_h = run_optimization_cycle(table, Q_THERM, P_SYS, L_LIM, L_LIM, L_LIM, K_COLD, K_HOT, G_CONST, EPSILON_3, find_L=True)
# Colebrook
v_od3_c, v_h3_c, min_od_3_c = run_optimization_cycle_colebrook(table, Q_THERM, P_SYS, L_LIM, L_LIM, L_LIM, K_COLD, K_HOT, G_CONST, EPSILON_3, find_L=True)

# --- VISUALIZZAZIONE GRAFICA ---
plt.figure(figsize=(10, 10))
plot_results(v_od, v_h, EPSILON, find_L=False)
plot_results(v_od2, v_h2, EPSILON, find_L=True)
plt.tight_layout()

# Grafico Punto 3
plot_results_comparison(v_od2, v_h2, v_od3_h, v_h3_h, v_od3_c, v_h3_c, EPSILON_3)


# --- CONFIGURAZIONE DIAMETRI SPECIFICI ---
# Lista dei diametri richiesti: (OD, Thickness)
target_diameters = [
    (8.625, 0.322),
    (10.750, 0.365),
    (12.750, 0.375),
    (14.000, 0.375)
]

# --- ESECUZIONE CALCOLI E TABELLA ---

print("\n" + "="*110)
print(f"{'TABELLA COMPARATIVA DIAMETRI SELEZIONATI':^110}")
print("="*110)

# Intestazione Colonne
header = f"{'OD [in]':<10} | {'Th [in]':<10} | {'h (Punto 1) [m]':<18} | {'L=h (Punto 2) [m]':<18} | {'L=h Rough (H) [m]':<18} | {'L=h Rough (C) [m]':<18}"
print(header)
print("-" * 110)

for od, th in target_diameters:
    # Punto 1: Smooth, L=10, trova h
    h1_val = solve_for_specific_diameter(od, th, Q_THERM, P_SYS, L_LIM, K_COLD, K_HOT, G_CONST, EPSILON, find_L=False)
    
    # Punto 2: Smooth, trova L=h
    h2_val = solve_for_specific_diameter(od, th, Q_THERM, P_SYS, L_LIM, K_COLD, K_HOT, G_CONST, EPSILON, find_L=True)
    
    # Punto 3: Rough (Haaland), trova L=h
    h3_h_val = solve_for_specific_diameter(od, th, Q_THERM, P_SYS, L_LIM, K_COLD, K_HOT, G_CONST, EPSILON_3, find_L=True)
    
    # Punto 3: Rough (Colebrook), trova L=h
    # Calcoliamo Colebrook direttamente per il diametro specifico
    D_spec, A_spec = get_pipe_geometry(od, th)
    m_spec, vl_spec, vv_spec, rl_spec, rv_spec = calculate_massflow_and_velocity(Q_THERM, P_SYS, D_spec, A_spec)
    
    mu_l = CP.PropsSI('V', 'P', P_SYS, 'Q', 0, 'Water')
    mu_v = CP.PropsSI('V', 'P', P_SYS, 'Q', 1, 'Water')
    re_l = (rl_spec * vl_spec * D_spec) / mu_l
    re_v = (rv_spec * vv_spec * D_spec) / mu_v
    
    fl_c = calculate_friction_factor_colebrook(re_l, EPSILON_3, D_spec)
    fv_c = calculate_friction_factor_colebrook(re_v, EPSILON_3, D_spec)
    
    h3_c_val = calculate_L_for_h_equal_L(rl_spec, rv_spec, vl_spec, vv_spec, fl_c, fv_c, K_COLD, K_HOT, D_spec, G_CONST)

    # Formattazione riga (N/A se il diametro non è fisicamente possibile)
    def fmt(val): return f"{val:18.3f}" if val and val > 0 else f"{'N/A':>18}"

    print(f"{od:<10.3f} | {th:<10.3f} | {fmt(h1_val)} | {fmt(h2_val)} | {fmt(h3_h_val)} | {fmt(h3_c_val)}")

print("="*110)
plt.show()