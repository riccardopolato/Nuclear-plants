import CoolProp.CoolProp as CP
import numpy as np
import matplotlib.pyplot as plt
import csv
import os

def get_pipe_geometry(outside_diameter_inch, thickness_inch):
    """Calcola diametro interno [m] e area [m2] da pollici."""
    d_inner_m = (outside_diameter_inch - 2 * thickness_inch) * 0.0254
    area_m2 = np.pi * (d_inner_m / 2)**2
    return d_inner_m, area_m2

def reynolds_number(m, D_int, p, T):
    ' Calcola il numero di Reynolds'
    mu = CP.PropsSI('V', 'T', T+273.15, 'P', p, 'Water') # viscosità dinamica
    Re = 4*m / (np.pi*D_int*mu)
    return Re

def reynolds_number_shell(m, p, T, d_p, od_HX, D_eq, A_shell):
    ' Calcola il numero di Reynolds per il lato shell (considero la total cross section dei tubi moltiplicata per il numero di tubi)'
    mu = CP.PropsSI('V', 'T', T+273.15, 'P', p, 'Water') # viscosità dinamica
    Re = m*D_eq / (A_shell*mu) # controllare
    return Re

def friction_factor(eps, D, Re):
    ' Calcola il fattore di attrito con Haaland '
    term = ((eps/D)/3.7)**1.11 + 6.9/Re
    return (-1.8 * np.log10(term))**(-2)

def global_heat_transfer_coefficient_HX(P, A_tot, D_ext, D_int, k, Re, p, T, f,  D_eq, Re_shell, Circuit, F, p_shell=None, T_shell=None):
    ' Calcola U e il delta T logaritmico per lo scambiatore'
    
    # Proprietà del fluido lato tubi (da CoolProp)
    Pr = CP.PropsSI('Prandtl', 'P', p, 'T', T + 273.15, 'Water')
    k_fluid = CP.PropsSI('L', 'P', p, 'T', T + 273.15, 'Water')

    # Convezione esterna
    if Circuit == 'PSC':
        # h_ext è lato shell (acqua ISC): usa le proprietà del fluido shell
        Pr_shell = CP.PropsSI('Prandtl', 'P', p_shell, 'T', T_shell + 273.15, 'Water')
        h_ext = 0.351*Re_shell**0.55*k/D_eq*Pr_shell**(1/3)
    elif Circuit == 'ISC':
        flux = P / A_tot
        h_ext = flux**(1-1/3.86) * 2.257**(1/3.86)

    # Conduzione
    cond = D_ext*np.log(D_ext/D_int) / (2*k)

    # Convezione interna: Correlazione di Gnielinski
    if Re < 2300:
        Nu = 3.66
    else:
        Nu = (f/8)*(Re-1000)*Pr / (1 + 12.7*(f/8)**0.5*(Pr**(2/3)-1))
    h_int = Nu * k_fluid / D_int
    
    # metodo delle resistenze termiche in serie
    U = 1 / (1/h_ext + cond + D_ext/(D_int*h_int))
    deltaT_lm = P / (A_tot * U * F)
    return U, deltaT_lm

def temperature_ISC(T, p, T_ext, deltaT_lm, m, P_term):
    ' Calcola la temperatura nel ramo caldo e nel ramo freddo '
    deltaT = P_term / (m*CP.PropsSI('C', 'T', T+273.15, 'P', p, 'Water'))

    T_h = (T_ext - (deltaT+T_ext)*np.exp(deltaT/deltaT_lm)) / (1 - np.exp(deltaT/deltaT_lm))
    T_c = T_h - deltaT
    return T_c, T_h

def pressure_drop_buoyancy(p, T_h, T_c, H, g=9.81):
    ' Calcola la perdita di carico per differenza di densità'
    rho_h = CP.PropsSI('D', 'T', T_h+273.15, 'P', p, 'Water')
    rho_c = CP.PropsSI('D', 'T', T_c+273.15, 'P', p, 'Water')
    deltaP_buoyancy = (rho_c - rho_h) * g * H
    return deltaP_buoyancy

def pressure_drop_buoyancy_core(p, T_h, T_c, H2, g=9.81):
    ' Calcola la perdita di carico nel interno del core'
    rho_h = CP.PropsSI('D', 'T', T_h+273.15, 'P', p, 'Water')
    rho_c = CP.PropsSI('D', 'T', T_c+273.15, 'P', p, 'Water')
    # la densità varia linearmente rho(z) = rho_c + (rho_h - rho_c)*z/H2
    # per il calcolo della caduta di pressione si integra: dp=g*int(rho_c-rho(z))dz da 0 a H2
    deltaP_buoyancy_core = g * (rho_c - rho_h) / 2 * H2
    return deltaP_buoyancy_core

def pressure_drop_friction(m, N_HX, p, D_pipe, D_HX, T_av, T_c, T_h, eps_pipe, eps_HX, A_pipe, A_HX, L_ISC, L_HX,
                            N_bends, k_bends, k_shell, A_shell, Circuit, A_headers = 0.883, k_valve = 0.12): 
    ' Calcola il fattore delle perdite di carico che moltiplica per la portata'
    # perdite distibuite lungo i tubi
    Re_h = reynolds_number(m, D_pipe, p, T_h)
    Re_c = reynolds_number(m, D_pipe, p, T_c)
    f_h = friction_factor(eps_pipe, D_pipe, Re_h)
    f_c = friction_factor(eps_pipe, D_pipe, Re_c)
    dp_dist_c = f_c * L_ISC / D_pipe * 1/(2*CP.PropsSI('D', 'T', T_c+273.15, 'P', p, 'Water') * A_pipe**2)
    dp_dist_h = f_h * L_ISC / D_pipe * 1/(2*CP.PropsSI('D', 'T', T_h+273.15, 'P', p, 'Water') * A_pipe**2)

    # perdite distribuite lungo lo scambiatore (considero la total cross section dei tubi moltiplicata per il numero di tubi)
    Re_HX = reynolds_number(m, D_HX, p, T_av)
    f_HX = friction_factor(eps_HX, D_HX, Re_HX)
    dp_dist_HX = f_HX * L_HX / D_HX * 1/(2*CP.PropsSI('D', 'T', T_av+273.15, 'P', p, 'Water') * (A_HX*N_HX)**2)

    # perdite localizzate
    # gomiti
    dp_loc_bends_c = N_bends * k_bends * 1/(2*CP.PropsSI('D', 'T', T_c+273.15, 'P', p, 'Water') * A_pipe**2)
    dp_loc_bends_h = N_bends * k_bends * 1/(2*CP.PropsSI('D', 'T', T_h+273.15, 'P', p, 'Water') * A_pipe**2)

    # entrata e uscita
    A_HX_tot = A_HX * N_HX
    rho_av = CP.PropsSI('D', 'T', T_av+273.15, 'P', p, 'Water')

    if Circuit == 'ISC':
        # Transizione unica: main pipe <-> HX tubes
        if A_HX_tot < A_pipe:
            sigma = A_HX_tot / A_pipe
            A_narrow = A_HX_tot
        else:
            sigma = A_pipe / A_HX_tot
            A_narrow = A_pipe
        k_in = 0.5*(1-sigma)
        k_out = (1-sigma)**2
        dp_loc_HX = (k_in + k_out) * 1/(2*rho_av * A_narrow**2)
        k_in_header, k_out_header, dp_loc_pipe_header = 0.0, 0.0, 0.0

    elif Circuit == 'PSC':
        # Transizione 1: main pipe <-> header
        if A_pipe < A_headers:
            sigma_ph = A_pipe / A_headers
            A_narrow_ph = A_pipe
        else:
            sigma_ph = A_headers / A_pipe
            A_narrow_ph = A_headers
        k_in_header = 0.5*(1-sigma_ph)
        k_out_header = (1-sigma_ph)**2
        dp_loc_pipe_header = (k_in_header + k_out_header) * 1/(2*rho_av * A_narrow_ph**2)

        # Transizione 2: header <-> HX tubes
        if A_HX_tot < A_headers:
            sigma = A_HX_tot / A_headers
            A_narrow = A_HX_tot
        else:
            sigma = A_headers / A_HX_tot
            A_narrow = A_headers
        k_in = 0.5*(1-sigma)
        k_out = (1-sigma)**2
        dp_loc_HX = (k_in + k_out) * 1/(2*rho_av * A_narrow**2)

    #perdite nel shell solo lato ISC
    if Circuit == 'ISC':
        dp_loc_shell = k_shell * 1/(2*CP.PropsSI('D', 'T', T_av+273.15, 'P', p, 'Water') * A_shell**2)
        dp_loc_valve = 0
        dp_core = 0
        dp_bends_HX = 0
    elif Circuit == 'PSC':
        dp_loc_shell = 0
        dp_loc_valve = k_valve * 1/(2*CP.PropsSI('D', 'T', T_av+273.15, 'P', p, 'Water') * (A_pipe)**2)
        dp_core = 1.2e5/3200**2  # perdita di carico nel core: deltaP_core [Pa] / m^2 [kg^2/s^2]
        dp_bends_HX = k_bends * 1/(2*CP.PropsSI('D', 'T', T_av+273.15, 'P', p, 'Water') * (A_HX*N_HX)**2)
        

    dp_dict = {
        'dist_cold': dp_dist_c,
        'dist_hot': dp_dist_h,
        'dist_HX': dp_dist_HX,
        'loc_bends_hot': dp_loc_bends_h,
        'loc_bends_cold': dp_loc_bends_c,
        'loc_HX': dp_loc_HX,
        'loc_pipe_header': dp_loc_pipe_header,
        'loc_shell': dp_loc_shell,
        'loc_valve': dp_loc_valve,
        'loc_core': dp_core,
        'loc_bends_HX': dp_bends_HX
    }
    k_dict = {
        'k_shell': k_shell,
        'k_in': k_in,
        'k_out': k_out,
        'k_in_header': k_in_header,
        'k_out_header': k_out_header
    }
    dist_total = dp_dist_c + dp_dist_h + dp_dist_HX
    loc_total = dp_loc_bends_h + dp_loc_bends_c + dp_loc_HX + dp_loc_pipe_header + dp_loc_shell + dp_loc_valve + dp_core + dp_bends_HX
    deltaP_friction = dist_total + loc_total
    return deltaP_friction, dp_dict, dist_total, loc_total, k_dict

def mass_flow_rate(deltaP_buoyancy, deltaP_friction):
    ' Calcola la portata massica per una potenza termica data'
    return np.sqrt(deltaP_buoyancy / deltaP_friction)

def iteration(config, m_init=100, T_av_init=120, tolerance=1e-6, max_iter=100, shell_params=None):
    """
    Itera il calcolo della portata aggiornando T_av come media tra T_h e T_c.
    Se shell_params è fornito (dict con chiavi 'A_shell', 'D_eq', 'Re_shell', 'k_shell'), li usa invece di ricalcolare.
    """
    m = m_init
    T_av = T_av_init
    # Inizializza variabili shell-side
    A_shell = D_eq = Re_shell = k_shell = None
    p_shell = T_shell = None
    if shell_params is not None:
        A_shell = shell_params.get('A_shell')
        D_eq = shell_params.get('D_eq')
        Re_shell = shell_params.get('Re_shell')
        k_shell = shell_params.get('k_shell')
        p_shell = shell_params.get('pressure')
        T_shell = shell_params.get('T_av')
    elif config['Circuit'] == 'ISC':
        # A_shell e D_eq dipendono solo dalla geometria, sono costanti: calcolo una volta sola
        A_shell = config['D_shell']/config['d_pitch']*(config['d_pitch']-config['od_hx1'])*config['L_baffles']
        D_eq = (2*np.sqrt(3)*config['d_pitch']**2)/(np.pi*config['od_hx1'])-config['od_hx1']
    
    for iter_num in range(max_iter):
        # Calcolo del numero di Reynolds e fattore di attrito per lo scambiatore
        Re_av = reynolds_number(m, config['id_hx'], config['pressure'], T_av)
        f_av = friction_factor(config['eps_hx'], config['id_hx'], Re_av)
        
        # Calcolo Re_shell e k_shell nel lato shell (dipendono da m e T_av, si aggiornano ad ogni iterazione)
        if config['Circuit'] == 'ISC':
            Re_shell = reynolds_number_shell(
                m, config['pressure'], T_av, config['d_pitch'], config['od_hx1'], D_eq, A_shell)
            k_shell = 8*0.227/Re_shell**0.193*config['D_shell']/D_eq*(config['N_baffles']+1)
        # Per PSC, se shell_params è fornito, i valori sono già stati assegnati all'inizio e non serve riassegnare
        
        # Calcolo U e deltaT_lm nello scambiatore
        U, deltaT_lm = global_heat_transfer_coefficient_HX(
            config['P_term'], config['A_hx_tot'], config['od_hx'], config['id_hx'], 
            config['k_hx'], Re_av, config['pressure'], T_av, f_av,  D_eq, Re_shell, 
            config['Circuit'], config['Correction_f'], p_shell=p_shell, T_shell=T_shell)
        
        # Calcolo temperature caldo e freddo
        T_c, T_h = temperature_ISC(T_av, config['pressure'], config['T_HX'], deltaT_lm, m, config['P_term'])
        
        # Calcolo perdite di carico per galleggiabilità
        if config['Circuit'] == 'ISC':
            deltaP_buoyancy = pressure_drop_buoyancy(config['pressure'], T_h, T_c, config['H'])
        elif config['Circuit'] == 'PSC':
            deltaP_buoyancy = pressure_drop_buoyancy_core(config['pressure'], T_h, T_c, config['H'][1]) + pressure_drop_buoyancy(config['pressure'], T_h, T_c, config['H'][0])
        
        # Calcolo perdite di carico per attrito nel circuito completo
        deltaP_friction, dp_dict, dist_total, loc_total, k_dict = pressure_drop_friction(
            m, config['N_tubes'], config['pressure'], config['D_internal'], config['id_hx'],
            T_av, T_c, T_h, config['eps_pipe'], config['eps_hx'], config['A_pipe'],
            config['A_hx'], config['L_circuit'], config['L_hx'], config['N_bends'],
            config['k_bends'], k_shell, A_shell, config['Circuit'],
            A_headers=config.get('A_headers', 0.883)
        )
        
        # Calcolo nuova portata
        m_new = mass_flow_rate(deltaP_buoyancy, deltaP_friction)
        T_av_new = (T_h + T_c) / 2

        # Verifica di convergenza
        error = abs(m_new - m) / m
        if error < tolerance:
            print(f"\n✓ Convergenza raggiunta in {iter_num+1} iterazioni ({config['Circuit']})\n")
            # Se ISC, restituisci anche i parametri shell-side (incluse p e T_av per calcolo Pr_shell in PSC)
            if config['Circuit'] == 'ISC':
                shell_results = {
                    'A_shell': A_shell,
                    'D_eq': D_eq,
                    'Re_shell': Re_shell,
                    'k_shell': k_shell,
                    'pressure': config['pressure'],
                    'T_av': T_av_new
                }
                return m_new, T_av_new, T_h, T_c, deltaP_buoyancy, dp_dict, dist_total * m_new**2, loc_total * m_new**2, k_dict, shell_results
            else:
                return m_new, T_av_new, T_h, T_c, deltaP_buoyancy, dp_dict, dist_total * m_new**2, loc_total * m_new**2, k_dict
        
        m = m_new
        T_av = T_av_new
        
        if (iter_num + 1) % 10 == 0 or iter_num == 0:
            print(f"Iterazione {iter_num+1} ({config['Circuit']}): m = {m:.2f} kg/s, T_av = {T_av:.2f} °C, errore = {error:.4f}")
    
    print(f"\n⚠ Avviso: max iterazioni ({max_iter}) raggiunto senza convergenza ({config['Circuit']})\n")
    if config['Circuit'] == 'ISC':
        shell_results = {
            'A_shell': A_shell,
            'D_eq': D_eq,
            'Re_shell': Re_shell,
            'k_shell': k_shell,
            'pressure': config['pressure'],
            'T_av': T_av
        }
        return m, T_av, T_h, T_c, deltaP_buoyancy, dp_dict, dist_total * m**2, loc_total * m**2, k_dict, shell_results
    else:
        return m, T_av, T_h, T_c, deltaP_buoyancy, dp_dict, dist_total * m**2, loc_total * m**2, k_dict

def plot_pressure_drops(dp_dict_ISC, m_ISC, dist_ISC, loc_ISC, dp_b_ISC,
                        dp_dict_PSC, m_PSC, dist_PSC, loc_PSC, dp_b_PSC, out_dir):
    """Genera due grafici affiancati: dettaglio per componente e riepilogo aggregato."""

    # Etichette leggibili per ogni componente
    component_labels = [
        ('dist_cold',       'Distributed - cold pipe'),
        ('dist_hot',        'Distributed - hot pipe'),
        ('dist_HX',         'Distributed - HX tubes'),
        ('loc_bends_hot',   'Localized - bends hot'),
        ('loc_bends_cold',  'Localized - bends cold'),
        ('loc_HX',          'Localized - HX inlet/outlet'),
        ('loc_pipe_header', 'Localized - pipe↔header'),
        ('loc_shell',       'Localized - shell'),
        ('loc_valve',       'Localized - valve'),
        ('loc_core',        'Localized - core'),
        ('loc_bends_HX',    'Localized - bends HX'),
    ]
    keys   = [k for k, _ in component_labels]
    labels = [l for _, l in component_labels]

    vals_ISC = [dp_dict_ISC[k] * m_ISC**2 for k in keys]
    vals_PSC = [dp_dict_PSC[k] * m_PSC**2 for k in keys]

    color_ISC, color_PSC = 'steelblue', 'tomato'

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))
    fig.suptitle('Pressure drop comparison: ISC vs PSC', fontsize=14, fontweight='bold')

    # --- Grafico 1: dettaglio per componente (barre orizzontali) ---
    y = np.arange(len(keys))
    h = 0.35
    ax1.barh(y + h/2, vals_ISC, h, label='ISC', color=color_ISC)
    ax1.barh(y - h/2, vals_PSC, h, label='PSC', color=color_PSC)
    ax1.set_yticks(y)
    ax1.set_yticklabels(labels)
    ax1.set_xlabel('Pressure drop [Pa]')
    ax1.set_title('Detail by component')
    ax1.axvline(0, color='black', linewidth=0.7)
    ax1.legend()

    # --- Grafico 2: riepilogo aggregato ---
    categories = ['Distributed', 'Localized', 'Buoyancy\n(driving force)']
    vals_sum_ISC = [dist_ISC, loc_ISC, dp_b_ISC]
    vals_sum_PSC = [dist_PSC, loc_PSC, dp_b_PSC]

    x = np.arange(len(categories))
    w = 0.35
    bars_ISC = ax2.bar(x - w/2, vals_sum_ISC, w, label='ISC', color=color_ISC)
    bars_PSC = ax2.bar(x + w/2, vals_sum_PSC, w, label='PSC', color=color_PSC)
    ax2.set_xticks(x)
    ax2.set_xticklabels(categories)
    ax2.set_ylabel('Pressure drop [Pa]')
    ax2.set_title('Summary: distributed vs localized vs driving force')
    ax2.legend()
    # Aggiunge i valori sopra le barre
    for bar in bars_ISC:
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() * 1.01,
                 f'{bar.get_height():.0f}', ha='center', va='bottom', fontsize=8)
    for bar in bars_PSC:
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() * 1.01,
                 f'{bar.get_height():.0f}', ha='center', va='bottom', fontsize=8)

    plt.tight_layout()
    out_path = os.path.join(out_dir, 'pressure_drops_comparison.png')
    plt.savefig(out_path, dpi=150)
    
def save_results_to_csv(dp_dict, m, buoyancy, dist_total, loc_total, filename, k_dict=None):

    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["VOCE", "VALORE [Pa]"])
        writer.writerow(["--- DISTRIBUTED PRESSURE DROPS ---", ""])
        for k in ['dist_cold', 'dist_hot', 'dist_HX']:
            writer.writerow([k.replace('_', ' ').title(), round(dp_dict[k] * m**2, 4)])
        writer.writerow(["Total Distributed", round(dist_total, 4)])
        
        writer.writerow(["", ""])
        writer.writerow(["--- LOCALIZED PRESSURE DROPS ---", ""])
        for k in ['loc_bends_hot', 'loc_bends_cold', 'loc_HX', 'loc_pipe_header', 'loc_shell', 'loc_valve', 'loc_core', 'loc_bends_HX']:
            writer.writerow([k.replace('_', ' ').title(), round(dp_dict[k] * m**2, 4)])
        writer.writerow(["Total Localized", round(loc_total, 4)])
        
        writer.writerow(["", ""])
        writer.writerow(["--- FINALE RESUME ---", ""])
        writer.writerow(["Driving Force (Buoyancy)", round(buoyancy, 4)])
        writer.writerow(["Total Friction Losses", round(dist_total + loc_total, 4)])
        writer.writerow(["Mass Flow Rate [kg/s]", round(m, 4)])
        if k_dict is not None:
            writer.writerow(["" , ""])
            writer.writerow(["--- K LOSS COEFFICIENTS ---", ""])
            for k, v in k_dict.items():
                writer.writerow([k, round(v, 6)])


# ---------- DATA ----------
P_nom = 600e6 
P_term = 0.01 * P_nom 

# --- ISC ---
# Pipe
L_ISC, H_ISC, p_ISC = 40.0, 10.0, 70e5 
T_sat_ISC = CP.PropsSI('T', 'P', p_ISC, 'Q', 0, 'Water') - 273.15
od_ISC, th_ISC = 16.0, 1.031 
id_ISC, A_ISC = get_pipe_geometry(od_ISC, th_ISC)
eps_ISC, k_bends_ISC, N_bends_ISC = 2e-4 * id_ISC, 0.45, 2

# HX2
L_HX2, od_HX2, th_HX2, k_HX2, N_tubes_HX2 = 7.0, 25.4e-3, 1.24e-3, 15.0, 770
id_HX2 = od_HX2 - 2 * th_HX2
eps_HX2 = 1e-4 * id_HX2
A_HX2 = np.pi * (id_HX2/2)**2 
A_HX2_tot = 430.0 

# Pool
T_pool = CP.PropsSI('T', 'P', 1e5, 'Q', 0, 'Water') - 273.15


# --- PSC ---
# Pipe
L_PSC, H1_PSC, H2_PSC, p_PSC = 16.0, 7.0, 3.0, 75e5
T_sat_PSC = CP.PropsSI('T', 'P', p_PSC, 'Q', 0, 'Water') - 273.15
od_PSC, th_PSC = 16.0, 1.031
id_PSC, A_PSC = get_pipe_geometry(od_PSC, th_PSC)
eps_PSC, k_bends_PSC, N_bends_PSC , k_valve = 2e-4 * id_PSC, 0.45, 4, 0.12

# HX1
# U_tubes
L_HX1, od_HX1, th_HX1, k_HX1, N_tubes_HX1 = 9.314, 19.05e-3, 1.24e-3, 15.0, 897
id_HX1 = od_HX1 - 2 * th_HX1
eps_HX1 = 1e-4 * id_HX1
A_HX1 = np.pi * (id_HX1/2)**2
A_HX1_tot = 500.0 
A_headers = 0.883
F_T = 0.7
# Shell side
d_p = 28.5e-3
D_shell = 1.5
N_baffles, l_baffles = 2, 1.6



# CONFIGURAZIONE PSC
config_PSC = {
    'Circuit': 'PSC',
    'D_internal': id_PSC,
    'pressure': p_PSC,
    'T_sat': T_sat_PSC,
    'T_HX': 0, # da aggiornare con la temperatura media ottenuta dall'iterazione ISC
    'eps_pipe': eps_PSC,
    'eps_hx': eps_HX1,
    'A_pipe': A_PSC,
    'A_hx': A_HX1,
    'A_headers': A_headers,
    'od_hx1': od_HX1,
    'od_hx': od_HX1,
    'id_hx': id_HX1,
    'k_hx': k_HX1,
    'N_tubes': N_tubes_HX1,
    'A_hx_tot': A_HX1_tot,
    'L_circuit': L_PSC,
    'L_hx': L_HX1,
    'H': [H1_PSC , H2_PSC],
    'N_bends': N_bends_PSC,
    'k_bends': k_bends_PSC,
    'P_term': P_term,
    'Correction_f': F_T,
    'd_pitch': d_p,
    'D_shell': D_shell,
    'N_baffles': N_baffles,
    'L_baffles': l_baffles
}

# CONFIGURAZIONE ISC
config_ISC = {
    'Circuit': 'ISC',
    'D_internal': id_ISC,
    'pressure': p_ISC,
    'T_sat': T_sat_ISC,
    'T_HX': T_pool,
    'eps_pipe': eps_ISC,
    'eps_hx': eps_HX2,
    'A_pipe': A_ISC,
    'A_hx': A_HX2,
    'od_hx1': od_HX1,
    'od_hx': od_HX2,
    'id_hx': id_HX2,
    'k_hx': k_HX2,
    'N_tubes': N_tubes_HX2,
    'A_hx_tot': A_HX2_tot,
    'L_circuit': L_ISC,
    'L_hx': L_HX2,
    'H': H_ISC,
    'N_bends': N_bends_ISC,
    'k_bends': k_bends_ISC,
    'P_term': P_term,
    'Correction_f': 1.0,
    'd_pitch': d_p,
    'D_shell': D_shell,
    'N_baffles': N_baffles,
    'L_baffles': l_baffles
}
# ---------- EXECUTION ----------
# ITERAZIONI per ISC
m_res_ISC, T_av_res_ISC, T_h_res_ISC, T_c_res_ISC, dp_b_res_ISC, dp_dict_res_ISC, dist_total_ISC, loc_total_ISC, k_dict_ISC, shell_results_ISC = iteration(config_ISC)

# Stampa risultati finali organizzati
print("\n" + "="*80)
print("RISULTATI FINALI ITERAZIONE".center(80))
print("="*80)
print("\n📊 PARAMETRI TERMICI:")
print(f"{'Portata massica':<35} {m_res_ISC:>40.2f} kg/s")
print(f"{'Temperatura media (T_av)':<35} {T_av_res_ISC:>40.2f} °C")
print(f"{'Temperatura ramo caldo (T_h)':<35} {T_h_res_ISC:>40.2f} °C")
print(f"{'Temperatura ramo freddo (T_c)':<35} {T_c_res_ISC:>40.2f} °C")
print(f"{'Temperatura saturazione (limite)':<35} {T_sat_ISC:>40.2f} °C")
print(f"{'Distributed pressure drop':<35} {dist_total_ISC:>40.2f} Pa")
print(f"{'Localized pressure drop':<35} {loc_total_ISC:>40.2f} Pa")
print(f"{'Driving force (Buoyancy)':<35} {dp_b_res_ISC:>40.2f} Pa")
# Stampa tabella perdite di carico
save_results_to_csv(dp_dict_res_ISC, m_res_ISC, dp_b_res_ISC, dist_total_ISC, loc_total_ISC, filename=os.path.join(os.path.dirname(__file__), "result_ISC.csv"), k_dict=k_dict_ISC)
print("\n" + "="*80)
# ITERAZIONI per PSC, riutilizzando i parametri shell-side calcolati in ISC
config_PSC['T_HX'] = T_av_res_ISC # aggiorno la temperatura media ottenuta dall'iterazione ISC
m_res_PSC, T_av_res_PSC, T_h_res_PSC, T_c_res_PSC, dp_b_res_PSC, dp_dict_res_PSC, dist_total_PSC, loc_total_PSC, k_dict_PSC = iteration(config_PSC, shell_params=shell_results_ISC)

# Stampa risultati finali organizzati
print("\n" + "="*80)
print("RISULTATI FINALI ITERAZIONE".center(80))
print("="*80)
print("\n📊 PARAMETRI TERMICI:")
print(f"{'Portata massica':<35} {m_res_PSC:>40.2f} kg/s")
print(f"{'Temperatura media (T_av)':<35} {T_av_res_PSC:>40.2f} °C")
print(f"{'Temperatura ramo caldo (T_h)':<35} {T_h_res_PSC:>40.2f} °C")
print(f"{'Temperatura ramo freddo (T_c)':<35} {T_c_res_PSC:>40.2f} °C")
print(f"{'Temperatura saturazione (limite)':<35} {T_sat_PSC:>40.2f} °C")
print(f"{'Distributed pressure drop':<35} {dist_total_PSC:>40.2f} Pa")
print(f"{'Localized pressure drop':<35} {loc_total_PSC:>40.2f} Pa")
print(f"{'Driving force (Buoyancy)':<35} {dp_b_res_PSC:>40.2f} Pa")
# Stampa tabella perdite di carico
save_results_to_csv(dp_dict_res_PSC, m_res_PSC, dp_b_res_PSC, dist_total_PSC, loc_total_PSC, filename=os.path.join(os.path.dirname(__file__), "result_PSC.csv"), k_dict=k_dict_PSC)

# ---------- PLOT ----------
plot_pressure_drops(
    dp_dict_res_ISC, m_res_ISC, dist_total_ISC, loc_total_ISC, dp_b_res_ISC,
    dp_dict_res_PSC, m_res_PSC, dist_total_PSC, loc_total_PSC, dp_b_res_PSC,
    out_dir=os.path.dirname(__file__)
)
