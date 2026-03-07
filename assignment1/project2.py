import CoolProp.CoolProp as CP
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
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

def friction_factor(eps, D, Re):
    ' Calcola il fattore di attrito con Haaland '
    term = ((eps/D)/3.7)**1.11 + 6.9/Re
    return (-1.8 * np.log10(term))**-2

def global_heat_transfer_coefficient_HX2(P, A_tot, D_ext, D_int, k, Re, p, T, f, Circuit):
    ' Calcola U e il delta T logaritmico per lo scambiatore'

    # Proprietà del fluido interno (da CoolProp)
    Pr = CP.PropsSI('Prandtl', 'P', p, 'T', T + 273.15, 'Water')
    k_fluid = CP.PropsSI('L', 'P', p, 'T', T + 273.15, 'Water')

    # Convezione esterna
    if Circuit == 'HX1':
        h_ext = 1000 # W/m^2K # da aggiornare con la formula del PSC del shell side
    elif Circuit == 'HX2':
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
    deltaT_lm = P / (A_tot * U)
    return U, deltaT_lm

def temperature_ISC(T, p, T_sat, deltaT_lm, m, P_term):
    ' Calcola la temperatura nel ramo caldo e nel ramo freddo '
    deltaT = P_term / (m*CP.PropsSI('C', 'T', T+273.15, 'P', p, 'Water'))

    T_h = (T_sat - (deltaT+T_sat)*np.exp(deltaT/deltaT_lm)) / (1 - np.exp(deltaT/deltaT_lm))
    T_c = T_h - deltaT
    return T_c, T_h

def pressure_drop_buoyancy(p, T_h, T_c, H, g=9.81):
    ' Calcola la perdita di carico per differenza di densità'
    rho_h = CP.PropsSI('D', 'T', T_h+273.15, 'P', p, 'Water')
    rho_c = CP.PropsSI('D', 'T', T_c+273.15, 'P', p, 'Water')
    deltaP_buoyancy = (rho_c - rho_h) * g * H
    return deltaP_buoyancy

def pressure_drop_friction(m, N_HX, p, D_pipe, D_HX, T_av, T_c, T_h, eps_pipe, eps_HX, A_pipe, A_HX, L_ISC, L_HX, N_bends, k_bends): 
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
    dp_loc_bends_c = 0.5*N_bends * k_bends * 1/(2*CP.PropsSI('D', 'T', T_c+273.15, 'P', p, 'Water') * A_pipe**2)
    dp_loc_bends_h = 0.5*N_bends * k_bends * 1/(2*CP.PropsSI('D', 'T', T_h+273.15, 'P', p, 'Water') * A_pipe**2)
    # entrata e uscita 
    sigma = (A_HX * N_HX) / A_pipe
    k_in = 0.5*(1-sigma)
    k_out = (1- sigma)**2
    dp_loc_HX = (k_in + k_out) * 1/(2*CP.PropsSI('D', 'T', T_av+273.15, 'P', p, 'Water') * (A_HX*N_HX)**2)

    dp_dict = {
        'dist_cold': dp_dist_c,
        'dist_hot': dp_dist_h,
        'dist_HX': dp_dist_HX,
        'loc_bends_hot': dp_loc_bends_h,
        'loc_bends_cold': dp_loc_bends_c,
        'loc_HX': dp_loc_HX
    }
    deltaP_friction = dp_dist_c + dp_dist_h + dp_dist_HX + dp_loc_bends_h + dp_loc_bends_c + dp_loc_HX
    return deltaP_friction, dp_dict

def mass_flow_rate(deltaP_buoyancy, deltaP_friction):
    ' Calcola la portata massica per una potenza termica data'
    return np.sqrt(deltaP_buoyancy / deltaP_friction)

def iteration(config, m_init=100, T_av_init=120, tolerance=1e-6, max_iter=100):
    """
    Itera il calcolo della portata aggiornando T_av come media tra T_h e T_c.
    
    Args:
        config (dict): Dizionario con i parametri del circuito contenente:
            - 'D_internal': diametro interno tubo [m]
            - 'pressure': pressione [Pa]
            - 'T_HX': temperatura lato scambiatore [°C]
            - 'T_sat': temperatura di saturazione nel circuito [°C]
            - 'eps_pipe': rugosità relativa tubo
            - 'A_pipe': area tubo [m²]
            - 'od_hx': diametro esterno HX [m]
            - 'id_hx': diametro interno HX [m]
            - 'k_hx': conduttività termica HX [W/m·K]
            - 'N_tubes': numero di tubi
            - 'A_hx_tot': area totale HX [m²]
            - 'L_circuit': lunghezza circuito [m]
            - 'L_hx': lunghezza HX [m]
            - 'H': altezza del circuito [m]
            - 'N_bends': numero di gomiti
            - 'k_bends': coefficiente di perdita gomiti
            - 'P_term': potenza termica [W]
            - 'Circuit': tipo di circuito ('HX1', 'HX2', ecc.)
        m_init (float): portata iniziale [kg/s]
        T_av_init (float): temperatura media iniziale [°C]
        tolerance (float): tolleranza di convergenza
        max_iter (int): numero massimo di iterazioni
    """
    m = m_init
    T_av = T_av_init
    
    for iter_num in range(max_iter):
        # Calcolo del numero di Reynolds e fattore di attrito per lo scambiatore
        Re_av = reynolds_number(m, config['id_hx'], config['pressure'], T_av)
        f_av = friction_factor(config['eps_hx'], config['id_hx'], Re_av)
        
        # Calcolo U e deltaT_lm nello scambiatore
        U, deltaT_lm = global_heat_transfer_coefficient_HX2(
            config['P_term'], config['A_hx_tot'], config['od_hx'], config['id_hx'], 
            config['k_hx'], Re_av, config['pressure'], T_av, f_av, Circuit=config['Circuit']
        )
        
        # Calcolo temperature caldo e freddo
        T_c, T_h = temperature_ISC(T_av, config['pressure'], config['T_HX'], deltaT_lm, m, config['P_term'])
        
            
        # Calcolo perdite di carico per galleggiabilità
        deltaP_buoyancy = pressure_drop_buoyancy(config['pressure'], T_h, T_c, config['H'])
        
        # Calcolo perdite di carico per attrito nel circuito completo
        deltaP_friction, dp_dict = pressure_drop_friction(
            m, config['N_tubes'], config['pressure'], config['D_internal'], config['id_hx'], 
            T_av, T_c, T_h, config['eps_pipe'], config['eps_hx'], config['A_pipe'], 
            config['A_hx'], config['L_circuit'], config['L_hx'], config['N_bends'], config['k_bends']
        )
        
        # Calcolo nuova portata
        m_new = mass_flow_rate(deltaP_buoyancy, deltaP_friction)
        T_av_new = (T_h + T_c) / 2

        # Verifica di convergenza
        error = abs(m_new - m) / m
        if error < tolerance:
            print(f"\n✓ Convergenza raggiunta in {iter_num+1} iterazioni ({config['Circuit']})\n")
            return m_new, T_av_new, T_h, T_c, deltaP_buoyancy, dp_dict
        
        m = m_new
        T_av = T_av_new
        
        if (iter_num + 1) % 10 == 0 or iter_num == 0:
            print(f"Iterazione {iter_num+1} ({config['Circuit']}): m = {m:.2f} kg/s, T_av = {T_av:.2f} °C, errore = {error:.4f}")
    
    print(f"\n⚠ Avviso: max iterazioni ({max_iter}) raggiunto senza convergenza ({config['Circuit']})\n")
    return m, T_av, T_h, T_c, deltaP_buoyancy, dp_dict


def save_results_to_csv(dp_dict, m, buoyancy, filename):
    dist_total = (dp_dict['dist_cold'] + dp_dict['dist_hot'] + dp_dict['dist_HX']) * m**2
    loc_total = (dp_dict['loc_bends_hot'] + dp_dict['loc_bends_cold'] + dp_dict['loc_HX']) * m**2
    
    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["VOCE", "VALORE [Pa]"])
        writer.writerow(["--- DISTRIBUTED PRESSURE DROPS ---", ""])
        for k in ['dist_cold', 'dist_hot', 'dist_HX']:
            writer.writerow([k.replace('_', ' ').title(), round(dp_dict[k] * m**2, 4)])
        writer.writerow(["Total Distributed", round(dist_total, 4)])
        
        writer.writerow(["", ""])
        writer.writerow(["--- LOCALIZED PRESSURE DROPS ---", ""])
        for k in ['loc_bends_hot', 'loc_bends_cold', 'loc_HX']:
            writer.writerow([k.replace('_', ' ').title(), round(dp_dict[k] * m**2, 4)])
        writer.writerow(["Total Localized", round(loc_total, 4)])
        
        writer.writerow(["", ""])
        writer.writerow(["--- FINALE RESUME ---", ""])
        writer.writerow(["Driving Force (Buoyancy)", round(buoyancy, 4)])
        writer.writerow(["Total Friction Losses", round(dist_total + loc_total, 4)])
        writer.writerow(["Mass Flow Rate [kg/s]", round(m, 4)])


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
T_sat = CP.PropsSI('T', 'P', 1e5, 'Q', 0, 'Water') - 273.15

# CONFIGURAZIONE ISC
config_ISC = {
    'D_internal': id_ISC,
    'pressure': p_ISC,
    'T_sat': T_sat_ISC,
    'T_HX': T_sat,
    'eps_pipe': eps_ISC,
    'eps_hx': eps_HX2,
    'A_pipe': A_ISC,
    'A_hx': A_HX2,
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
    'Circuit': 'HX2'
}

# --- PSC ---
# Pipe


# HX1


# CONFIGURAZIONE PSC


# ---------- EXECUTION ----------
# ITERAZIONI per ISC
m_res, T_av_res, T_h_res, T_c_res, dp_b_res, dp_dict_res = iteration(config_ISC)

# Stampa risultati finali organizzati
print("\n" + "="*80)
print("RISULTATI FINALI ITERAZIONE".center(80))
print("="*80)
print("\n📊 PARAMETRI TERMICI:")
print(f"{'Portata massica':<35} {m_res:>40.2f} kg/s")
print(f"{'Temperatura media (T_av)':<35} {T_av_res:>40.2f} °C")
print(f"{'Temperatura ramo caldo (T_h)':<35} {T_h_res:>40.2f} °C")
print(f"{'Temperatura ramo freddo (T_c)':<35} {T_c_res:>40.2f} °C")
print(f"{'Temperatura saturazione (limite)':<35} {T_sat_ISC:>40.2f} °C")
print(f"{'Distributed pressure drop':<35} {(dp_dict_res['dist_cold'] + dp_dict_res['dist_hot'] + dp_dict_res['dist_HX']) * m_res**2:>40.2f} Pa")
print(f"{'Localized pressure drop':<35} {(dp_dict_res['loc_bends_hot'] + dp_dict_res['loc_bends_cold'] + dp_dict_res['loc_HX']) * m_res**2:>40.2f} Pa")
print(f"{'Driving force (Buoyancy)':<35} {dp_b_res:>40.2f} Pa")
# Stampa tabella perdite di carico

m_res, T_av_res, T_h_res, T_c_res, dp_b_res, dp_dict_res = iteration(config_ISC)
save_results_to_csv(dp_dict_res, m_res, dp_b_res, filename=os.path.join(os.path.dirname(__file__), "result_ISC.csv"))

# ITERAZIONI per PSC