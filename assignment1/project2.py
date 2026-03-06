import CoolProp.CoolProp as CP
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

def get_pipe_geometry(outside_diameter_inch, thickness_inch):
    """Calcola diametro interno [m] e area [m2] da pollici."""
    d_inner_m = (outside_diameter_inch - 2 * thickness_inch) * 0.0254
    area_m2 = np.pi * (d_inner_m / 2)**2
    return d_inner_m, area_m2

# calcolo preliminare ISC
def global_heat_transfer_coefficient(P, A_tot, D_ext, D_int, k, Re, p, T, f):
    ' Calcola U per lo scambiatore nel ISC'
    # convezione esterna
    flux = P / A_tot
    h_pool = flux**(1-1/3.86) * 2.257**(1/3.86)

    # conduzione
    cond = D_ext*np.log(D_ext/D_int) / (2*k)

    # convezione interna: Correlazione di gnieliski
    Pr = CP.PropsSI('Prandtl', 'P', p, 'T', T+273.15, 'Water')
    if Re < 2300:
        Nu = 3.66
    else:
        Nu = (f/8)*(Re-1000)*Pr / (1 + 12.7*(f/8)**0.5*(Pr**(2/3)-1))
    h_internal = Nu * k / D_int
    
    U = 1 / (1/h_pool + cond + 1/h_internal)

    deltaT_lm = P / (A_tot*U)
    return U, deltaT_lm

def reynolds_number(m, D_int, p, T):
    ' Calcola il numero di Reynolds per lo scambiatore nel ISC'
    mu = CP.PropsSI('V', 'T', T+273.15, 'P', p, 'Water') # viscosità dinamica
    Re = 4*m / (np.pi*D_int*mu)
    return Re

def friction_factor(eps, D, Re):
    ' Calcola il fattore di attrito per lo scambiatore nel ISC'
    term = ((eps/D)/3.7)**1.11 + 6.9/Re
    return (-1.8 * np.log10(term))**-2


def temperature_ISC(T, p, T_sat, deltaT_lm, U, P_term):
    ' Calcola la temperatura nel ramo caldo e nel ramo freddo nel ISC'
    deltaT = P_term / (U*CP.PropsSI('C', 'T', T+273.15, 'P', p, 'Water'))

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
    Re_pipe_h = reynolds_number(m, D_pipe, p, T_h)
    Re_pipe_c = reynolds_number(m, D_pipe, p, T_c)
    f_pipes_h = friction_factor(eps_pipe, D_pipe, Re_pipe_h)
    f_pipes_c = friction_factor(eps_pipe, D_pipe, Re_pipe_c)
    dp_dist_c = f_pipes_c * L_ISC / D_pipe * 1/(2*CP.PropsSI('D', 'T', T_c+273.15, 'P', p, 'Water') * A_pipe**2)
    dp_dist_h = f_pipes_h * L_ISC / D_pipe * 1/(2*CP.PropsSI('D', 'T', T_h+273.15, 'P', p, 'Water') * A_pipe**2)
    # perdite distribuite lungo lo scambiatore
    Re_HX = reynolds_number(m, D_HX, p, T_av)
    f_HX = friction_factor(eps_HX, D_HX, Re_HX)
    dp_dist_HX = f_HX * L_HX / D_HX * 1/(2*CP.PropsSI('D', 'T', T_av+273.15, 'P', p, 'Water') * (A_HX*N_HX)**2)
    # perdite localizzate
    # gomiti
    dp_loc_bends_c = 0.5*N_bends * k_bends * 1/(2*CP.PropsSI('D', 'T', T_c+273.15, 'P', p, 'Water') * A_pipe**2)
    dp_loc_bends_h = 0.5*N_bends * k_bends * 1/(2*CP.PropsSI('D', 'T', T_h+273.15, 'P', p, 'Water') * A_pipe**2)
    # entrata e uscita
    k_in = 0.5*(1-(A_HX*N_HX)/A_pipe)
    k_out = (1- A_pipe/(A_HX*N_HX))**2
    dp_loc_HX = (k_in + k_out) * 1/(2*CP.PropsSI('D', 'T', T_av+273.15, 'P', p, 'Water') * (A_HX*N_HX)**2)

    deltaP_friction = dp_dist_c + dp_dist_h + dp_dist_HX + dp_loc_bends_h + dp_loc_bends_c + dp_loc_HX
    return deltaP_friction

def mass_flow_rate(deltaP_buoyancy, deltaP_friction):
    ' Calcola la portata massica per una potenza termica data'
    return np.sqrt(deltaP_buoyancy / deltaP_friction)

def iteration(m_init=100, T_av_init=120, tolerance=1e-6, max_iter=100):
    """Itera il calcolo della portata aggiornando T_av come media tra T_h e T_c."""
    m = m_init
    T_av = T_av_init
    
    for iter_num in range(max_iter):
        # Calcolo del numero di Reynolds e fattore di attrito per il circuito ISC
        Re_av = reynolds_number(m, id_ISC, p_ISC, T_av)
        f_av = friction_factor(eps_ISC, id_ISC, Re_av)
        
        # Calcolo U e deltaT_lm nello scambiatore
        U, deltaT_lm = global_heat_transfer_coefficient(P_term, A_HX2_tot, od_HX2, id_HX2, k_HX2, Re_av, p_ISC, T_av, f_av)
        
        # Calcolo temperature caldo e freddo
        T_c, T_h = temperature_ISC(T_av, p_ISC, T_sat, deltaT_lm, U, P_term)
        
        # Aggiornamento T_av come media aritmetica
        T_av_new = (T_h + T_c) / 2
        
        # Calcolo perdite di carico per galleggiabilità
        deltaP_buoyancy = pressure_drop_buoyancy(p_ISC, T_h, T_c, H_ISC)
        
        # Calcolo perdite di carico per attrito nel circuito completo
        
        deltaP_friction = pressure_drop_friction(
            m, N_tubes_HX2, p_ISC, id_ISC, id_HX2, T_av_new, T_c, T_h, 
            eps_ISC, eps_HX2, A_ISC, A_HX2, L_ISC, L_HX2, N_bends_ISC, k_bends_ISC
        )
        
        # Calcolo nuova portata
        m_new = mass_flow_rate(deltaP_buoyancy, deltaP_friction)
        
        # Verifica di convergenza
        error = abs(m_new - m)
        if error < tolerance:
            print(f"Convergenza raggiunta in {iter_num+1} iterazioni")
            print(f"m finale = {m_new:.2f} kg/s, T_av finale = {T_av_new:.2f} °C")
            print(f"T_h = {T_h:.2f} °C, T_c = {T_c:.2f} °C")
            return m_new, T_av_new, T_h, T_c
        
        m = m_new
        T_av = T_av_new
        
        if (iter_num + 1) % 10 == 0:
            print(f"Iterazione {iter_num+1}: m = {m:.2f} kg/s, T_av = {T_av:.2f} °C, errore = {error:.4f}")
    
    print(f"Avviso: max iterazioni ({max_iter}) raggiunto senza convergenza")
    return m, T_av, T_h, T_c


# Dati
# Dati design
P_nom = 600e6 #W
P_term = 0.01*P_nom #W

# Dati ISC
L_ISC = 40 #m
H_ISC = 10 #m
p_ISC = 70e5 #Pa
od_ISC = 16 #inch
thickness_ISC = 0.5 #inch
id_ISC, A_ISC = get_pipe_geometry(od_ISC, thickness_ISC) #m
eps_ISC = 2e-4 * id_ISC #m

k_bends_ISC = 0.45 # concentrated loss
N_bends_ISC = 2 # number of bends in the ISC

# Dati HX2
L_HX2 = 7 #m
od_HX2 = 25.4e-3 #m
thickness_HX2 = 1.24e-3 #m
id_HX2 = od_HX2 - 2*thickness_HX2 #m
eps_HX2 = 1e-4 * id_HX2 #m
A_HX2 = np.pi * id_HX2 * L_HX2 #m^2 per tubo
A_HX2_tot = 430 #m^2

p_pool = 1e5 #Pa
T_sat = CP.PropsSI('T', 'P', p_pool, 'Q', 0, 'Water') - 273.15 #°C

k_HX2 = 15 #W/mK
N_tubes_HX2 = 770

m_final, T_av_final, T_h, T_c = iteration()