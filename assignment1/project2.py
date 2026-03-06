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

def reynolds_number(m, D_int, mu, p, T):
    ' Calcola il numero di Reynolds per lo scambiatore nel ISC'
    mu = CP.PropsSI('V', 'T', T+273.15, 'P', p) # viscosità dinamica
    Re = 4*m / (np.pi*D_int*mu)
    return Re

def friction_factor(eps, D, Re):
    ' Calcola il fattore di attrito per lo scambiatore nel ISC'
    term = ((eps/D)/3.7)**1.11 + 6.9/Re
    return (-1.8 * np.log10(term))**-2


def temperature_ISC(T, p, T_sat, deltaT_lm):
    ' Calcola la temperatura nel ramo caldo e nel ramo freddo nel ISC'
    deltaT = P_term / (m*CP.PropsSI('C', 'T', T+273.15, 'P', p))

    T_h = (T_sat - (deltaT+T_sat)*np.exp(deltaT/deltaT_lm)) / (1 - np.exp(deltaT/deltaT_lm))
    T_c = T_h - deltaT
    return T_c, T_h

def pressure_drop_buoyancy(p, T_h, T_c, H, g=9.81):
    ' Calcola la perdita di carico per differenza di densità'
    rho_h = CP.PropsSI('D', 'T', T_h+273.15, 'P', p)
    rho_c = CP.PropsSI('D', 'T', T_c+273.15, 'P', p)
    deltaP_buoyancy = (rho_c - rho_h) * g * H
    return deltaP_buoyancy

def pressure_drop_friction(m, D_int, mu, p, T_av, eps_pipe, eps_HX): 
    ' Calcola la perdita di carico per attrito'
    # perdite distibuite lungo i tubi
    Re_pipe = reynolds_number(m, D_int, mu, p, T_av)
    f_pipes = friction_factor(eps_pipe, D_int, Re_pipe)
    dp_dist_pipe = f_pipes * (m**2) / (2*CP.PropsSI('D', 'T', T_av+273.15, 'P', p) * D_int)

    # perdite distribuite lungo lo scambiatore
    Re_HX = reynolds_number(m, D_int, mu, p, T_av)
    f_HX = friction_factor(eps_HX, D_int, Re_HX)

    # perdite localizzate
    # gomiti

    return deltaP_friction




# Dati
# Dati assumed for first iteration
m = 100 #kg/s
T_av = 120 #°C

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

k_loss_ISC = 0.45 # concentrated loss
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