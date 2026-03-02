import CoolProp.CoolProp as CP
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np



def get_pipe_geometry(outside_diameter_inch, thickness_inch):
    """
    Calcola diametro interno e area a partire dalle dimensioni nominali in pollici.
    Restituisce i valori in metri e metri quadri.
    """
    # 1. Calcolo del diametro interno in pollici (OD - 2*t)
    # 2. Conversione immediata in metri (* 0.0254)
    d_inner_m = (outside_diameter_inch - 2 * thickness_inch) * 0.0254
    
    # 3. Calcolo dell'area della sezione trasversale (A = pi * r^2)
    area_m2 = np.pi * (d_inner_m / 2)**2
    
    return d_inner_m, area_m2

def calculate_massflow_and_velocity(Q_total, p, D):
    # Calcola la portata massica e le velocità del liquido e del vapore
    # Entalpia del liquido saturo (inizio ebollizione)
    h_liq = CP.PropsSI('H', 'P', p, 'Q', 0, 'Water')

    # Entalpia del vapore saturo (fine ebollizione)
    h_vap = CP.PropsSI('H', 'P', p, 'Q', 1, 'Water')

    # Salto entalpico di evaporazione (Calore latente)
    delta_h_evap = h_vap - h_liq

    # Portata massica necessaria per far evaporare tutto il fluido
    m_dot = Q_total / delta_h_evap

    # Area della sezione trasversale del tubo
    A = np.pi * (D / 2) ** 2

    # Velocità del liquido
    v_liq = m_dot / (CP.PropsSI('D', 'P', p, 'Q', 0, 'Water') * A)

    # Velocità del vapore
    v_vap = m_dot / (CP.PropsSI('D', 'P', p, 'Q', 1, 'Water') * A)

    return m_dot, v_liq, v_vap





# datas
p = 70 *10**5  # Pa
Q_thermal = 34.8 * 10**6  # W
L1 = 10  # m
# concentrated pressure lost
k_h = 40
k_c = 20
# choose the diameter from the table by doing
# Outside diameter - thickness (STD)
Outside_diameter = 1.9  # inches
thickness = 0.145  # inches  
D = calculate_diameter(Outside_diameter, thickness)

# Esempio di utilizzo della funzione
m_dot, v_liq, v_vap = calculate_massflow_and_velocity(Q_thermal, p, D)
