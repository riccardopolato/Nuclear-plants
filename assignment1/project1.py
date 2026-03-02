import CoolProp.CoolProp as CP
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

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
D, A = calculate_diameter_and_area(Outside_diameter, thickness)

def calculate_diameter_and_area(outside_diameter, thickness):
    # Calcola il diametro in metri e l'area della sezione trasversale del tubo
    # Formula diametro: (Diametro esterno - Spessore) * Fattore di conversione (0.0254)
    D = (outside_diameter - 2 * thickness) * 0.0254
    # Formula area: π * (D / 2)^2
    A = np.pi * (D / 2) ** 2
    return D, A

#ricavo le velocità
# def calculate_massflow_and_velocity(Q_total, p, D):
#     # Calcola la portata massica e le velocità del liquido e del vapore
#     # Entalpia del liquido saturo (inizio ebollizione)
#     h_liq = CP.PropsSI('H', 'P', p, 'Q', 0, 'Water')

#     # Entalpia del vapore saturo (fine ebollizione)
#     h_vap = CP.PropsSI('H', 'P', p, 'Q', 1, 'Water')

#     # Salto entalpico di evaporazione (Calore latente)
#     delta_h_evap = h_vap - h_liq

#     # Portata massica necessaria per far evaporare tutto il fluido
#     m_dot = Q_total / delta_h_evap

#     # Area della sezione trasversale del tubo
#     A = np.pi * (D / 2) ** 2

#     # Velocità del liquido
#     v_liq = m_dot / (CP.PropsSI('D', 'P', p, 'Q', 0, 'Water') * A)

#     # Velocità del vapore
#     v_vap = m_dot / (CP.PropsSI('D', 'P', p, 'Q', 1, 'Water') * A)

#     return m_dot, v_liq, v_vap

# Esempio di utilizzo della funzione
# m_dot, v_liq, v_vap = calculate_massflow_and_velocity(Q_thermal, p, D)
