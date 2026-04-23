import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP

# DATI DA TABELLA
P_nom=3400e6 #W
p_sys=15.51e6 #Pa
n_rods=157*264
m_flow_tot=1.4559e4 #kg/s

F_q=2.6 #heat flux hot channel factor

m_flow_eff=  13456 #kg/s (portata effettiva considerando 5,9% di bypass flow)
A_flow_eff = 3.883 #m^2 (area di flusso effettiva, considerando il 92.5% dell'area totale)



# DATI GEOMETRICI (square array)
Dout_clad=9.5e-3 #m
H_active=4.2672 #m
w=12.6e-3 #m (pitch, passo tra le barre)
s_guaina=0.57e-3 #m (spessore guaina)
D_in_clad=Dout_clad-2*s_guaina #m (diametro interno)
D_pellet=8.22e-3 #m (diametro pellet)

# discretizzazione di z (assiale)
z=np.linspace(-H_active/2, H_active/2, 100) #m (asse z, da -H/2 a H/2)


## 1) Volumetric heat generation rate
Tot_power = P_nom*0.974 #W (potenza totale generata, considerando il 97.4% della potenza nominale)
q_avg=Tot_power / (n_rods * np.pi * (D_in_clad/2)**2 * H_active) #W/m^3 (tasso di generazione di calore volumetrico medio)
q_v_max=q_avg * F_q #W/m^3 (tasso di generazione di calore volumetrico massimo)
lambda_tr=0.0029 # m (lunghezza di trasporto)
D_c=lambda_tr/3 # diffusion coeff. in the core
D_r=0.16 # diffusion coeff. in the reflector
L_r=2.85 # m (lunghezza del riflettore)
delta=D_c/D_r*L_r # transport lenght
H_e=H_active+1.42*lambda_tr+2*delta # m altezza estrapolata
qv_profile=q_v_max*np.cos(np.pi*z/H_e) #W/m^3 (profilo di generazione di calore volumetrico lungo l'asse z)

## 2) Average mass velocity
G_avg = m_flow_eff / A_flow_eff # kg/(m^2 s) (velocità di massa media)
print(f"Average mass velocity G_avg: {G_avg:.2f} kg/(m^2 s)")

## 3) Coolant specific enthalpy profile
A_c = w**2-np.pi/4*Dout_clad**2 # m^2 (area di flusso per canale)
W_hc = G_avg*A_c # kg/s (portata per canale)
A_fuel = np.pi/4*D_pellet**2 # m^2 (area del combustibile)

T_in = 279.44 + 273.15 # K (temperatura di ingresso del refrigerante da tabella)
h_in = CP.PropsSI('H', 'T', T_in, 'P', p_sys, 'Water') # J/kg (entalpia specifica di ingresso)

h_profile  = h_in +1.0267*(q_v_max*A_fuel*H_e)/(W_hc*np.pi)*(np.sin(np.pi*z/H_e)+np.sin(np.pi*H_active/2/H_e)) # J/kg (profilo di entalpia specifica lungo l'asse z)

## 4) Coolant Temperature profile
T_profile = CP.PropsSI('T', 'H', h_profile, 'P', p_sys, 'Water') - 273.15 # °C (profilo di temperatura del refrigerante lungo l'asse z)

plt.figure(figsize=(10, 6))
plt.plot(qv_profile, z)
plt.title('Volumetric Heat Generation Rate along the z-axis')
plt.ylabel('z (m)')
plt.xlabel('q_v (W/m^3)')
plt.grid()


plt.figure(figsize=(10, 6))
plt.plot(h_profile, z)
plt.title('Coolant Specific Enthalpy Profile along the z-axis')
plt.ylabel('z (m)')
plt.xlabel('h (J/kg)')
plt.grid()


plt.figure(figsize=(10, 6))
plt.plot(T_profile, z)
plt.title('Coolant Temperature Profile along the z-axis')
plt.ylabel('z (m)')
plt.xlabel('T (°C)')
plt.grid()
plt.show()