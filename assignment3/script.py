import sys
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

from inputs_loader import load_inputs

# CSV con i profili di temperatura della guaina, prodotto dall'Assignment 2
# (eseguire prima assignment2.py per (ri)generarlo).
A2_DIR = Path(__file__).resolve().parent.parent / "assignment2"
CLADDING_T_CSV = A2_DIR / "output_cladding_T.csv"

# Cartella di output per i grafici (assignment3/plots)
PLOTS_DIR = Path(__file__).resolve().parent / "plots"


def load_cladding_temperatures(csv_path=CLADDING_T_CSV):
    """
    Legge il CSV prodotto dall'Assignment 2 (colonne: z_m, T_ci_C, T_co_C).
    Ritorna (z, T_ci, T_co): posizione assiale [m] e temperature [deg C]
    della parete interna ed esterna della guaina.
    """
    data = np.genfromtxt(csv_path, delimiter=",", names=True)
    return data["z_m"], data["T_ci_C"], data["T_co_C"]

# properties
def properties(T):
    E = (9.9e3 - 5.669 * (T - 273)) * 9.81   # MPa
    nu = 0.3303 + 8.376e-5 * (T - 273)
    alpha = (6.72e-6 * T - 2.07e-3)/(T - 308)
    return E, nu, alpha

# 1) BUCKLING ANALYSIS
def buckling(T, r_avg, t):
    """
    Pressione critica di collasso (buckling), Eq.(1) del PDF:
        p_cr = E / (4 (1 - nu^2)) * (t / r_avg)^3
    T: temperatura della guaina in K (scalare o vettore); E e nu (Zircaloy-4,
    Kye-Ho Tab.2) sono valutati a quella T. Ritorna p_cr in Pa.
    """

    # proprieta' (T in K)
    E, nu, alpha = properties(T)

    p_cr = E / (4 * (1 - nu**2)) * (t / r_avg)**3   # MPa
    return p_cr * 1e6  # Pa                              


# 2) INTERNAL PRESSURE ANALYSIS
def internal_pressure_analysis(sigma_y, r_avg, t):
    p_i = sigma_y * t / r_avg  # MPa
    return p_i * 1e6  # Pa


# 3) STRESS ANALYSIS
def stress_analysis(r_in, r_out, p_i, p_o):
    sigma_h_in = (p_i*(r_in**2 + r_out**2) - 2*p_o*r_out**2)/(r_out**2 - r_in**2)  # MPa, tensione circonferenziale interna
    sigma_r_in = - p_i
    sigma_l_in = (p_i*r_in**2 - p_o*r_out**2)/(r_out**2 - r_in**2)  # MPa, tensione longitudinale interna

    sigma_h_out = (2*p_i*r_in**2 - p_o*(r_in**2 + r_out**2))/(r_out**2 - r_in**2)  # MPa, tensione circonferenziale interna
    sigma_r_out = - p_o
    sigma_l_out = (p_i*r_in**2 - p_o*r_out**2)/(r_out**2 - r_in**2)  # MPa, tensione longitudinale interna

    sigma_avg_h = (sigma_h_in + sigma_h_out) / 2
    sigma_avg_r = (sigma_r_in + sigma_r_out) / 2   
    sigma_avg_l = (sigma_l_in + sigma_l_out) / 2

    # Meccanici secondari
    sigma_h_in_sec = sigma_h_in - sigma_avg_h  # MPa, tensione circonferenziale secondaria interna
    sigma_r_in_sec = sigma_r_in - sigma_avg_r  # MPa, tensione radiale secondaria interna
    sigma_l_in_sec = sigma_l_in - sigma_avg_l  # MPa, tensione longitudinale secondaria interna

    sigma_h_out_sec = sigma_h_out - sigma_avg_h  # MPa, tensione circonferenziale secondaria esterna
    sigma_r_out_sec = sigma_r_out - sigma_avg_r  # MPa, tensione radiale secondaria esterna
    sigma_l_out_sec = sigma_l_out - sigma_avg_l  # MPa, tensione longitudinale secondaria esterna

    sigma_max_avg = max(abs(sigma_avg_h-sigma_avg_r), abs(sigma_avg_r-sigma_avg_l), abs(sigma_avg_l-sigma_avg_h))  # MPa, tensione massima di Von Mises interna
    
    sigma_in_diz = {"h": sigma_h_in, "r": sigma_r_in, "l": sigma_l_in}
    sigma_out_diz = {"h": sigma_h_out, "r": sigma_r_out, "l": sigma_l_out}
    sigma_in_sec_diz = {"h": sigma_h_in_sec, "r": sigma_r_in_sec, "l": sigma_l_in_sec}
    sigma_out_sec_diz = {"h": sigma_h_out_sec, "r": sigma_r_out_sec, "l": sigma_l_out_sec}
    sigma_avg_diz = {"h": sigma_avg_h, "r": sigma_avg_r, "l": sigma_avg_l}
    
    return sigma_max_avg, sigma_in_diz, sigma_out_diz, sigma_in_sec_diz, sigma_out_sec_diz, sigma_avg_diz

# 4) THERMAL STRESS
def thermal_stress(T_ci, T_co, r_in, r_out):
    # proprieta' (T in K)
    E, nu, alpha = properties((T_ci + T_co)/2)

    term1 = (E*alpha*(T_ci - T_co))/(2*(1 - nu))/np.log(r_out/r_in)  # MPa, termine comune
    sigma_h_in = term1*(1-(2*r_out**2)/(r_out**2-r_in**2)* np.log(r_out/r_in))  # MPa, tensione termica interna
    sigma_h_out = term1*(1-(2*r_in**2)/(r_out**2-r_in**2)* np.log(r_out/r_in))  # MPa, tensione termica interna
    sigma_th_in = {"h": sigma_h_in, "r": np.zeros_like(sigma_h_in), "l": sigma_h_in}
    sigma_th_out = {"h": sigma_h_out, "r": np.zeros_like(sigma_h_out), "l": sigma_h_out}
    return sigma_th_in, sigma_th_out


# 5) GAS PLENUM ANALYSIS
def gas_plenum_analysis(T_gas, D_f, D_in, H, p_i_max):
    # DATI
    # gas prodotti solo da Xe e Kr e possono essere trattati come gas perfetti
    burn_up_MWd = 60e3                            # MWd/tU
    burn_up_J = burn_up_MWd * 1e6 * 86400 / 1e3   # J/kg_U  (1 MWd = 1e6 W * 86400 s; per tU -> per 1e3 kg)
    N2 = 25   # ppm in peso (impurita')
    H2O = 75  # ppm in peso (impurita')
    M_N2, M_H2O = 0.028, 0.018               # kg/mol
    Y_XeKr = 0.28                # resa di fissione COMBINATA Xe+Kr (atomi di gas nobile per fissione)
    R_r=0.4 
    rho_teorica = 10960                   # frazione di gas rilasciata dalla pastiglia
    rho_fuel = 0.95 * rho_teorica     # kg/m3, densita' UO2 (95% del teorico)
    E_fiss = 200 * 1.60218e-13   # J, energia rilasciata per fissione (200 MeV)
    e_U235 = 0.0445              # arricchimento (frazione in peso di U-235); coerente con enr_max in main (AP1000, verificare su ML11171A443)
    M_235, M_238, M_O = 235, 238, 16.00            # g/mol
    M_U = 1 / (e_U235 / M_235 + (1 - e_U235) / M_238)    # g/mol, massa molare media dell'U arricchito (media in peso)
    M_UO2 = M_U + 2 * M_O        # g/mol
    N_A = 6.022e23              # 1/mol, numero di Avogadro
    # massa di combustibile e di uranio nella barretta
    V_fuel = np.pi * (D_f / 2)**2 * H   # m3, volume della colonna di pastiglie
    m_UO2 = V_fuel * rho_fuel           # kg, massa di UO2
    f_uranio = M_U / M_UO2                   # frazione in peso di uranio nell'UO2
    m_U = m_UO2 * M_U / M_UO2           # kg, massa di uranio
    R = 8.314  # J/(mol K), costante dei gas
    # moli di gas di fissione rilasciati
    N_fiss = burn_up_J * m_U / E_fiss   # adimensionale (numero di fissioni)
    n_fiss = N_fiss / N_A                   # mol di fissioni
    n_fg = N_fiss * R_r * Y_XeKr / N_A   # mol di gas di fissione (Xe+Kr) rilasciate


    # moli di impurità
    n_N2  = (N2  * 1e-6 * m_UO2) / M_N2       # mol  (N2  = 25 ppm)
    n_H2O = (H2O * 1e-6 * m_UO2) / M_H2O      # mol  (H2O = 75 ppm)

    # calcolo le moli totali
    n_tot = n_fg + n_N2 + n_H2O

    # temperatura del gas = T interna guaina alla quota del plenum (in cima al pin)
    T = T_gas[-1] + 273.15  # K  (T_gas e' in degC)

    # legge dei gas perfetti -> volume del plenum a p_i_max
    V = n_tot * R * T / p_i_max          # m3
    # il plenum e' delimitato dalla parete INTERNA della guaina (non dal pellet)
    A_plenum = np.pi * (D_in / 2)**2     # m2, sezione interna guaina
    H_plenum = V / A_plenum              # m, altezza del plenum
    H_tot_plenum = H_plenum + 0.15 #m, add 15 cm to account for the spring

    # ---- CHECK: validita' del gas ideale per il vapore d'acqua ----
    p_H2O = n_H2O / n_tot * p_i_max      # Pa, pressione parziale H2O (legge di Dalton)
    Tc_H2O, pc_H2O = 647.1, 22.06e6      # K, Pa (coordinate critiche acqua)
    Tr, pr = T / Tc_H2O, p_H2O / pc_H2O
    try:
        import CoolProp.CoolProp as CP
        Z = CP.PropsSI('Z', 'T', T, 'P', p_H2O, 'Water')   # fattore di compressibilita'
    except Exception:
        Z = float('nan')

    # ---- RISULTATO FINALE ----
    print()
    print("=" * 52)
    print(" GAS PLENUM SIZING")
    print("=" * 52)
    print(f" n (Xe+Kr) rilasciate  = {n_fg * 1e3:8.3f} mmol")
    print(f" n N2                  = {n_N2 * 1e3:8.3f} mmol")
    print(f" n H2O                 = {n_H2O * 1e3:8.3f} mmol")
    print(f" n_tot                 = {n_tot * 1e3:8.3f} mmol")
    print("-" * 52)
    print(f" T_gas (cima, interna) = {T - 273.15:8.1f} degC ({T:7.1f} K)")
    print(f" p_i,max (Mariotte)    = {p_i_max / 1e6:8.2f} MPa")
    print(f" V_plenum              = {V * 1e6:8.2f} cm^3")
    print(f" H_plenum              = {H_plenum * 100:8.2f} cm")
    print("-" * 52)
    gas_ok = "OK (Z~1)" if abs(Z - 1) < 0.05 else "non ideale, verificare"
    print(f" H2O ideal gas: Tr={Tr:.3f}  pr={pr:.3f}  Z={Z:.3f}  -> {gas_ok}")
    print("=" * 52)

    return V, H_plenum, n_tot,H_tot_plenum



# ------- \MAIN -------
def main():
    # importare dati da file csv
    inputs, _ = load_inputs()  # legge inputs.csv in questa stessa cartella

    D_out_clad = inputs["D_out_clad"]    # m, diametro esterno guaina
    thickness  = inputs["s_clad"]        # m, spessore guaina
    D_in = D_out_clad - 2 * thickness    # m, diametro interno guaina
    r_avg = (D_out_clad + D_in) / 4          # raggio medio = (r_in + r_out)/2
    p_sys = inputs["p_sys"]              # Pa, pressione di sistema
    D_fuel = inputs["D_pellet"]            # m, diametro del pellet di combustibile
    H_active = inputs["H_active"]          # m, altezza attiva della barra
    
    # mechanical properties    
    sigma_y = 241  # MPa, tensione di snervamento 
    sigma_u = 413  # MPa, tensione di rottura
    S = min(2/3 * sigma_y, 1/3 * sigma_u)  # MPa, tensione ammissibile (Kye-Ho Tab.2)

    # Temperature della guaina rilette dall'output dell'Assignment 2
    z, T_ci, T_co = load_cladding_temperatures()
    T_gas = T_ci
    # temperatura media nella parete (media interna/esterna)
    T_c_avg = (T_co + T_ci) / 2          # deg C, profilo lungo z

    # 1) BUCKLING ANALYSIS  -- p_cr(z) sul profilo; il caso peggiore e' il minimo
    p_cr = buckling(T_c_avg + 273.15, r_avg, thickness)  # Pa
    i_min = int(np.argmin(p_cr))
    # check (collasso: la p critica deve superare la pressione esterna)
    if p_cr[i_min] > p_sys:
        print("The cladding is safe against buckling.")
    else:
        print("The cladding is NOT safe against buckling.")
    
    # 2) INTERNAL PRESSURE ANALYSIS
    p_i = internal_pressure_analysis(sigma_y, r_avg, thickness)  # Pa

    # 3) STRESS ANALYSIS
    sigma_max_avg, sigma_in_diz, sigma_out_diz,sigma_in_sec_diz,sigma_out_sec_diz,sigma_avg_diz = stress_analysis(D_in/2, D_out_clad/2, p_i, p_sys)  # Pa
    P_avg = sigma_max_avg / 1e6     # MPa  
    # componenti meccaniche in MPa (stress_analysis lavora in Pa: p_i, p_sys sono in Pa)
    sigma_in_diz = {k: v / 1e6 for k, v in sigma_in_diz.items()}    # MPa
    sigma_out_diz = {k: v / 1e6 for k, v in sigma_out_diz.items()}  # MPa
    sigma_avg_diz = {k: v / 1e6 for k, v in sigma_avg_diz.items()}
    sigma_in_sec_diz = {k: v / 1e6 for k, v in sigma_in_sec_diz.items()}  # MPa
    sigma_out_sec_diz = {k: v / 1e6 for k, v in sigma_out_sec_diz.items()}  # MPa


    # Tabella di verifica ASME (stress primari, condizioni di design): sigma_max <= Sm
    print()
    print("=" * 52)
    print(" ASME verification - primary stresses (sigma_max <= Sm)")
    print("=" * 52)
    print(f" Sm = min(2/3 sy, 1/3 su) = {S:.2f} MPa")
    print("-" * 52)
    print(f" {'Wall':<8}{'sigma_max [MPa]':>18}{'Sm [MPa]':>12}{'Check':>10}")
    # ASME check: use average stress value `P_avg` (one value instead of inner/outer)
    name = "Average"
    smax = P_avg
    verdict = "OK" if smax < S else "NOT OK"
    print(f" {name:<8}{smax:>18.2f}{S:>12.2f}{verdict:>10}")
    print("=" * 52)

    # 4) THERMAL STRESS
    sigma_th_in_diz, sigma_th_out_diz = thermal_stress(T_ci + 273.15, T_co + 273.15, D_in/2, D_out_clad/2)  # MPa 
    
    sigma_h_in_tot = sigma_avg_diz["h"]+sigma_in_sec_diz["h"] + sigma_th_in_diz["h"] # primarimeccanici+ secondari meccanici + termici
    sigma_h_out_tot = sigma_avg_diz["h"]+sigma_out_sec_diz["h"] + sigma_th_out_diz["h"]
    sigma_l_in_tot = sigma_avg_diz["l"]+sigma_in_sec_diz["l"] + sigma_th_in_diz["l"]
    sigma_l_out_tot = sigma_avg_diz["l"]+sigma_out_sec_diz["l"] + sigma_th_out_diz["l"]
    sigma_r_in_tot = sigma_avg_diz["r"]+sigma_in_sec_diz["r"] + sigma_th_in_diz["r"]
    sigma_r_out_tot = sigma_avg_diz["r"]+sigma_out_sec_diz["r"] + sigma_th_out_diz["r"]
    
    # Tresca (element-wise) per in e per out: massimo fra le tre differenze in ogni punto z
    arr1_in = np.abs(sigma_l_in_tot - sigma_h_in_tot)
    arr2_in = np.abs(sigma_l_in_tot - sigma_r_in_tot)
    arr3_in = np.abs(sigma_h_in_tot - sigma_r_in_tot)
    sigma_max_in_tot = np.maximum.reduce([arr1_in, arr2_in, arr3_in])

    arr1_out = np.abs(sigma_l_out_tot - sigma_h_out_tot)
    arr2_out = np.abs(sigma_l_out_tot - sigma_r_out_tot)
    arr3_out = np.abs(sigma_h_out_tot - sigma_r_out_tot)
    sigma_max_out_tot = np.maximum.reduce([arr1_out, arr2_out, arr3_out])

    # confronto scalare sul massimo lungo il profilo z (peggiore caso)
    if np.max(sigma_max_in_tot) < 3 * S and np.max(sigma_max_out_tot) < 3 * S:
        print("The cladding is safe against combined mechanical and thermal stresses.")
    else:
        print("The cladding is NOT safe against combined mechanical and thermal stresses.")

    # calcolo delle tensioni totali (meccaniche + termiche) per ogni componente (h, r, l) e per interno/esterno
    

    # 5) GAS PLENUM SIZING  -- p_i_max = pressione massima ammissibile (Mariotte, sigma_h = sigma_y)
    V_plenum, H_plenum, n_tot, H_tot_plenum = gas_plenum_analysis(T_gas, D_fuel, D_in, H_active, p_i)
    print(f"Total plenum height (including spring): {H_tot_plenum:.2f} m")
    print(f"Plenum height (without spring): {H_plenum:.2f} m")

# ----- PLOT ------
# 1) Buckling analysis: p_cr(z)

    plt.figure()
    plt.plot(p_cr / 1e6, z,  label="Critical pressure (MPa)")
    plt.axvline(p_sys / 1e6, color="red", linestyle="--", label="System pressure (MPa)")
    plt.ylabel("z (m)")
    plt.xlabel("Pressure (MPa)")
    plt.title("Buckling analysis")
    plt.legend()
    plt.grid()
    PLOTS_DIR.mkdir(exist_ok=True)
    plt.savefig(PLOTS_DIR / "buckling_analysis_pcr.png", dpi=300)
    
    plt.figure()
    plt.plot(sigma_max_in_tot,  z, label="Max stress inner (MPa)")
    plt.plot(sigma_max_out_tot, z, label="Max stress outer (MPa)")
    plt.axvline(3*S, color="red", linestyle="--", label="Sm (MPa)")
    plt.ylabel("z (m)")
    plt.xlabel("Pressure (MPa)")
    plt.title("Tresca vs 3*Sm")
    plt.legend()
    plt.grid()
    PLOTS_DIR.mkdir(exist_ok=True)
    plt.savefig(PLOTS_DIR / "Mechanical_stresses.png", dpi=300)

    # 2) Tensioni totali (meccaniche + termiche) per componente (h, r, l), pareti interna/esterna
    plt.figure()
    plt.plot(sigma_h_in_tot,  z, label=r"$\sigma_h$ inner")
    plt.plot(sigma_h_out_tot, z, label=r"$\sigma_h$ outer")
    plt.plot(sigma_r_in_tot,  z, label=r"$\sigma_r$ inner")
    plt.plot(sigma_r_out_tot, z, label=r"$\sigma_r$ outer")
    plt.plot(sigma_l_in_tot,  z, label=r"$\sigma_l$ inner")
    plt.plot(sigma_l_out_tot, z, label=r"$\sigma_l$ outer")
    plt.xlabel("Stress (MPa)")
    plt.ylabel("z (m)")
    plt.title("Total stresses (mechanical + thermal)")
    plt.legend()
    plt.grid()
    plt.savefig(PLOTS_DIR / "total_stresses.png", dpi=300)



if __name__ == "__main__":
    main()
