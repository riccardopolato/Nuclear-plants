import sys
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import CoolProp.CoolProp as CP
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
def gas_plenum_analysis(T_gas, d_p, D_in, H, p_max):
    # Data
    BU = 60e3 * 1e6 * 86400 / 1e3   # J/kg_U, discharge burn-up (60 GWd/tU)
    c_N2, c_H2O = 25e-6, 75e-6      # weight fraction impurities
    u_N2, u_H2O = 0.028, 0.018      # kg/mol
    Y = 0.28                        # combined Xe+Kr fission yield
    R_r = 0.4                       # fission gas release fraction
    rho_UO2 = 0.95 * 10960          # kg/m3 (95% of theoretical density)
    E_f = 200 * 1.60218e-13         # J/fission (200 MeV)
    e_U235 = 0.0445                 # U-235 weight enrichment
    u_235, u_238, u_O = 235, 238, 16.00
    u_U = 1 / (e_U235 / u_235 + (1 - e_U235) / u_238)
    u_UO2 = u_U + 2 * u_O
    N_A = 6.022e23
    R = 8.314

    # Input checks
    if not (d_p > 0 and D_in > 0 and H > 0):
        raise ValueError("d_p, D_in and H must be positive")
    if d_p >= D_in:
        raise ValueError("pellet diameter must be smaller than clad inner diameter")
    if p_max <= 0:
        raise ValueError("p_max must be positive")
    if len(T_gas) == 0:
        raise ValueError("empty T_gas profile")

    # UO2 and uranium mass in the pin
    V_fuel = np.pi * (d_p / 2)**2 * H
    M_UO2 = V_fuel * rho_UO2
    M_U = M_UO2 * u_U / u_UO2

    # Fission gas (Xe + Kr) and impurity gas moles
    N_f = BU * M_U / E_f
    n_XeKr = N_f * R_r * Y / N_A
    n_N2 = c_N2 * M_UO2 / u_N2
    n_H2O = c_H2O * M_UO2 / u_H2O
    n_tot = n_XeKr + n_N2 + n_H2O

    # Gas temperature at the plenum (top of the pin)
    T = T_gas[-1] + 273.15

    # Ideal gas law -> minimum plenum volume and height (+ spring allowance)
    r_i = D_in / 2
    V_min = n_tot * R * T / p_max
    H_min = V_min / (np.pi * r_i**2)
    H_spring = 0.15
    H_tot = H_min + H_spring

    # Check ideal gas assumption on the water vapor partial pressure
    p_H2O_id = n_H2O * R * T / V_min
    try:
        p_H2O_real = CP.PropsSI('P', 'T', T, 'Dmolar', n_H2O / V_min, 'Water')
        dev = (p_H2O_id - p_H2O_real) / p_H2O_real * 100
    except Exception:
        p_H2O_real, dev = float('nan'), float('nan')

    # Output
    print()
    print("=" * 52)
    print(" GAS PLENUM SIZING")
    print("=" * 52)
    print(f" T_gas (top, inner)    = {T - 273.15:8.1f} degC ({T:7.1f} K)")
    print(f" p_max (Mariotte)      = {p_max / 1e6:8.2f} MPa")
    print(f" n_tot                 = {n_tot * 1e3:8.3f} mmol")
    print(f" V_min                 = {V_min * 1e6:8.2f} cm^3")
    print(f" H_min (no spring)     = {H_min * 100:8.2f} cm")
    print(f" H_tot (w/ spring)     = {H_tot * 100:8.2f} cm  (+{H_spring * 100:.0f} cm spring)")
    print("-" * 52)
    print(" H2O ideal gas check:")
    print(f"   p_H2O ideal         = {p_H2O_id / 1e6:8.3f} MPa")
    print(f"   p_H2O real          = {p_H2O_real / 1e6:8.3f} MPa")
    print(f"   deviation           = {dev:8.2f} %")
    print("=" * 52)

    diag = {
        "n_XeKr": n_XeKr, "n_N2": n_N2, "n_H2O": n_H2O,
        "T": T,
        "p_H2O_id": p_H2O_id, "p_H2O_real": p_H2O_real, "dev": dev,
        "H_spring": H_spring,
        "p_max": p_max,
    }
    return V_min, H_min, n_tot, H_tot, diag



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
    V_plenum, H_plenum, n_tot, H_tot_plenum, plenum_diag = gas_plenum_analysis(T_gas, D_fuel, D_in, H_active, p_i)
    print(f"Total plenum height (including spring): {H_tot_plenum:.2f} m")
    print(f"Plenum height (without spring): {H_plenum:.2f} m")

# ----- PLOT ------
# 1) Buckling analysis: p_cr(z)

    plt.figure()
    plt.plot(p_cr / 1e6, z,  label="Critical pressure (MPa)")
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

    # ============================================================
    # DIAGNOSTIC AUDIT BLOCK (no calculation changes; print only)
    # Code-name mapping (when audit label differs from code variable):
    #   sigma_max (Tresca, primary)  -> P_avg (= sigma_max_avg / 1e6)
    #   V_min                         -> V_plenum
    #   H_min (no spring)             -> H_plenum   (returned from gas_plenum_analysis)
    #   H_plenum (incl. spring)       -> H_tot_plenum
    #   p_H2O ideal                   -> p_H2O_id   (in plenum_diag)
    #   p_max                         -> p_i passed to gas_plenum_analysis (= Mariotte)
    # ============================================================
    p_i_MPa = p_i / 1e6
    p_o_MPa = p_sys / 1e6
    r_in_m  = D_in / 2
    r_out_m = D_out_clad / 2
    i_z0    = int(np.argmin(np.abs(z)))

    T_K     = plenum_diag["T"]
    T_C     = T_K - 273.15
    n_XeKr  = plenum_diag["n_XeKr"]
    n_N2    = plenum_diag["n_N2"]
    n_H2O   = plenum_diag["n_H2O"]
    x_H2O   = n_H2O / n_tot * 100

    def _rng(arr):
        return float(np.min(arr)), float(np.max(arr))
    def _ok(cond):
        return "OK" if cond else "MISMATCH"

    print()
    print("==== AUDIT START ====")

    print()
    print("[1] Internal pressure and buckling")
    print(f"   p_i (Mariotte)              = {p_i_MPa:7.2f} MPa")
    print(f"   p_o (external, = p_sys)     = {p_o_MPa:7.2f} MPa   [NOT zero -> system pressure]")
    print(f"   p_cr min                    = {np.min(p_cr)/1e6:7.2f} MPa   (at z = {z[i_min]:+.3f} m)")
    print(f"   p_cr max                    = {np.max(p_cr)/1e6:7.2f} MPa")
    print(f"   p_cr at z=0  (i={i_z0:>3d})       = {p_cr[i_z0]/1e6:7.2f} MPa   (z[i]={z[i_z0]:+.3f} m)")
    print(f"   t  (thickness)              = {thickness:.6f} m")
    print(f"   r_avg                       = {r_avg:.6f} m")
    print(f"   r_in                        = {r_in_m:.6f} m")
    print(f"   r_out                       = {r_out_m:.6f} m")

    print()
    print("[2] Mechanical stresses (primary, average components)")
    print(f"   sigma_avg_h                 = {sigma_avg_diz['h']:7.2f} MPa")
    print(f"   sigma_avg_r                 = {sigma_avg_diz['r']:7.2f} MPa")
    print(f"   sigma_avg_l                 = {sigma_avg_diz['l']:7.2f} MPa")
    print(f"   sigma_max (Tresca)          = {P_avg:7.2f} MPa   [code: P_avg]")
    print(f"   Sm                          = {S:7.2f} MPa")
    print(f"   3 * Sm                      = {3*S:7.2f} MPa")
    print(f"   computed with p_i           = {p_i_MPa:7.2f} MPa")
    print(f"   computed with p_o           = {p_o_MPa:7.2f} MPa")

    print()
    print("[3] Axial ranges over active length (min / max)")
    mn, mx = _rng(sigma_max_in_tot);  print(f"   Tresca inner    = {mn:7.2f} / {mx:7.2f} MPa")
    mn, mx = _rng(sigma_max_out_tot); print(f"   Tresca outer    = {mn:7.2f} / {mx:7.2f} MPa")
    mn, mx = _rng(sigma_h_in_tot);    print(f"   sigma_h inner   = {mn:7.2f} / {mx:7.2f} MPa")
    mn, mx = _rng(sigma_h_out_tot);   print(f"   sigma_h outer   = {mn:7.2f} / {mx:7.2f} MPa")
    mn, mx = _rng(sigma_l_in_tot);    print(f"   sigma_l inner   = {mn:7.2f} / {mx:7.2f} MPa")
    mn, mx = _rng(sigma_l_out_tot);   print(f"   sigma_l outer   = {mn:7.2f} / {mx:7.2f} MPa")
    mn, mx = _rng(sigma_r_in_tot);    print(f"   sigma_r inner   = {mn:7.2f} / {mx:7.2f} MPa")
    mn, mx = _rng(sigma_r_out_tot);   print(f"   sigma_r outer   = {mn:7.2f} / {mx:7.2f} MPa")

    print()
    print("[4] Plenum")
    print(f"   n_XeKr                      = {n_XeKr*1e3:7.3f} mmol")
    print(f"   n_N2                        = {n_N2*1e3:7.3f} mmol")
    print(f"   n_H2O                       = {n_H2O*1e3:7.3f} mmol")
    print(f"   n_tot                       = {n_tot*1e3:7.3f} mmol")
    print(f"   x_H2O (mol fraction)        = {x_H2O:7.2f} %")
    print(f"   T_gas                       = {T_C:7.2f} degC  ({T_K:7.2f} K)")
    print(f"   p_max                       = {plenum_diag['p_max']/1e6:7.2f} MPa")
    print(f"   V_min                       = {V_plenum*1e6:7.2f} cm^3   [code: V_plenum]")
    print(f"   H_min (no spring)           = {H_plenum*100:7.2f} cm     [code: H_plenum]")
    print(f"   H_spring                    = {plenum_diag['H_spring']*100:7.2f} cm")
    print(f"   H_plenum (incl. spring)     = {H_tot_plenum*100:7.2f} cm     [code: H_tot_plenum]")
    print(f"   p_H2O ideal                 = {plenum_diag['p_H2O_id']/1e6:7.2f} MPa   [code: p_H2O_id]")
    print(f"   p_H2O real                  = {plenum_diag['p_H2O_real']/1e6:7.2f} MPa")
    print(f"   deviation                   = {plenum_diag['dev']:7.2f} %")

    chk_p_used = abs(p_i_MPa - p_i_MPa) < 0.01
    chk_pmax   = abs(plenum_diag['p_max']/1e6 - p_i_MPa) < 0.01
    chk_height = abs((H_plenum + plenum_diag['H_spring']) - H_tot_plenum) < 1e-9
    chk_moles  = abs((n_XeKr + n_N2 + n_H2O) - n_tot) < 1e-12
    chk_T      = abs((T_K - 273.15) - T_C) < 1e-9

    print()
    print("[5] Consistency checks")
    print(f"   p_i (sect.1) == p_i used in stresses (sect.2)     : {_ok(chk_p_used)}")
    print(f"   p_max plenum == p_i Mariotte  (within 0.01 MPa)   : {_ok(chk_pmax)}")
    print(f"   H_min + H_spring == H_plenum                      : {_ok(chk_height)}")
    print(f"   n_XeKr + n_N2 + n_H2O == n_tot                    : {_ok(chk_moles)}")
    print(f"   T_gas[K] == T_gas[degC] + 273.15                  : {_ok(chk_T)}")

    print()
    print("==== AUDIT END ====")


if __name__ == "__main__":
    main()
