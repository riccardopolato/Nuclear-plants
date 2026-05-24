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

    sigma_max_in = max(abs(sigma_h_in-sigma_r_in), abs(sigma_r_in-sigma_l_in), abs(sigma_l_in-sigma_h_in))  # MPa, tensione massima di Von Mises interna
    sigma_max_out = max(abs(sigma_h_out-sigma_r_out), abs(sigma_r_out-sigma_l_out), abs(sigma_l_out-sigma_h_out))  # MPa, tensione massima di Von Mises esterna

    sigma_in = {"h": sigma_h_in, "r": sigma_r_in, "l": sigma_l_in}
    sigma_out = {"h": sigma_h_out, "r": sigma_r_out, "l": sigma_l_out}
    return sigma_max_in, sigma_max_out, sigma_in, sigma_out

# 4) THERMAL STRESS
def thermal_stress(T_ci, T_co, r_in, r_out):
    # proprieta' (T in K)
    E, nu, alpha = properties((T_ci + T_co)/2)

    term1 = (E*alpha*(T_ci - T_co))/(2*(1 - nu))/np.log(r_out/r_in)  # MPa, termine comune
    sigma_h_in = term1*(1-(2*r_out**2)/(r_out**2-r_in**2)* np.log(r_out/r_in))  # MPa, tensione termica interna
    sigma_h_out = term1*(1-(2*r_in**2)/(r_out**2-r_in**2)* np.log(r_out/r_in))  # MPa, tensione termica interna
    sigma_th_in = {"h": sigma_h_in, "r": np.zeros_like(sigma_h_in), "l": sigma_h_in}
    sigma_th_out = {"h": sigma_h_out, "r": np.zeros_like(sigma_h_out), "l": sigma_h_out}
    return sigma_h_in, sigma_h_out, sigma_th_in, sigma_th_out


# 5) GAS PLENUM ANALYSIS
def gas_plenum_analysis(T_gas, d_p, D_in, H, p_max):
    # Data
    BU_MWd = 60e3
    BU = BU_MWd * 1e6 * 86400 / 1e3
    c_N2 = 25e-6
    c_H2O = 75e-6
    u_N2, u_H2O = 0.028, 0.018
    Y = 0.28
    R_r = 0.4
    rho_th = 10960
    rho_UO2 = 0.95 * rho_th
    E_f = 200 * 1.60218e-13
    e_U235 = 0.0445
    u_235, u_238, u_O = 235, 238, 16.00
    u_U = 1 / (e_U235 / u_235 + (1 - e_U235) / u_238)
    u_UO2 = u_U + 2 * u_O
    N_A = 6.022e23
    R = 8.314

    # UO2 and uranium mass in the pin
    V_fuel = np.pi * (d_p / 2)**2 * H
    M_UO2 = V_fuel * rho_UO2
    M_U = M_UO2 * u_U / u_UO2

    # Fission gas moles (Xe + Kr)
    N_f = BU * M_U / E_f
    n_XeKr = N_f * R_r * Y / N_A

    # Impurity gas moles
    n_N2 = c_N2 * M_UO2 / u_N2
    n_H2O = c_H2O * M_UO2 / u_H2O
    
    # Total moles
    n_tot = n_XeKr + n_N2 + n_H2O

    # Gas temperature at the plenum (top of the pin)
    T = T_gas[-1] + 273.15

    # Ideal gas law -> plenum volume and height (+ spring allowance)
    r_i = D_in / 2
    V_min = n_tot * R * T / p_max
    H_min = V_min / (np.pi * r_i**2)
    H_spring = 0.15
    H_plenum = H_min + H_spring

    # Check ideal gas assumption for water vapor (slide step 5)
    p_H2O_id = n_H2O * R * T / V_min
    import CoolProp.CoolProp as CP
    p_H2O_real = CP.PropsSI('P', 'T', T, 'Dmolar', n_H2O / V_min, 'Water')
    dev = (p_H2O_id - p_H2O_real) / p_H2O_real * 100

    # Output
    print()
    print("=" * 52)
    print(" GAS PLENUM SIZING")
    print("=" * 52)
    print(f" T_gas (top, inner)    = {T - 273.15:8.1f} degC ({T:7.1f} K)")
    print(f" p_max (Mariotte)      = {p_max / 1e6:8.2f} MPa")
    print(f" V_min                 = {V_min * 1e6:8.2f} cm^3")
    print(f" H_min (no spring)     = {H_min * 100:8.2f} cm")
    print(f" H_plenum (w/ spring)  = {H_plenum * 100:8.2f} cm  (+{H_spring * 100:.0f} cm spring)")
    print("-" * 52)
    print(" H2O ideal gas check:")
    print(f"   p_H2O ideal         = {p_H2O_id / 1e6:8.3f} MPa")
    print(f"   p_H2O real          = {p_H2O_real / 1e6:8.3f} MPa")
    print(f"   deviation           = {dev:8.2f} %")
    print("=" * 52)

    return V_min, H_plenum, n_tot



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
    sigma_max_in, sigma_max_out, sigma_in_diz, sigma_out_diz = stress_analysis(D_in/2, D_out_clad/2, p_i, p_sys)  # Pa
    P_in = sigma_max_in / 1e6     # MPa
    P_out = sigma_max_out / 1e6    # MPa
    # componenti meccaniche in MPa (stress_analysis lavora in Pa: p_i, p_sys sono in Pa)
    sigma_in_diz = {k: v / 1e6 for k, v in sigma_in_diz.items()}    # MPa
    sigma_out_diz = {k: v / 1e6 for k, v in sigma_out_diz.items()}  # MPa

    # Tabella di verifica ASME (stress primari, condizioni di design): sigma_max <= Sm
    print()
    print("=" * 52)
    print(" ASME verification - primary stresses (sigma_max <= Sm)")
    print("=" * 52)
    print(f" Sm = min(2/3 sy, 1/3 su) = {S:.2f} MPa")
    print("-" * 52)
    print(f" {'Wall':<8}{'sigma_max [MPa]':>18}{'Sm [MPa]':>12}{'Check':>10}")
    for name, smax in (("Inner", P_in), ("Outer", P_out)):
        verdict = "OK" if smax < S else "NOT OK"
        print(f" {name:<8}{smax:>18.2f}{S:>12.2f}{verdict:>10}")
    print("=" * 52)

    # 4) THERMAL STRESS
    Q_in, Q_out, sigma_th_in_diz, sigma_th_out_diz = thermal_stress(T_ci + 273.15, T_co + 273.15, D_in/2, D_out_clad/2)  # MPa 
    sigma_max_in_tot = P_in + Q_in
    sigma_max_out_tot = P_out + Q_out
    if max(sigma_max_in_tot) < 3*S and max(sigma_max_out_tot) < 3*S:
        print("The cladding is safe against combined mechanical and thermal stresses.")
    else:
        print("The cladding is NOT safe against combined mechanical and thermal stresses.")

    # calcolo delle tensioni totali (meccaniche + termiche) per ogni componente (h, r, l) e per interno/esterno
    sigma_h_in_tot = sigma_in_diz["h"] + sigma_th_in_diz["h"]
    sigma_h_out_tot = sigma_out_diz["h"] + sigma_th_out_diz["h"]
    sigma_l_in_tot = sigma_in_diz["l"] + sigma_th_in_diz["l"]
    sigma_l_out_tot = sigma_out_diz["l"] + sigma_th_out_diz["l"]
    sigma_r_in_tot = sigma_in_diz["r"] + sigma_th_in_diz["r"]
    sigma_r_out_tot = sigma_out_diz["r"] + sigma_th_out_diz["r"]

    # 5) GAS PLENUM SIZING  -- p_i_max = pressione massima ammissibile (Mariotte, sigma_h = sigma_y)
    V_plenum, H_plenum, n_tot = gas_plenum_analysis(T_gas, D_fuel, D_in, H_active, p_i)


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
