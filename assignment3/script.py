import numpy as np

from inputs_loader import load_inputs




# 1) BUCKLING ANALYSIS
def buckling(T, D_out, t):
    D_in = D_out - 2*t
    r_avg = (D_out + D_in) / 2
    E =  (9.9e3 - 5.669*(T-273))*9.81 #MPa
    nu = 0.3303 + 8.376e-5*(T-273)
    p_cr = E/(4*(1-nu**2))*()



def main():
    # importare dati da file csv
    inputs, _ = load_inputs()  # legge inputs.csv in questa stessa cartella

    D_out_clad = inputs["D_out_clad"]  # m, diametro esterno guaina
    thickness  = inputs["s_clad"]      # m, spessore guaina

    print(f"D_out_clad = {D_out_clad:.6e} m")
    print(f"thickness  = {thickness:.6e} m")


if __name__ == "__main__":
    main()
