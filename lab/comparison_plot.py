import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# 1. Caricamento dei dati
# Nota: uso sep=';' perché Excel spesso esporta i CSV in italiano usando il punto e virgola
df = pd.read_csv('risultati_analisi.csv', sep=';')

# Imposto una dimensione del font leggermente più grande per la relazione
plt.rcParams.update({'font.size': 12})

# ==========================================
# GRAFICO 1: VOID FRACTION
# ==========================================
plt.figure(figsize=(8, 8))

# Disegno i punti sperimentali vs calcolati per ogni modello
plt.scatter(df['alpha_exp'], df['alpha_hom'], marker='s', color='black', label='Homogeneous', s=60, zorder=3)
plt.scatter(df['alpha_exp'], df['alpha_zivi'], marker='o', color='darkorange', label='Zivi', s=60, zorder=3)
plt.scatter(df['alpha_exp'], df['alpha_chisholm'], marker='^', color='silver', label='Chisholm', s=60, zorder=3)
plt.scatter(df['alpha_exp'], df['alpha_drift_flux'], marker='D', edgecolor='royalblue', facecolor='none', linewidth=1.5, label='Drift Flux', s=60, zorder=3)

# Creo i valori per le linee di riferimento (diagonale e +/- 20%)
x_vals = np.linspace(0, 1, 100)
plt.plot(x_vals, x_vals, color='crimson', label='diagonal', zorder=2)
plt.plot(x_vals, 1.2 * x_vals, color='crimson', linestyle='-.', label='+20%', zorder=2)
plt.plot(x_vals, 0.8 * x_vals, color='indianred', linestyle='--', label='-20%', zorder=2)

# Formattazione degli assi e della griglia
plt.xlim(0, 1.0)
plt.ylim(0, 1.0)
plt.xticks(np.arange(0, 1.1, 0.2))
plt.yticks(np.arange(0, 1.1, 0.1))
plt.grid(True, which='major', linestyle='-', linewidth=1.2, color='lightgray', zorder=0)

# Etichette
plt.xlabel('Experimental Void Fraction', fontsize=12)
plt.ylabel('Calculated Void Fraction', fontsize=12)
plt.title('Void Fraction Plot', fontsize=14)

# Aggiungo la legenda fuori dal grafico
plt.legend(bbox_to_anchor=(1.05, 0.5), loc='center left', frameon=False)

# Aggiusto i margini e salvo/mostro l'immagine
plt.tight_layout()
plt.savefig('void_fraction_plot.png', dpi=300, bbox_inches='tight')
plt.show()


# ==========================================
# GRAFICO 2: PRESSURE DROP (Tutte le 5 combinazioni)
# ==========================================
plt.figure(figsize=(8, 8))

# Trovo il valore massimo includendo TUTTE le colonne calcolate
max_val = max(
    df['dp_exp'].max(), 
    df[['dp_tot_hom_hom', 'dp_tot_friedel_zivi', 'dp_tot_friedel_chisholm', 'dp_tot_friedel_cise', 'dp_tot_friedel_drift_flux']].max().max()
)
# Arrotondo al multiplo di 5.000 superiore per "tagliare" gli assi su misura
max_val = np.ceil(max_val / 5000) * 5000

# Disegno i punti per TUTTE le 5 combinazioni
plt.scatter(df['dp_exp'], df['dp_tot_hom_hom'], marker='o', color='cornflowerblue', label='Homogeneous (Frict+Elev)', s=70, zorder=3)
plt.scatter(df['dp_exp'], df['dp_tot_friedel_zivi'], marker='s', color='darkorange', label='Friedel + Zivi', s=60, zorder=3)
plt.scatter(df['dp_exp'], df['dp_tot_friedel_chisholm'], marker='^', color='forestgreen', label='Friedel + Chisholm', s=70, zorder=3)
plt.scatter(df['dp_exp'], df['dp_tot_friedel_cise'], marker='D', color='crimson', label='Friedel + CISE', s=50, zorder=3)
plt.scatter(df['dp_exp'], df['dp_tot_friedel_drift_flux'], marker='*', color='purple', label='Friedel + Drift Flux', s=100, zorder=3)

# Creo i valori per le linee di riferimento (diagonale e +/- 20%)
x_vals2 = np.linspace(0, max_val, 100)
plt.plot(x_vals2, x_vals2, color='crimson', label='diagonal', zorder=2)
plt.plot(x_vals2, 1.2 * x_vals2, color='crimson', linestyle='-.', label='+20%', zorder=2)
plt.plot(x_vals2, 0.8 * x_vals2, color='indianred', linestyle='--', label='-20%', zorder=2)

# Formattazione assi e griglia
plt.xlim(0, max_val)
plt.ylim(0, max_val)
plt.grid(True, which='major', linestyle='-', linewidth=0.5, color='gray', zorder=0)

# Etichette
plt.xlabel('Experimental pressure drop (Pa)', fontsize=12)
plt.ylabel('Calculated pressure drop (Pa)', fontsize=12)
plt.title('Total Pressure Drop Plot', fontsize=14)

# Aggiungo la legenda fuori dal grafico
plt.legend(bbox_to_anchor=(1.05, 0.5), loc='center left', frameon=False)

# Aggiusto i margini e salvo/mostro l'immagine
plt.tight_layout()
plt.savefig('pressure_drop_all_models_plot.png', dpi=300, bbox_inches='tight')
plt.show()