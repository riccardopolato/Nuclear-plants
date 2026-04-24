import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re
import runpy
from pathlib import Path


FIGURES_DIR = Path(__file__).resolve().parent / 'Figures'
FIGURES_DIR.mkdir(exist_ok=True)

# 1. Caricamento dei dati
# Nota: uso sep=';' perché Excel spesso esporta i CSV in italiano usando il punto e virgola
df = pd.read_csv('lab/risultati_analisi.csv', sep=';')

# Imposto una dimensione del font leggermente più grande per la relazione
plt.rcParams.update({'font.size': 12})

# %% Plot 1: Void Fraction
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
plt.savefig(FIGURES_DIR / 'void_fraction_plot.png', dpi=300, bbox_inches='tight')
#plt.show()


# %% Plot 2: Pressure Drop
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
plt.savefig(FIGURES_DIR / 'pressure_drop_all_models_plot.png', dpi=300, bbox_inches='tight')



# %% Plot con Diagram plotter
# Hewitt-Roberts 
# Eseguo lo script originale
runpy.run_path('lab/Diagram plotter/Hewitt_Roberts.py')

# Estraggo i dati
df_data = pd.read_csv('lab/risultati_analisi.csv', sep=';')

# Uso le densità calcolate per ogni shot (rho_l e rho_g) dal CSV risultati
df_data['x_hewitt'] = df_data['rho_l'] * (df_data['j_l'] ** 2)
df_data['y_hewitt'] = df_data['rho_g'] * (df_data['j_g'] ** 2)


def pattern_to_sigla(pattern):
    mapping = {
        'slug': 'S',
        'churn': 'C',
        'annular': 'A',
        'anular': 'A',
        'bubble': 'B',
        'bubbly': 'B'
    }
    parts = [p.strip() for p in str(pattern).split('/')]
    sigle = []
    for part in parts:
        key = part.lower()
        if key in mapping:
            sigle.append(mapping[key])
        elif part:
            sigle.append(part[0].upper())
    return '/'.join(sigle)

# Aggiungo i miei punti in rosso sul grafico esistente
for idx, row in df_data.iterrows():
    x_point = row['x_hewitt']
    y_point = row['y_hewitt']
    shot = row['shot_id']
    pattern = row['flow_pattern']
    
    # Plot punto rosso
    plt.scatter(x_point, y_point, color='red', s=80, marker='o', zorder=3, edgecolors='darkred', linewidth=1.5)
    
    # Etichetta compatta: shot:sigla_pattern (es. 45:S/C)
    label_text = f"{int(shot)}:{pattern_to_sigla(pattern)}"
    plt.annotate(label_text, 
                xy=(x_point, y_point),
                xytext=(5, 5),
                textcoords='offset points',
                fontsize=7,
                color='darkred',
                weight='bold',
                zorder=4)

plt.savefig(FIGURES_DIR / 'hewitt_roberts_with_data.png', dpi=300, bbox_inches='tight')


# %% Plot con Diagram plotter
# Taitel-Dukler: Slug/Churn e Bubbly/Slug-Churn

G = 9.81
D_TUBE = 26e-3  # 26 mm in metri
SIGMA_WATER = 0.072  # N/m, valore tipico in condizioni ambiente


def normalize_pattern_tokens(pattern):
    text = str(pattern).strip().lower()
    text = text.replace('anular', 'annular')
    tokens = [t for t in re.split(r'[^a-z]+', text) if t]
    return set(tokens)


def classify_td_target(pattern):
    # Mantengo la funzione per retrocompatibilita, ma i filtri TD sono ora indipendenti.
    return None


def canonical_sigla(pattern):
    """Converte il flow pattern in una sigla canonica (es. S/C, B/S, C/A)."""
    tokens = normalize_pattern_tokens(pattern)
    has_slug = 'slug' in tokens or 'slugs' in tokens
    has_churn = 'churn' in tokens
    has_bubble = 'bubble' in tokens or 'bubbly' in tokens or 'bubbles' in tokens
    has_annular = 'annular' in tokens

    if has_bubble and has_slug:
        return 'B/S'
    if has_slug and has_churn:
        return 'S/C'
    if has_churn and has_annular:
        return 'C/A'
    if has_bubble:
        return 'B'
    if has_slug:
        return 'S'
    if has_churn:
        return 'C'
    if has_annular:
        return 'A'
    return pattern_to_sigla(pattern)


def include_td_slug_churn(pattern):
    # Accetta: S, S/C, C
    return canonical_sigla(pattern) in {'S', 'S/C', 'C'}


def include_td_bubbly_slugchurn(pattern):
    # Accetta: B, S, C, B/S, S/C
    return canonical_sigla(pattern) in {'B', 'S', 'C', 'B/S', 'S/C'}


def include_td_annular_slugchurn(pattern):
    # Accetta: S/C, C, C/A, A
    return canonical_sigla(pattern) in {'S/C', 'C', 'C/A', 'A'}


td_data = pd.read_csv('lab/risultati_analisi.csv', sep=';')
td_data = td_data.dropna(subset=['j_l', 'j_g', 'rho_l', 'rho_g', 'alpha_exp', 'flow_pattern'])

# --- TD_Sl+Ch: x = J/sqrt(gD), y = alpha_exp ---
runpy.run_path('lab/Diagram plotter/TD_Sl+Ch.py')

for _, row in td_data.iterrows():
    if not include_td_slug_churn(row['flow_pattern']):
        continue

    j_total = row['j_l'] + row['j_g']
    x_td_slch = j_total / np.sqrt(G * D_TUBE)
    y_td_slch = row['alpha_exp']
    label_text = f"{int(row['shot_id'])}:{pattern_to_sigla(row['flow_pattern'])}"

    plt.scatter(
        x_td_slch,
        y_td_slch,
        color='red',
        s=80,
        marker='o',
        zorder=3,
        edgecolors='darkred',
        linewidth=1.5
    )

    plt.annotate(
        label_text,
        xy=(x_td_slch, y_td_slch),
        xytext=(5, 5),
        textcoords='offset points',
        fontsize=7,
        color='darkred',
        weight='bold',
        zorder=4,
    )

plt.savefig(FIGURES_DIR / 'td_slug_churn_with_data.png', dpi=300, bbox_inches='tight')

# --- TD_Bb+SlCh: x = jg*sqrt(rho_g)/(g*(rho_l-rho_g)*sigma)^(1/4), y = jl/jg ---
runpy.run_path('lab/Diagram plotter/TD_Bb+SlCh.py')

for _, row in td_data.iterrows():
    if not include_td_bubbly_slugchurn(row['flow_pattern']):
        continue

    if row['j_g'] <= 0:
        continue

    rho_diff = row['rho_l'] - row['rho_g']
    if rho_diff <= 0:
        continue

    x_td_bb = (row['j_g'] * np.sqrt(row['rho_g'])) / ((G * rho_diff * SIGMA_WATER) ** 0.25)
    y_td_bb = row['j_l'] / row['j_g']
    label_text = f"{int(row['shot_id'])}:{pattern_to_sigla(row['flow_pattern'])}"

    plt.scatter(
        x_td_bb,
        y_td_bb,
        color='red',
        s=80,
        marker='o',
        zorder=3,
        edgecolors='darkred',
        linewidth=1.5
    )

    plt.annotate(
        label_text,
        xy=(x_td_bb, y_td_bb),
        xytext=(5, 5),
        textcoords='offset points',
        fontsize=7,
        color='darkred',
        weight='bold',
        zorder=4,
    )

plt.savefig(FIGURES_DIR / 'td_bubbly_slugchurn_with_data.png', dpi=300, bbox_inches='tight')

# --- TD_An+SlCh: x = X_martinelli, y = jg*sqrt(rho_g)/(g*(rho_l-rho_g)*sigma)^(1/4) ---
runpy.run_path('lab/Diagram plotter/TD_An+SlCh.py')

for _, row in td_data.iterrows():
    if not include_td_annular_slugchurn(row['flow_pattern']):
        continue

    if pd.isna(row.get('X_martinelli')) or row['X_martinelli'] <= 0:
        continue

    rho_diff = row['rho_l'] - row['rho_g']
    if rho_diff <= 0:
        continue

    x_td_an = row['X_martinelli']
    y_td_an = (row['j_g'] * np.sqrt(row['rho_g'])) / ((G * rho_diff * SIGMA_WATER) ** 0.25)
    label_text = f"{int(row['shot_id'])}:{pattern_to_sigla(row['flow_pattern'])}"

    plt.scatter(
        x_td_an,
        y_td_an,
        color='red',
        s=80,
        marker='o',
        zorder=3,
        edgecolors='darkred',
        linewidth=1.5
    )

    plt.annotate(
        label_text,
        xy=(x_td_an, y_td_an),
        xytext=(5, 5),
        textcoords='offset points',
        fontsize=7,
        color='darkred',
        weight='bold',
        zorder=4,
    )

plt.savefig(FIGURES_DIR / 'td_annular_slugchurn_with_data.png', dpi=300, bbox_inches='tight')
plt.show()