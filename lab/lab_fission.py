from numpy import record
import pandas as pd
import math
import csv
import CoolProp.CoolProp as CP
from CoolProp.CoolProp import PropsSI

def leggi_e_prepara_dati(file_path, m_ref):
    """
    Legge un file CSV, lo pulisce, rinomina le colonne, aggiunge un valore di riferimento
    e restituisce i dati raggruppati in un dizionario.
    """
    # Leggi il file CSV
    df = pd.read_csv(file_path, sep=';', na_values=['NULL', 'null'])

    # Assegna i nuovi nomi alle colonne
    new_column_names = [
        'test', 'diaphragm', 'Q_s', 'p_e_rel', 'dp_dia_water',
        'p_rel_test', 'type_dp', 'dp_transd', 'T_water', 'm_l', 'flow_pattern'
    ]
    df.columns = new_column_names

    # Forzare la conversione a numerico delle colonne rilevanti, gestendo errori
    cols_to_numeric = ['Q_s', 'p_e_rel', 'dp_dia_water', 'p_rel_test', 'dp_transd', 'T_water', 'm_l']
    for col in cols_to_numeric:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    # Pulisci i dati
    df = df.dropna(how='all')
    df = df.dropna(subset=['test'])

    # Aggiungi M_ref alla colonna m_l
    df['m_l'] = df['m_l']/1000 + m_ref
    df['dp_transd'] = df['dp_transd']*100 # conversione da mbar a Pa
    df['p_rel_test'] = df['p_rel_test']*100000 # conversione da bar a Pa
    df['dp_dia_water'] = df['dp_dia_water']*100 # conversione da mbar a Pa
    # Crea il dizionario raggruppando per 'test'
    all_data = {}
    for shot_id, group in df.groupby('test'):
        all_data[shot_id] = group.to_dict('records')
    
    return all_data


def calculate_water_flow_rate(diaphragm_type, delta_p, diaphragm_data, temperature_c):
    """
    Calcola la portata massica dell'acqua in kg/s.
    
    Parametri:
    - diaphragm_type: stringa, 'S' per small o 'M' per medium.
    - delta_p_mbar: float, caduta di pressione in mbar misurata dal trasduttore differenziale.
    - diaphragm_data: dizionario con i dati dei diaframmi.
    - temperature_c: float, temperatura dell'acqua in gradi Celsius.
    """
    
    # Controllo che il tipo di diaframma sia tra quelli salvati ('S' o 'M')
    if diaphragm_type not in diaphragm_data:
        raise ValueError("Tipo di diaframma non valido. Usa 'S' o 'M'.")
        
    # Estrazione dei parametri corretti dal dizionario
    params = diaphragm_data[diaphragm_type]
    d = params['d']
    alpha_mq = params['alpha_mq']
    B = params['B']

    rho = PropsSI('D', 'T', temperature_c + 273.15, 'P', 101325, 'Water')
    
    # Calcolo dell'area della sezione di passaggio del diaframma
    area = (math.pi * d**2) / 4
    
    # Calcolo della portata massica (W) in kg/s
    W = alpha_mq * area * math.sqrt(2 * rho * delta_p) + B
    
    return W


def calculate_air_flow_rate(Qs, p_rel):
    """
    Calcola la portata massica dell'aria in kg/s.
    
    Parametri:
    - Qs: float, portata volumetrica letta dal rotametro in Nm^3/h.
    - p_rel: float, pressione relativa misurata a monte del rotametro in barg.
    """
    
    # Calcolo della portata massica in kg/h usando la formula del manuale
    W_air_kg_h = 1.204 * math.sqrt((p_rel + 1.013) / 1.013) * Qs
    
    # Conversione in kg/s per avere la stessa unità di misura dell'acqua
    W_air_kg_s = W_air_kg_h / 3600.0
    
    return W_air_kg_s


def analisi_exp(all_data, diaphragm_data, M_l0, rho, h, g, A_cross):
    """
    Itera attraverso il dizionario di dati, calcola la portata per ogni record
    e la aggiunge al record stesso.
    """
    for shot_id, records in all_data.items():
        for record in records:
            # calcolo densità
            record['rho_l'] = CP.PropsSI('D', 'T', record['T_water'] + 273.15, 'P', 101325, 'Water')
            record['rho_g'] = CP.PropsSI('D', 'T', 20 + 273.15, 'P', 101325 + record['p_rel_test'], 'Air')
            # Calcola la portata (W) e aggiungila al record
            W = calculate_water_flow_rate(record['diaphragm'], record['dp_dia_water'], diaphragm_data, record['T_water'])
            
            W_air = calculate_air_flow_rate(record['Q_s'], record['p_e_rel'])
            void_fraction = 1 - record['m_l']/M_l0 

            # Aggiungi i risultati al record corrente
            record['W_water'] = W #[kg/s]
            record['W_air'] = W_air #[kg/s]
            record['void_fraction_exp'] = void_fraction

            # calcolo caduta di pressione
            if record['type_dp'] == 'pD-pC':
                record['dp_exp'] = rho*g*h - record['dp_transd'] #[Pa]
            elif record['type_dp'] == 'pC-pD':
                record['dp_exp'] = rho*g*h + record['dp_transd'] #[Pa]
            elif record['type_dp'] == 'Null':
                record['dp_exp'] = None
            
            # calcolo titolo
            record['x_exp'] = record['W_air']/(record['W_air'] + record['W_water'])

            # calcolo mass flux
            record['G_exp'] = (record['W_air'] + record['W_water'])/A_cross #[kg/m^2/s]

            #calcolo velocità superficiali aria: T_amb, p_atm; water: T_water, p_atm
            record['j_g'] = record['W_air']/A_cross/record['rho_g'] #[m/s]
            record['j_l'] = record['W_water']/A_cross/record['rho_l'] #[m/s]
            record['j_tot'] = record['j_g'] + record['j_l'] #[m/s]

# calcolo void fraction con i modelli
def void_fraction_homogeneous(x, rho_l, rho_g):
    """Calcola il void fraction con il modello omogeneo."""
    # S = 1 per il modello omogeneo
    alpha = 1 / (1 + ((1 - x) / x) * (rho_g / rho_l))
    return alpha

def void_fraction_zivi(x, rho_l, rho_g):
    """Calcola il void fraction con la correlazione di Zivi."""
    S = (rho_l / rho_g)**(1/3)
    alpha = 1 / (1 + S * ((1 - x) / x) * (rho_g / rho_l))
    return alpha

def void_fraction_chisholm(x, rho_l, rho_g):
    """Calcola il void fraction con la correlazione di Chisholm."""
    # Attenzione: la radice quadrata potrebbe avere un argomento negativo se rho_l > rho_g
    # e x è vicino a 1. Aggiungiamo un controllo.
    arg = 1 - x * (1 - rho_l / rho_g)
    if arg < 0:
        return float('nan') # Restituisce Not a Number se l'argomento è negativo
    S = math.sqrt(arg)
    alpha = 1 / (1 + S * ((1 - x) / x) * (rho_g / rho_l))
    return alpha

def void_fraction_cise(x, rho_l, rho_g, G, D, mu_l, sigma):
    """Calcola il void fraction con la correlazione CISE."""
    # Calcolo di beta
    beta = (rho_l * x) / (rho_l * x + rho_g * (1 - x))
    
    # Calcolo di y
    y = beta / (1 - beta)
    
    # Calcolo dei numeri adimensionali Re e We
    Re = (G * D) / mu_l
    We = (G**2 * D) / (sigma * rho_l)
    
    # Calcolo di E1 e E2
    E1 = 1.578 * (Re**(-0.19)) * ((rho_l / rho_g)**0.22)
    E2 = 0.0273 * We * (Re**(-0.51)) * ((rho_l / rho_g)**(-0.08))
    
    # Calcolo di S
    term = (y / (1 + y * E2)) - y * E2
    S = 1 + E1 * (term**0.5)

    # Calcolo di alpha
    alpha = 1 / (1 + S * ((1 - x) / x) * (rho_g / rho_l))
    return alpha

def void_fraction_drift_flux(W_g, W_l, rho_g, rho_l, sigma, D, A, flow_pattern, g=9.81):
    """Calcola il void fraction con il Drift Flux Model."""
    
    # Calcolo delle velocità superficiali
    j_g = W_g / (rho_g * A)
    j_l = W_l / (rho_l * A)
    j = j_g + j_l

    # # Selezione dei parametri in base al flow pattern
    # if 'Bubble' in flow_pattern:
    #     C0 = 1.13
    #     # Calcolo velocità di deriva per bubbly flow
    #     term = (sigma * g * (rho_l - rho_g)) / (rho_l**2)
    #     u_gj = 1.41 * (term**0.25)
    # elif 'Slug' in flow_pattern:
    #     C0 = 1.2
    #     # Calcolo velocità di deriva per plug flow
    #     term = ((rho_l - rho_g) * g * D) / rho_l
    #     u_gj = 0.35 * math.sqrt(term)
    # else:
    #     # Se il flow pattern non è riconosciuto, non possiamo calcolare alpha
    #     return float('nan')

    C0 = 1.13
    # Calcolo velocità di deriva per bubbly flow
    term = (sigma * g * (rho_l - rho_g)) / (rho_l**2)
    u_gj = 1.41 * (term**0.25)
    # Calcolo del void fraction alpha
    denominator = C0 * j + u_gj
    if denominator == 0:
        return 0.0
        
    alpha = j_g / denominator
    return alpha


def dp_el(alpha, rho_l, rho_g, g, h):
    """Calcola la caduta di pressione idrostatica."""
    return (alpha * rho_g + (1 - alpha) * rho_l) * g * h

def get_friction_factor(Re):
    """
    Calcola il fattore di attrito di Darcy per tubo liscio.
    """
    if Re < 2000:
        return 64.0 / Re  # Moto laminare Darcy
    else:
        return 0.316 * (Re ** -0.25)  # Moto turbolento (Blasius)
    
def calculate_homogeneous_friction_drop(x, G, D, h, rho_h, mu_l, mu_g):
    mu_h = x * mu_g + (1 - x) * mu_l
    Re_h = (G * D) / mu_h
    f_hom = get_friction_factor(Re_h)
    return f_hom * (h / D) * G**2 / (2 * rho_h) 


def calculate_friedel_friction_drop(x, G, D, h, rho_l, rho_g, rho_h, mu_l, mu_g, sigma):
    """
    Calcola la caduta di pressione per attrito usando il modello di Friedel.
    
    Parametri:
    - x: titolo termodinamico (mass quality)
    - G: flusso di massa (kg/m^2 s)
    - D: diametro interno del tubo (m)
    - h: lunghezza della sezione di test (m)
    - rho_l, rho_g, rho_h: densità liquido, gas e omogenea (kg/m^3)
    - mu_l, mu_g: viscosità dinamica liquido e gas (Pa*s)
    - sigma: tensione superficiale (N/m)
    """
    
    # 1. Numeri adimensionali (Froude e Weber omogenei)
    Fr_h = (G**2) / (g * D * rho_h**2)
    We_h = (G**2 * D) / (sigma * rho_h)
    
    # 2. Calcolo dei Reynolds assumendo che tutto il flusso sia liquido (l0) o gas (g0)
    Re_l0 = (G * D) / mu_l
    Re_g0 = (G * D) / mu_g
    
    # 3. Calcolo dei fattori di attrito per fase singola
    f_l0 = get_friction_factor(Re_l0)
    f_g0 = get_friction_factor(Re_g0)
    
    # 4. Parametri di Friedel (E, F, H)
    E = (1 - x)**2 + (x**2 * (rho_l / rho_g) * (f_g0 / f_l0))
    F = (x**0.78) * ((1 - x)**0.224)
    H = ((rho_l / rho_g)**0.91) * ((mu_g / mu_l)**0.19) * ((1 - (mu_g / mu_l))**0.7)
    
    # 5. Moltiplicatore bifase di Friedel
    phi_l0_sq = E + ((3.24 * F * H) / ((Fr_h**0.045) * (We_h**0.035)))
    
    # 6. Caduta di pressione assumendo solo liquido (Delta p_l0)
    dp_l0 = f_l0 * (h / D) * ((G**2) / (2 * rho_l))
    
    # 7. Caduta di pressione totale per attrito
    dp_frict_friedel = dp_l0 * phi_l0_sq
    
    return dp_frict_friedel


def calculate_single_phase_dp_dz(W_phase, rho_phase, mu_phase, D, h, A, g=9.81):
    """
    Calcola (dp/dz) monofase [Pa/m] come:
    (Delta p idrostatica + Delta p attrito) / h.
    """
    if h <= 0 or A <= 0 or D <= 0 or rho_phase <= 0 or mu_phase <= 0:
        return float('nan')

    G_phase = W_phase / A
    Re_phase = (G_phase * D) / mu_phase
    if Re_phase <= 0:
        return float('nan')

    f_phase = get_friction_factor(Re_phase)
    dp_fric = f_phase * (h / D) * (G_phase**2) / (2 * rho_phase)
    dp_hydro = rho_phase * g * h
    return (dp_hydro + dp_fric) / h


def calcola_valori_derivati(all_data, diaphragm_data, M_l0, D, A):
    """
    Itera attraverso i dati, calcola le portate, le proprietà termofisiche
    e i void fraction secondo vari modelli, arricchendo i record.
    """
    for shot_id, records in all_data.items():
        for record in records:

            mu_l = CP.PropsSI('V', 'T', record['T_water'] + 273.15, 'P', 101325, 'Water')
            mu_g = CP.PropsSI('V', 'T', 20 + 273.15, 'P', 101325 + record['p_rel_test'], 'Air')
            sigma = PropsSI('I', 'T', record['T_water'] + 273.15, 'Q', 0, 'Water') # Use Q=0 for saturated liquid property
            # --- 2. Calcolo Void Fraction con i modelli ---
            
            # Modello Omogeneo
            alpha_hom = void_fraction_homogeneous(record['x_exp'], record['rho_l'], record['rho_g'])
            
            # Correlazione di Zivi
            alpha_zivi = void_fraction_zivi(record['x_exp'], record['rho_l'], record['rho_g'])
            
            # Correlazione di Chisholm
            alpha_chisholm = void_fraction_chisholm(record['x_exp'], record['rho_l'], record['rho_g'])
            
            # Correlazione CISE
            alpha_cise = void_fraction_cise(record['x_exp'], record['rho_l'], record['rho_g'], record['G_exp'], D, mu_l, sigma)

            # Correlazione Drift Flux
            alpha_drift_flux = void_fraction_drift_flux(record['W_air'], record['W_water'], record['rho_g'], record['rho_l'], sigma, D, A, record['flow_pattern'])


            # Aggiungi i void fraction calcolati al record
            record['alpha_hom'] = alpha_hom
            record['alpha_zivi'] = alpha_zivi
            record['alpha_chisholm'] = alpha_chisholm
            record['alpha_cise'] = alpha_cise
            record['alpha_drift_flux'] = alpha_drift_flux

            # Aggiungi idp elevazione calcolati al record
            record['dp_hom'] = dp_el(alpha_hom, record['rho_l'], record['rho_g'], g, h)
            record['dp_zivi'] = dp_el(alpha_zivi, record['rho_l'], record['rho_g'], g, h)
            record['dp_chisholm'] = dp_el(alpha_chisholm, record['rho_l'], record['rho_g'], g, h)
            record['dp_cise'] = dp_el(alpha_cise, record['rho_l'], record['rho_g'], g, h)
            record['dp_drift_flux'] = dp_el(alpha_drift_flux, record['rho_l'], record['rho_g'], g, h)

            # delta p attrito omogeneo
            rho_h = ((record['x_exp'] / record['rho_g']) + ((1 - record['x_exp']) / record['rho_l']))**(-1)
            dp_fric_hom = calculate_homogeneous_friction_drop(record['x_exp'], record['G_exp'], D, h, rho_h, mu_l, mu_g)
            record['dp_fric_hom'] = dp_fric_hom

            # delta p attrito friedel
            dp_fric_friedel = calculate_friedel_friction_drop(record['x_exp'], record['G_exp'], D, h, record['rho_l'], record['rho_g'], rho_h, mu_l, mu_g, sigma)
            record['dp_fric_friedel'] = dp_fric_friedel

            # Parametro di Martinelli X = sqrt((dp/dz)_l / (dp/dz)_g)
            dp_dz_l = calculate_single_phase_dp_dz(record['W_water'], record['rho_l'], mu_l, D, h, A, g)
            dp_dz_g = calculate_single_phase_dp_dz(record['W_air'], record['rho_g'], mu_g, D, h, A, g)
            if dp_dz_g <= 0 or math.isnan(dp_dz_l) or math.isnan(dp_dz_g):
                record['X_martinelli'] = float('nan')
            else:
                record['X_martinelli'] = math.sqrt(dp_dz_l / dp_dz_g)

            record['dp_dz_l_single'] = dp_dz_l
            record['dp_dz_g_single'] = dp_dz_g


def crea_dizionario_finale(all_data):
    """
    Costruisce un dizionario finale con, per ogni record:
      - tutti gli alpha (sperimentale + correlazioni)
      - dp_exp (totale misurata)
      - dp totali predette:
          * hom (friction) + hom (elevation)
          * friedel (friction) + {zivi, chisholm, cise, drift_flux} (elevation)
    Arrotonda tutti i valori numerici a 4 cifre decimali.
    """
    risultati = {}

    for shot_id, records in all_data.items():
        risultati[shot_id] = []

        for record in records:
            # Funzione helper per arrotondare i valori se sono numerici
            def round_if_numeric(value, decimals=4):
                if isinstance(value, (int, float)):
                    try:
                        return round(value, decimals)
                    except (TypeError, ValueError):
                        # Restituisce il valore originale se non può essere arrotondato
                        return value
                return value

            # Calcola le somme prima per chiarezza
            dp_tot_hom_hom = record.get('dp_hom', 0) + record.get('dp_fric_hom', 0)
            dp_tot_friedel_zivi = record.get('dp_zivi', 0) + record.get('dp_fric_friedel', 0)
            dp_tot_friedel_chisholm = record.get('dp_chisholm', 0) + record.get('dp_fric_friedel', 0)
            dp_tot_friedel_cise = record.get('dp_cise', 0) + record.get('dp_fric_friedel', 0)
            dp_tot_friedel_drift_flux = record.get('dp_drift_flux', 0) + record.get('dp_fric_friedel', 0)

            nuovo = {
                # Identificativi utili (facoltativi, rimuovili se non ti servono)
                'test': record.get('test'),
                'flow_pattern': record.get('flow_pattern'),
                'x_exp': round_if_numeric(record.get('x_exp')),
                'G_exp': round_if_numeric(record.get('G_exp')),
                'j_g': round_if_numeric(record.get('j_g')),
                'j_l': round_if_numeric(record.get('j_l')),
                'rho_l': round_if_numeric(record.get('rho_l')),
                'rho_g': round_if_numeric(record.get('rho_g')),
                'W_water': round_if_numeric(record.get('W_water')),
                'W_air': round_if_numeric(record.get('W_air')),
                'dp_dz_l_single': round_if_numeric(record.get('dp_dz_l_single')),
                'dp_dz_g_single': round_if_numeric(record.get('dp_dz_g_single')),
                'X_martinelli': round_if_numeric(record.get('X_martinelli')),

                # --- Void fraction ---
                'alpha_exp':        round_if_numeric(record.get('void_fraction_exp')),
                'alpha_hom':        round_if_numeric(record.get('alpha_hom')),
                'alpha_zivi':       round_if_numeric(record.get('alpha_zivi')),
                'alpha_chisholm':   round_if_numeric(record.get('alpha_chisholm')),
                'alpha_cise':       round_if_numeric(record.get('alpha_cise')),
                'alpha_drift_flux': round_if_numeric(record.get('alpha_drift_flux')),

                # --- Caduta di pressione sperimentale (totale misurata) ---
                'dp_exp': round_if_numeric(record.get('dp_exp')),

                # --- Cadute di pressione predette (elevation + friction) ---
                'dp_tot_hom_hom': round_if_numeric(dp_tot_hom_hom),
                'dp_tot_friedel_zivi': round_if_numeric(dp_tot_friedel_zivi),
                'dp_tot_friedel_chisholm': round_if_numeric(dp_tot_friedel_chisholm),
                'dp_tot_friedel_cise': round_if_numeric(dp_tot_friedel_cise),
                'dp_tot_friedel_drift_flux': round_if_numeric(dp_tot_friedel_drift_flux),
            }
            risultati[shot_id].append(nuovo)

    return risultati



def esporta_csv(risultati_finali, file_output='risultati.csv'):
    """
    Esporta il dizionario finale in un unico CSV (long format).
    Ogni riga = un record. Il campo 'shot_id' identifica il test.
    """
    # Colonne nell'ordine in cui le vuoi nel CSV
    colonne = [
        'shot_id', 'test', 'flow_pattern',
        'x_exp', 'G_exp', 'j_g', 'j_l', 'rho_l', 'rho_g', 'W_water', 'W_air',
        'dp_dz_l_single', 'dp_dz_g_single', 'X_martinelli',
        # Void fraction
        'alpha_exp', 'alpha_hom', 'alpha_zivi',
        'alpha_chisholm', 'alpha_cise', 'alpha_drift_flux',
        # Cadute di pressione [Pa]
        'dp_exp',
        'dp_tot_hom_hom',
        'dp_tot_friedel_zivi',
        'dp_tot_friedel_chisholm',
        'dp_tot_friedel_cise',
        'dp_tot_friedel_drift_flux',
    ]

    with open(file_output, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=colonne, delimiter=';')
        writer.writeheader()
        for shot_id, records in risultati_finali.items():
            for rec in records:
                riga = {'shot_id': shot_id}
                riga.update({k: rec.get(k) for k in colonne if k != 'shot_id'})
                writer.writerow(riga)

    print(f"CSV salvato in: {file_output}")

def esporta_portate_csv(risultati_finali, file_output='portate.csv'):
    """
    Esporta le portate di aria e acqua in un file CSV separato.
    """
    colonne = [
        'shot_id', 'test', 'W_air', 'W_water'
    ]

    with open(file_output, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=colonne, delimiter=';')
        writer.writeheader()
        for shot_id, records in risultati_finali.items():
            for rec in records:
                riga = {'shot_id': shot_id}
                riga.update({k: rec.get(k) for k in colonne if k != 'shot_id'})
                writer.writerow(riga)

    print(f"CSV con portate salvato in: {file_output}")

# --- Blocco di esecuzione principale ---
if __name__ == "__main__":
    # Definizioni iniziali
    file_da_leggere = 'lab/tab_dat_flowpat.csv'
    # Dati della tabella per i diaframmi estratti dal manuale
    diaphragm_data = {
        'S': {
            'd': 0.0084,          # Diametro (m)
            'alpha_mq': 0.650729, # Costante a_MQ
            'B': 0.01415274       # Costante B
        },
        'M': {
            'd': 0.0153,          # Diametro (m)
            'alpha_mq': 0.728193, # Costante a_MQ
            'B': 0.00375718       # Costante B
        }
    }
    M_res = 0.116 #[kg]
    M_l0 = 1.014 #[kg]
    h = 1.5 #[m]
    g = 9.81 #[m/s^2]
    rho_water = 1000 #[kg/m^3]
    Diameter = 0.026 #[m]
    A_cross = math.pi * (Diameter/2)**2 #[m^2]

    # 1. Leggi e prepara i dati
    dati_esperimento = leggi_e_prepara_dati(file_da_leggere, M_res)

    # 2. calcolo della void fraction sperimentale e aggiunta dei risultati al dizionario
    analisi_exp(dati_esperimento, diaphragm_data, M_l0, rho_water, h, g, A_cross)

    # 3. Calcolo dei risultati con correlazioni
    calcola_valori_derivati(dati_esperimento, diaphragm_data, M_l0, Diameter, A_cross)

    # 4. Dizionario finale "pulito"
    risultati_finali = crea_dizionario_finale(dati_esperimento)

    # 5. Export CSV
    esporta_csv(risultati_finali, 'lab/risultati_analisi.csv')
    esporta_portate_csv(risultati_finali, 'lab/portate.csv')
