import pandas as pd
import math
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
        'p_testsection', 'type_dp', 'dp_transd', 'T_water', 'm_l', 'flow_pattern'
    ]
    df.columns = new_column_names

    # Pulisci i dati
    df = df.dropna(how='all')
    df = df.dropna(subset=['test'])

    # Aggiungi M_ref alla colonna m_l
    df['m_l'] = df['m_l']/1000 + m_ref

    # Crea il dizionario raggruppando per 'test'
    all_data = {}
    for shot_id, group in df.groupby('test'):
        all_data[shot_id] = group.to_dict('records')
    
    return all_data


def calculate_water_flow_rate(diaphragm_type, delta_p_mbar, diaphragm_data, temperature_c):
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
    
    temperature_k = temperature_c + 273.15

    rho = PropsSI('D', 'T', temperature_k, 'P', 101325, 'Water')

    # Conversione della caduta di pressione da mbar a Pascal (1 mbar = 100 Pa)
    delta_p_pa = delta_p_mbar * 100
    
    # Calcolo dell'area della sezione di passaggio del diaframma
    area = (math.pi * d**2) / 4
    
    # Calcolo della portata massica (W) in kg/s
    W = alpha_mq * area * math.sqrt(2 * rho * delta_p_pa) + B
    
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


def analisi_exp(all_data, diaphragm_data, M_l0, rho, h, g):
    """
    Itera attraverso il dizionario di dati, calcola la portata per ogni record
    e la aggiunge al record stesso.
    """
    for shot_id, records in all_data.items():
        for record in records:
            # Calcola la portata (W) e aggiungila al record
            W = calculate_water_flow_rate(record['diaphragm'], record['dp_dia_water'], diaphragm_data, record['T_water'])
            
            W_air = calculate_air_flow_rate(record['Q_s'], record['p_e_rel'])
            void_fraction = 1 - record['m_l']/M_l0 

            # Aggiungi i risultati al record corrente
            record['W_water'] = W
            record['W_air'] = W_air
            record['void_fraction_exp'] = void_fraction

            # calcolo caduta di pressione
            if record['type_dp'] == 'pD-pC':
                record['dp_exp'] = rho*g*h - record['dp_transd']
            elif record['type_dp'] == 'pC-pD':
                record['dp_exp'] = rho*g*h + record['dp_transd']




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

    # 1. Leggi e prepara i dati
    dati_esperimento = leggi_e_prepara_dati(file_da_leggere, M_res)

    # 2. calcolo della void fraction sperimentale e aggiunta dei risultati al dizionario
    analisi_exp(dati_esperimento, diaphragm_data, M_l0, rho_water, h, g)

    # 3. Calcolo caduta di pressione

    # Se vuoi esportare il risultato pulito in un nuovo file CSV (richiede modifiche alle funzioni):
    # df.to_csv('tab_dat_flowpat_pulito.csv', sep=';', index=False, na_rep='NaN')
