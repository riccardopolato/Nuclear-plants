"""
Loader dei parametri di input da CSV per assignment2.

Il CSV ha colonne: name, value, unit, value_SI, unit_SI, description.
La sorgente di verità è (value, unit) in unità "originarie": il loader
converte sempre in SI usando le funzioni di unit_conversions.py.
Le colonne value_SI/unit_SI del CSV sono solo documentazione (utili per
la relazione) e vengono verificate per coerenza ma non usate dal codice.
"""

import csv
from pathlib import Path

from unit_conversions import (
    psia_to_pa,
    lbm_per_hr_to_kg_per_s,
    fahrenheit_to_celsius,
    inches_to_meters,
    square_feet_to_square_meters,
)


_CONVERTERS = {
    # potenze
    'MW':     (lambda x: x * 1e6,            'W'),
    'kW':     (lambda x: x * 1e3,            'W'),
    'W':      (lambda x: x,                  'W'),
    # pressioni
    'psia':   (psia_to_pa,                   'Pa'),
    'MPa':    (lambda x: x * 1e6,            'Pa'),
    'kPa':    (lambda x: x * 1e3,            'Pa'),
    'Pa':     (lambda x: x,                  'Pa'),
    # portate
    'lbm/hr': (lbm_per_hr_to_kg_per_s,       'kg/s'),
    'kg/s':   (lambda x: x,                  'kg/s'),
    # temperature
    'F':      (fahrenheit_to_celsius,        'C'),
    'C':      (lambda x: x,                  'C'),
    'K':      (lambda x: x - 273.15,         'C'),
    # lunghezze
    'in':     (inches_to_meters,             'm'),
    'ft':     (lambda x: x * 0.3048,         'm'),
    'm':      (lambda x: x,                  'm'),
    # aree
    'ft^2':   (square_feet_to_square_meters, 'm^2'),
    'm^2':    (lambda x: x,                  'm^2'),
    # adimensionali / counts
    '-':      (lambda x: x,                  '-'),
}


def load_inputs(csv_path='inputs.csv', check_si=True, rtol=1e-3):
    """
    Carica il CSV di input e converte ogni valore in unità SI interne.

    Parametri:
        csv_path: path al CSV (relativo a questo modulo se non assoluto).
        check_si: se True, verifica che value_SI nel CSV corrisponda
                  al SI calcolato (entro tolleranza rtol).
        rtol:     tolleranza relativa per il check di consistenza SI.

    Ritorna:
        values: dict {name: valore_SI}
        meta:   dict {name: {value_orig, unit_orig, value_SI, unit_SI, description}}
    """
    csv_path = Path(csv_path)
    if not csv_path.is_absolute():
        csv_path = Path(__file__).parent / csv_path

    values = {}
    meta = {}
    with open(csv_path, encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row['name'].strip()
            if not name or name.startswith('#'):
                continue

            value_orig = float(row['value'])
            unit_orig = row['unit'].strip()
            desc = row.get('description', '').strip()

            if unit_orig not in _CONVERTERS:
                raise ValueError(f"Unità '{unit_orig}' non riconosciuta per parametro '{name}'")

            converter, unit_SI_expected = _CONVERTERS[unit_orig]
            value_SI = converter(value_orig)

            if check_si and 'value_SI' in row and row['value_SI'].strip():
                value_SI_csv = float(row['value_SI'])
                if value_SI_csv != 0 and abs(value_SI - value_SI_csv) / abs(value_SI_csv) > rtol:
                    print(f"[WARN] '{name}': SI in CSV ({value_SI_csv:g}) ≠ SI calcolato ({value_SI:g})")

            values[name] = value_SI
            meta[name] = {
                'value_orig': value_orig,
                'unit_orig': unit_orig,
                'value_SI': value_SI,
                'unit_SI': unit_SI_expected,
                'description': desc,
            }

    return values, meta


def print_inputs_table(meta):
    """Stampa una tabella leggibile dei parametri."""
    header = f"{'Name':<18} {'Original':>14} {'Unit':<8} {'SI':>14} {'Unit_SI':<6} Description"
    print(header)
    print('-' * len(header) * 2)
    for name, m in meta.items():
        print(f"{name:<18} {m['value_orig']:>14.6g} {m['unit_orig']:<8} "
              f"{m['value_SI']:>14.6g} {m['unit_SI']:<6} {m['description']}")


if __name__ == '__main__':
    values, meta = load_inputs()
    print_inputs_table(meta)
