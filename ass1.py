# se entrate sappiate che vi devo accettare 
import CoolProp.CoolProp as CP
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import numpy as np

# Definiamo un range di pressioni da 1 bar a 220 bar (vicino al punto critico)
pressioni = np.linspace(1e5, 220e5, 100)
temperature = [CP.PropsSI('T', 'P', p, 'Q', 0, 'Water') - 273.15 for p in pressioni]

def complement_color(color):
	"""Restituisce il colore complementare, lasciando intatti valori speciali (es. 'inherit')."""
	try:
		r, g, b = mcolors.to_rgb(color)
		return (1 - r, 1 - g, 1 - b)
	except ValueError:
		# Per valori come 'inherit', 'auto', ecc. ritorna il valore originale
		return color

base = mpl.rcParamsDefault
mpl.rcParams['font.family'] = 'Comic Sans MS'
mpl.rcParams['figure.facecolor'] = complement_color(base['figure.facecolor'])
mpl.rcParams['axes.facecolor'] = complement_color(base['axes.facecolor'])
mpl.rcParams['axes.edgecolor'] = complement_color(base['axes.edgecolor'])
mpl.rcParams['text.color'] = complement_color(base['text.color'])
mpl.rcParams['axes.labelcolor'] = complement_color(base['axes.labelcolor'])
mpl.rcParams['xtick.color'] = complement_color(base['xtick.color'])
mpl.rcParams['ytick.color'] = complement_color(base['ytick.color'])
mpl.rcParams['grid.color'] = complement_color(base['grid.color'])
mpl.rcParams['legend.facecolor'] = complement_color(base['legend.facecolor'])
mpl.rcParams['legend.edgecolor'] = complement_color(base['legend.edgecolor'])

plt.figure(figsize=(10, 6))
default_line = base['axes.prop_cycle'].by_key()['color'][0]
plt.plot(pressioni/1e5, temperature, label='Curva di Saturazione', color=complement_color(default_line))
plt.xlabel('Pressione [bar]')
plt.ylabel('Temperatura [°C]')
plt.title('Diagramma P-T Saturazione Acqua (Progetto Nucleare)')
plt.grid(True)
plt.legend()
plt.show()
