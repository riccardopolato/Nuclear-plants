# se entrate sappiate che vi devo accettare 
import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
import numpy as np

# Definiamo un range di pressioni da 1 bar a 220 bar (vicino al punto critico)
pressioni = np.linspace(1e5, 220e5, 200)

# Calcolo della temperatura di saturazione per l'acqua
temperature_acqua = [CP.PropsSI('T', 'P', p, 'Q', 0, 'Water') - 273.15 for p in pressioni]

# Calcolo della temperatura di saturazione per l'olio (esempio: usiamo 'Oil' come fluido generico)
temperature_olio = [CP.PropsSI('T', 'P', p, 'Q', 0.5, 'Water') - 273.15 for p in pressioni]

# Creazione del grafico combinato
plt.figure(figsize=(10, 6))
plt.plot(pressioni/1e5, temperature_acqua, label='Acqua - Curva di Saturazione', color='blue')
plt.plot(pressioni/1e5, temperature_olio, label='Olio - Curva di Saturazione', color='orange')
plt.xlabel('Pressione [bar]')
plt.ylabel('Temperatura [°C]')
plt.title('Diagramma P-T Saturazione Acqua e Olio')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()