# se entrate sappiate che vi devo accettare 
import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
import numpy as np

# Definiamo un range di pressioni da 1 bar a 220 bar (vicino al punto critico)
pressioni = np.linspace(1e5, 220e5, 100)
temperature = [CP.PropsSI('T', 'P', p, 'Q', 0, 'Water') - 273.15 for p in pressioni]

plt.figure(figsize=(10, 6))
plt.plot(pressioni/1e5, temperature, label='Curva di Saturazione')
plt.xlabel('Pressione [bar]')
plt.ylabel('Temperatura [°C]')
plt.title('Diagramma P-T Saturazione Acqua (Progetto Nucleare)')
plt.grid(True)
plt.legend()
plt.show()

a = 0.5  # Coefficiente di attrito
L = 100  # Lunghezza del tubo in metri
D = 0.1  # Diametro del tubo in metri