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
# vorrei runnare ma non posso
# riesci a scrivere ludo?
# dove si runna?

# sei da browser?
# si, devo scaricare qualcosa?
# non lo so, sto ancora cercando di capire come funzioni, ho paura che solo l'host possa 

# la AI ti funziona da browser?
#da cosa lo capisco
# che prova a scrivere lui per te in continuazione, quindi credo di no
#no anzi mi va su ogni volta che scrivo, sembra abbia crisi epilettich
# perfetto, ora provo a smanettare un po per capire come migliorare sta roba
# , per ora ho solo scaricato la roba dell'acqua che dice nelle slide il
#