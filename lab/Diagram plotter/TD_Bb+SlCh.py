import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

file_csv = Path(__file__).resolve().parent / "TD_Bb+SlCh.csv"

df = pd.read_csv(file_csv, header=None, names=["x", "y"])
df["x"] = pd.to_numeric(df["x"], errors="coerce")
df["y"] = pd.to_numeric(df["y"], errors="coerce")
df = df.dropna()

fig, ax = plt.subplots(figsize=(8, 6))

# punti digitalizzati
ax.scatter(df["x"], df["y"], color="black", s=18)

# scala logaritmica come nella figura
ax.set_xscale("log")

# etichette corrette
ax.set_xlabel(r"$ \dfrac{j_g \rho_G^{1/2}}{\left[g(\rho_L-\rho_G)\sigma\right]^{1/4}}$")
ax.set_ylabel(r"$j_L / j_G$")

ax.set_title("Taitel-Dukler: Bubbly / Slug-Churn transition")
ax.grid(True, which="both", linestyle="--", alpha=0.35)

# scritte nelle regioni
ax.text(0.9, 2.1, "Bubbles", fontsize=18, ha="center", va="center")
ax.text(8, 1.2, "Slugs\nChurn flow", fontsize=18, ha="center", va="center")

plt.tight_layout()
#plt.show()