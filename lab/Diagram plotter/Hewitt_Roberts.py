import pandas as pd
import matplotlib.pyplot as plt

file_csv = "lab/Diagram Plotter/Hewitt-Roberts.csv"

df = pd.read_csv(file_csv, header=None, names=["x", "y"])
df["x"] = pd.to_numeric(df["x"], errors="coerce")
df["y"] = pd.to_numeric(df["y"], errors="coerce")
df = df.dropna()

fig, ax = plt.subplots(figsize=(8, 7))

# punti del diagramma
ax.scatter(df["x"], df["y"], color="black", s=18)

# scale logaritmiche
ax.set_xscale("log")
ax.set_yscale("log")

# etichette assi e titolo
ax.set_xlabel(r"$\rho_L j_L^2$ [m/s]")
ax.set_ylabel(r"$\rho_G j_G^2$ [m/s]")
ax.set_title("Hewitt-Roberts map")

# griglia
ax.grid(True, which="both", linestyle="--", alpha=0.35)

# scritte nelle regioni
ax.text(2e1, 2e3, "Annular", fontsize=13)
ax.text(7e3, 2e3, "Wispy-annular", fontsize=13)
ax.text(2e1, 4e1, "Churn", fontsize=13)
ax.text(1.2e4, 2e1, "Bubbly", fontsize=13)
ax.text(1e3, 1.2, "Bubbles\nslugs", fontsize=13, ha="center")
ax.text(2e1, 0.35, "Slugs", fontsize=13)

plt.tight_layout()
#plt.show()
