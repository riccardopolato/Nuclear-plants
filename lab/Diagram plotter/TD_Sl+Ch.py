import pandas as pd
import matplotlib.pyplot as plt

file_csv = "TD_Sl+Ch.csv"

df = pd.read_csv(file_csv, header=None, names=["x", "y"])
df["x"] = pd.to_numeric(df["x"], errors="coerce")
df["y"] = pd.to_numeric(df["y"], errors="coerce")
df = df.dropna()

fig, ax = plt.subplots(figsize=(8, 5))

# punti digitalizzati
ax.scatter(df["x"], df["y"], color="black", s=18)

# scala dell'asse x come nella figura
ax.set_xscale("log")

# etichette assi
ax.set_xlabel(r"$\frac{J}{\sqrt{gD}}$")
ax.set_ylabel(r"$\beta$")

ax.set_title("Taitel-Dukler: Slug / Churn transition")
ax.grid(True, which="both", linestyle="--", alpha=0.35)

# scritte nelle regioni
ax.text(2e2, 0.92, "Churn flow", fontsize=18, ha="center", va="center")
ax.text(8e2, 0.74, "Slugs", fontsize=18, ha="center", va="center")

ax.set_ylim(0.5, 1.0)

plt.tight_layout()
#plt.show()