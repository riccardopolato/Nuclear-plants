import pandas as pd
import matplotlib.pyplot as plt

file_csv = "TD_An+SlCh.csv"

df = pd.read_csv(file_csv, header=None, names=["x", "y"])
df["x"] = pd.to_numeric(df["x"], errors="coerce")
df["y"] = pd.to_numeric(df["y"], errors="coerce")
df = df.dropna()

fig, ax = plt.subplots(figsize=(8, 5))

ax.scatter(df["x"], df["y"], color="black", s=18)

ax.set_xscale("log")
ax.set_yscale("log")

ax.set_xlabel(r"$X_1 \equiv \left[\dfrac{(dp/dz)_L}{(dp/dz)_G}\right]^{1/2}$")
ax.set_ylabel(r"$ \dfrac{j_g \rho_G^{1/2}}{\left[g(\rho_L-\rho_G)\sigma\right]^{1/4}}$")

ax.set_title("Taitel-Dukler: Annular / Slug-Churn transition")
ax.grid(True, which="both", linestyle="--", alpha=0.35)

ax.text(5e-1, 8e-2, "Slug\nChurn", fontsize=18, ha="center", va="center")
ax.text(3e2, 1.2, "Annular", fontsize=18, ha="center", va="center")

plt.tight_layout()
#plt.show()