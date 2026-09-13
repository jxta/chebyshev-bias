"""Rebuild the table and figure of mod4_Dirichlet_char from partial_euler_product_mod4.csv."""
import math, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
d = pd.read_csv("../mod4_Dirichlet_char/partial_euler_product_mod4.csv")
x = d.x.values.astype(float); P = d.partial_euler_product.values.astype(float)
target = 0.944258314238200
print("| x | partial Euler product |\n|---|---|")
for k in range(5, 16): print(f"| 10^{k} | {P[np.where(x==10.0**k)[0][0]]:.12f} |")
m = x >= 1e3
fig, ax = plt.subplots(figsize=(7.4, 4.3))
ax.semilogx(x[m], P[m], color="k", lw=0.9); ax.axhline(target, color="0.45", ls="--")
ax.set_xlabel("x"); ax.set_ylabel("partial Euler product"); fig.tight_layout(); fig.savefig("fig_mod4.png")
