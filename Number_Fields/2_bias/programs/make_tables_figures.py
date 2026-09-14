"""Rebuild the table and figure of mod4_Dirichlet_char from bias_mod4.csv."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
d = pd.read_csv("../mod4_Dirichlet_char/bias_mod4.csv")
x = d.x.values.astype(float)
c1, c2 = 0.6355768227, 0.7786367396
print("| x | LHS s1 | data const s1 | LHS s-1 | data const s-1 |\n|---|---|---|---|---|")
for k in range(5, 16):
    r = d[d.x == 10**k].iloc[0]
    print(f"| 10^{k} | {r.lhs_sigma1:.12f} | {r.data_constant_sigma1:.12f} | {r.lhs_sigma_minus1:.12f} | {r.data_constant_sigma_minus1:.12f} |")
m = x >= 1e3
fig, ax = plt.subplots(figsize=(7.4, 4.6))
ax.semilogx(x[m], d.data_constant_sigma_minus1.values[m], color="0.45", lw=0.9)
ax.semilogx(x[m], d.data_constant_sigma1.values[m], color="k", lw=0.9)
ax.axhline(c1, color="k", ls="--"); ax.axhline(c2, color="0.45", ls="--")
ax.set_xlabel("x"); ax.set_ylabel("data constant"); fig.tight_layout(); fig.savefig("fig_constants_mod4.png")
