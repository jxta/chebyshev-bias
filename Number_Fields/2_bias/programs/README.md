# Programs (folder 2)

| file | purpose |
|---|---|
| `drh_scan_ak.c` | serial scanner for the Gaussian field: accumulates Σ χ(p)/√p (χ mod 4), Σ 1/p, Σ (1/p + log(1 − 1/p)), the prime-power sum Σ_(k≥3) Σ_p χ(p^k)/(k p^(k/2)) and, as a check, log ∏(1 − χ(p)p^(−1/2))^(−1). `gcc -O2 -o drh_scan_ak drh_scan_ak.c -lm`. |
| `drh_scan_ak_mp.c` | parallel version (OpenMP, splittable across machines), same output. |
| `merge_traj.py` | merges partial outputs and checks π(10^k). |
| `make_tables_figures.py` | rebuilds the tables and figures of the field directories from the `.csv` files. |

Output CSV columns: `x, pi_x, SB4, ME, M2, K3, S4chk` with
SB4 = Σ_(p≤x, p≠2) χ(p)/√p, ME = Σ_(p≤x) 1/p, M2 = Σ_(p≤x) (1/p + log(1 − 1/p)), K3 = Σ_(k≥3) Σ_(p≤x) χ(p^k)/(k p^(k/2)), S4chk = log of the partial Euler product.
For **Q**(i): π_(1/2)(x) − 2π_(1/2)(x; σ1) = 1/√2 − SB4 and π_(1/2)(x) − 2π_(1/2)(x; σ−1) = √2 − (1/√2 − SB4).
