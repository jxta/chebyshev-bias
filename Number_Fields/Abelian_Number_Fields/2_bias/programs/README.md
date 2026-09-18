# Programs (folder 2)

| file | purpose |
|---|---|
| `drh_scan_ak.c` | serial scanner for the Gaussian field: accumulates Σ χ(p)/√p (χ mod 4), Σ 1/p, Σ (1/p + log(1 − 1/p)), the prime-power sum Σ_(k≥3) Σ_p χ(p^k)/(k p^(k/2)) and, as a check, log ∏(1 − χ(p)p^(−1/2))^(−1). `gcc -O2 -o drh_scan_ak drh_scan_ak.c -lm`. |
| `drh_scan_ak_mp.c` | parallel version (OpenMP, splittable across machines), same output. |
| `scan_mod58_mp.c` | parallel scanner for **Q**(ζ5) and **Q**(ζ8): accumulates, prime by prime, the class sums Σ_(p ≡ a) p^(−1/2) for a mod 5 and a mod 8, Σ 1/p, Σ (1/p + log(1 − 1/p)), Σ χ₄(p)/p (χ₄ mod 5), the prime-power sums for the characters mod 5 and mod 8, and the Euler-product logs for the characters mod 8 (used for `1_partial_euler_products/mod8_Dirichlet_char/`). `gcc -O3 -march=native -fopenmp -o scan_mod58_mp scan_mod58_mp.c -lm`; `./scan_mod58_mp N per node num_nodes threads`; the parts are merged with `merge_traj.py`. |
| `merge_traj.py` | merges partial outputs and checks π(10^k). |
| `make_tables_figures.py` | rebuilds the tables and figures of `mod4_Dirichlet_char/` from the `.csv` files. |
| `build_mod58.py` | rebuilds the tables, figures and README files of `mod5_Dirichlet_char/` and `mod8_Dirichlet_char/` (in both folders) from the merged scan output `mod58_1e15.csv` and the mod-5 Euler-product logs; uses `minisvg.py`. |
| `minisvg.py` | compact SVG plotting used for the figures. |

Output CSV columns of `drh_scan_ak*.c`: `x, pi_x, SB4, ME, M2, K3, S4chk` with
SB4 = Σ_(p≤x, p≠2) χ(p)/√p, ME = Σ_(p≤x) 1/p, M2 = Σ_(p≤x) (1/p + log(1 − 1/p)), K3 = Σ_(k≥3) Σ_(p≤x) χ(p^k)/(k p^(k/2)), S4chk = log of the partial Euler product.
For **Q**(i): π_(1/2)(x) − 2π_(1/2)(x; σ1) = 1/√2 − SB4 and π_(1/2)(x) − 2π_(1/2)(x; σ−1) = √2 − (1/√2 − SB4).

Output CSV columns of `scan_mod58_mp.c`: `x, pi_x, P5_1, P5_2, P5_3, P5_4, P8_1, P8_3, P8_5, P8_7, ME, M2, H4, K3_5, K3_2re, K3_2im, K3_m4, K3_m8, K3_8, S_m8, S_8, S_m4` with P5_a = Σ_(p≤x, p≡a mod 5) p^(−1/2), P8_a likewise mod 8, H4 = Σ χ₄(p)/p (mod 5), K3_* the prime-power sums Σ_(k≥3) Σ_p χ(p^k)/(k p^(k/2)) (real and imaginary parts for the complex χ₂ mod 5), S_* = log of the partial Euler products for χ₃, χ₄ mod 8 and for the character mod 4.
