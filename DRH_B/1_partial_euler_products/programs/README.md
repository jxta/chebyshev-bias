# Programs (folder 1)

| file | purpose |
|---|---|
| `drh_scan.c` | serial scanner: segmented sieve of Eratosthenes; accumulates log ∏(1 − χ(p)p^(−1/2))^(−1) for the characters mod 3, 4, 5 (real and complex), 7 at 128 checkpoints per decade. `gcc -O2 -o drh_scan drh_scan.c -lm`, `./drh_scan 10000000000 128 > traj.csv` (10^10 in about 30 s on one core). |
| `drh_scan_mp.c` | parallel version (OpenMP; the range of x can also be split across machines). `gcc -O3 -fopenmp -o drh_scan_mp drh_scan_mp.c -lm`, `./drh_scan_mp N per node num_nodes threads`. Output is bit-for-bit independent of the number of threads. |
| `merge_traj.py` | merges the partial outputs of `drh_scan_mp` and checks π(10^k) against known values. |
| `make_tables_figures.py` | rebuilds the tables and figures of the field directories from the `.csv` files. |

Output CSV columns (`traj.csv`): `x, pi_x, S4, S3, S5, S5c_re, S5c_im, S7` where `S*` = log of the partial Euler product (for the complex character mod 5, real and imaginary parts). The product itself is `exp(S)`.
