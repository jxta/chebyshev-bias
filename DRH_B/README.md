# Number_Fields — the Deep Riemann Hypothesis (A), (B) and Chebyshev's bias: numerical verification

This part of the repository accompanies the work of **Miho Aoki** (Shimane University) and **Shin-ya Koyama** (Toyo University) on Chebyshev's bias for Galois extensions of number fields and on the Deep Riemann Hypothesis (A), (B), and contains the large-scale computations (primes up to 10^15) carried out by **Shigetoshi Yokoyama** (National Institute of Informatics). It is written so that a reader who does not use computers can follow everything except the program files.

## What is here

| Folder | Content |
|---|---|
| [`1_partial_euler_products/`](1_partial_euler_products/) | Verification of the Deep Riemann Hypothesis DRH (A), (B) for Artin L-functions: partial Euler products on the critical line, calculation results for all primes up to 10^15, verification of DRH (A) (convergence) and DRH (B) (convergence value). |
| [`2_bias/`](2_bias/) | Verification of Chebyshev's bias: the main term (M(σ) + m(σ)) log log x (Aoki–Koyama, 2023) and the explicit **constant term** (Aoki, 2026) involving L(1/2, χ), the Meissel–Mertens constant and prime-power sums. |

Each folder has its own `README.md` (the mathematical setting, what was computed, how to read the tables and figures), a `programs/` directory (the source code that produced the data) and one directory per field:

* `mod4_Dirichlet_char/` — the Gaussian field **Q**(i), character mod 4
* `mod5_Dirichlet_char/` — the field **Q**(ζ5), characters mod 5 (to be added)
* `mod8_Dirichlet_char/` — the field **Q**(ζ8), characters mod 8 (to be added)

Every such directory contains a short `README.md` with the tables and figures, the figure files (`.svg`, `.png`) and `.csv` files with the full data (128 checkpoints per decade of x).

## The underlying statements (Case of abelian extensions over Q)

* **DRH (A), (B)** (Kimura–Koyama–Kurokawa): for a primitive Dirichlet character χ with m = ord_(s=1/2) L(s, χ), (A) the limit lim_(x→∞) (log x)^m ∏_(p≤x) (1 − χ(p) p^(−1/2))^(−1) exists and is not 0, and (B) it equals √2^(ν(χ)) · L^((m))(1/2, χ) / (e^(mγ) m!), where ν(χ) = 1 if χ² = 1 and 0 otherwise.
* **Chebyshev's bias (Aoki–Koyama, 2023) and with explicit constant (Aoki, 2026)**: for an abelian extension L/**Q** with Galois group G and σ ∈ G, under DRH (A),
  π_(1/2)(x) − |G| π_(1/2)(x; σ) = (M(σ) + m(σ)) log log x + c + o(1)   (Theorem 2.2 of [1]),
  and the constant is
  c = (M(σ) + m(σ)) γ + Σ_(p | D_L) p^(−1/2) − M(σ)(log 2 + c_**Q**) − Σ_(χ≠1) χ̄(σ) ( log( L^((m))(1/2, χ) / m! ) − c(χ) )   (Aoki, 2026)
  under DRH (A), (B). See `2_bias/README.md` for the notation.

## How to read the numbers

The tabulated data — the partial Euler products and the sums π_(1/2)(x) − |G| π_(1/2)(x; σ) — are finite sums over the primes p ≤ x; they are computed exactly (about 13 significant digits are reliable) and printed with 12 decimals so that anyone can reproduce and compare them. The theoretical quantities they are compared with — L(1/2, χ), the constants c_**Q** and c(χ), which are given by convergent infinite sums or special values — are computed independently to high precision. The **convergence** of the data to the theoretical values is slow — the error decays only like 1/log x and oscillates — so the agreement one should expect at x = 10^15 is a few parts in a thousand, not twelve digits. The figures show this: the curves oscillate around the dashed theoretical values with a slowly narrowing envelope.

## Reproducibility and checks

The prime sums were computed with a segmented sieve of Eratosthenes in C, with compensated (Kahan) summation, on the mdx research cloud (one node, 152 threads, about 28 hours for x = 10^15). The prime-counting function π(10^k) was checked against the known values for every k ≤ 15 (π(10^15) = 29,844,570,422,669); the sums were checked against brute force at x = 10^6 and two independent implementations agree to 3·10^(−13) at every checkpoint. Details are given in the folder READMEs.

## Licenses

* Text, tables, figures and data files: [CC BY 4.0](../LICENSE-DATA.md).
* Program source code: [MIT License](../LICENSE-CODE.md).

## References

[1] M. Aoki and S. Koyama, *Chebyshev's bias against splitting and principal primes in global fields*, J. Number Theory 245 (2023), 233–262.

[2] K. Conrad, *Partial Euler products on the critical line*, Canad. J. Math. 57 (2005), 267–297.

[3] K. Kimura, S. Koyama and N. Kurokawa, *Euler products beyond the boundary*, Lett. Math. Phys. 104 (2014), 1–19.
