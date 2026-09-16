# 1. Verification of the Deep Riemann Hypothesis (B) for Dirichlet L-functions

The fields in this folder are abelian extensions of **Q** and Dirichlet L-functions.

## 1.1 Setting

Let χ be a primitive Dirichlet character modulo N and

L(s, χ) = ∏_p (1 − χ(p) p^(−s))^(−1)   (Re s > 1),  with χ(p) = 0 for p | N,

its Dirichlet L-function. Put m = ord_(s=1/2) L(s, χ) (for Dirichlet characters m = 0 is conjectured (Chowla conjecture) and is true in every case treated here). The **Deep Riemann Hypothesis** (Kimura–Koyama–Kurokawa [3], Conjecture 2; see also Aoki–Koyama [1], Conjecture 1.1) asserts:

* **(A)** the limit lim_(x→∞) (log x)^m ∏_(p≤x) (1 − p^(−1/2) χ(p))^(−1) exists and is not 0;
* **(B)** the limit equals √2^(ν(χ)) · L^((m))(1/2, χ) / (e^(mγ) m!), where ν(χ) = 1 if χ² = 1 and ν(χ) = 0 otherwise, and γ is Euler's constant.

When m = 0 this reads

lim_(x→∞) ∏_(p≤x) (1 − χ(p) p^(−1/2))^(−1) = √2 · L(1/2, χ)   (χ real),   = L(1/2, χ)   (χ complex).

The Euler product is taken *on the critical line* s = 1/2, where it does not converge absolutely; DRH says that the partial products nevertheless converge, and to a value that differs from L(1/2, χ) by the factor √2 exactly when χ is real. (Convergence in the ordinary sense is known to be strictly stronger than the Riemann Hypothesis for L(s, χ); see Conrad [2] and Kimura–Koyama–Kurokawa [3].)

## 1.2 What was computed

For each character χ we computed the **partial Euler product**

P_χ(x) = ∏_(p≤x) (1 − χ(p) p^(−1/2))^(−1)

for x up to **10^15**, i.e. over the first 29,844,570,422,669 primes. In practice the logarithm log P_χ(x) = Σ_(p≤x) −log(1 − χ(p) p^(−1/2)) is accumulated prime by prime (this avoids the loss of accuracy that comes from multiplying 3·10^13 factors), and the product is recovered by exponentiating. The value is recorded at 128 checkpoints per decade of x, from x = 100 to x = 10^15; the tables show the checkpoints x = 10^5, 10^6, …, 10^15 and the figures show all of them.

The primes were generated with a segmented sieve of Eratosthenes (program `drh_scan.c`; the parallel version `drh_scan_mp.c` splits the range of x among the cores of one or several machines and reproduces the serial result exactly). The reference values L(1/2, χ) were computed independently, to 30 digits, from the Hurwitz zeta function: L(s, χ) = N^(−s) Σ_(a=1)^(N) χ(a) ζ(s, a/N).

## 1.3 How to read the tables and figures

* The tables list P_χ(x) at x = 10^5, …, 10^15 together with the conjectured limit. The values are finite products over the primes p ≤ x, computed exactly (about 13 significant digits are reliable; 12 decimals are printed so that they can be reproduced and compared). They should be read as a *sequence approaching the limit*, not as approximations of the limit to 12 digits.
* The approach to the limit is **slow and oscillating**: under DRH the relative error is assumed to decay like 1/log x, so between x = 10^5 and x = 10^15 the size of the error only shrinks by a factor of about 3, and it does not decrease monotonically. At x = 10^15 an agreement of a few parts in a thousand is what this assumption suggests. The figures show this behaviour: the curve oscillates around the dashed line (the conjectured limit) inside a slowly narrowing envelope.
* For a complex character the product is a complex number, and the table lists real and imaginary parts.

## 1.4 Checks

* π(10^k) (the number of primes ≤ 10^k) produced by the sieve agrees with the known values for every k ≤ 15; π(10^15) = 29,844,570,422,669.
* At x = 10^6 all sums agree with a direct term-by-term computation in Python to 10^(−14).
* Two independent programs (this one and the program of folder 2, which also records log P_χ(x)) agree at every checkpoint to 3·10^(−13).
* The parallel program gives results identical bit for bit to the serial program, for any number of threads and any splitting of the range.

## 1.5 Fields

| Directory | Field / characters | Status |
|---|---|---|
| [`mod4_Dirichlet_char/`](mod4_Dirichlet_char/) | **Q**(i): the character mod 4 | complete (x ≤ 10^15) |
| `mod5_Dirichlet_char/` | **Q**(ζ5): characters mod 5 | to be added |
| `mod8_Dirichlet_char/` | **Q**(ζ8): characters mod 8 | to be added |

## References

[1] M. Aoki and S. Koyama, Chebyshev's bias against splitting and principal primes in global fields, J. Number Theory 245 (2023), 233–262.

[2] K. Conrad, Partial Euler products on the critical line, Canad. J. Math. 57 (2005), 267–297.

[3] K. Kimura, S. Koyama and N. Kurokawa, Euler products beyond the boundary, Lett. Math. Phys. 104 (2014), 1–19.
