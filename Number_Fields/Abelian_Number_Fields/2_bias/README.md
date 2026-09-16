# 2. Verification of Chebyshev's bias: main term and constant term

## 2.1 Setting

Let L/**Q** be a finite abelian extension with Galois group G, and for a prime p unramified in L let Frob_p ∈ G be its Frobenius element. Following Aoki–Koyama [1] put

π_(1/2)(x) = Σ_(p≤x) p^(−1/2),   π_(1/2)(x; σ) = Σ_(p≤x, Frob_p = σ) p^(−1/2)   (σ ∈ G).

**Main term (Aoki–Koyama, 2023).** Theorem 2.2 of [1] states that under the Deep Riemann Hypothesis, for every σ ∈ G there is a constant c such that

π_(1/2)(x) − |G| π_(1/2)(x; σ) = (M(σ) + m(σ)) log log x + c + o(1)   (x → ∞),

where M(σ) = ½ Σ_(ρ≠1) χ_ρ(σ) ν(ρ) is computed from the character table of G and m(σ) from the orders of vanishing of the L-functions at s = 1/2 (m(σ) = 0 in the fields treated here). The primes with Frob_p = σ appear "late" when M(σ) + m(σ) > 0 and "early" when it is negative.

**Constant term (Aoki, 2026).** For abelian G, all representations are irreducible and 1-dimensional.

**π_(1/2)(x) − |G| π_(1/2)(x; σ) = (M(σ) + m(σ)) log log x + c + o(1),** and the constant is

**c = (M(σ) + m(σ)) γ + R − M(σ)(log 2 + c_Q) − Σ_(ρ≠1) ρ̄(σ) ( log( L^((m))(1/2, ρ) / m! ) − c(ρ) ),**

and this statement (for all σ ∈ G) is *equivalent* to DRH (A) and (B) for all representations ρ ≠ 1 of G. The notation:

* γ = 0.5772156649… is Euler's constant; R = Σ_(p | D_L) p^(−1/2), the sum over the primes ramified in L (the ramified primes belong to no class σ);
* **c_Q** = − Σ_p ( 1/p + log(1 − 1/p) ) = γ − M_Mertens = 0.3157184520539…, where M_Mertens = 0.2614972128… is the Meissel–Mertens constant; this constant does not depend on the field;
* m = m(ρ) = ord_(s=1/2) L(s, ρ), and L^((m))(1/2, ρ) is the m-th derivative at s = 1/2;
* **c(ρ)** is determined by ρ; this includes infinite sums which converge absolutely;
* ρ̄(σ) is the complex conjugate of ρ(σ).

The theoretical constant is therefore

**c(σ) := (M(σ) + m(σ)) γ + R − M(σ)(log 2 + c_Q) − Σ_(ρ≠1) ρ̄(σ) ( log( L^((m))(1/2, ρ) / m! ) − c(ρ) ),**

known to many digits, and the theorem says that the **data constant**

**D_σ(x) := π_(1/2)(x) − |G| π_(1/2)(x; σ) − (M(σ) + m(σ)) log log x**

converges to c(σ) as x → ∞.

## 2.2 What was computed

For each field we computed, for all primes up to **10^15**, the class sums π_(1/2)(x; σ) for every σ ∈ G, and the sums entering c_Q and c(ρ) (which are already stable to 10^(−11) far below 10^15). For every σ we tabulate

* the **left-hand side** π_(1/2)(x) − |G| π_(1/2)(x; σ), and
* the **data constant** D_σ(x) = left-hand side − (M(σ) + m(σ)) log log x,

at x = 10^5, …, 10^15, and give two figures:

1. **the trajectories** π_(1/2)(x) − |G| π_(1/2)(x; σ) for all σ, together with the theoretical law (M(σ) + m(σ)) log log x + c(σ) drawn with the *theoretical* constant c(σ) (no fitting) — this displays the main term;
2. **the data constants** D_σ(x) for all σ, together with the theoretical constants c(σ) as horizontal lines — this displays the constant term.

Each field's README also lists the values of L(1/2, ρ), the orders of vanishing m(ρ) = ord_(s=1/2) L(s, ρ), and the lowest non-trivial zero of each L(s, ρ) (its imaginary part γ1), which governs the slowest oscillation of the data around the law. Values of L(1/2, ρ) and of the zeros were computed to 30 digits from the Hurwitz zeta function; c_Q was taken from the known constants γ and M_Mertens (our own prime sum reproduces it to 13 digits).

## 2.3 How to read the tables and figures

* The left-hand side grows like (M(σ) + m(σ)) log log x; in figure 1 the theoretical law (dashed) runs through the data without any adjustable parameter. The data constant removes this growth and should settle at c(σ) (figure 2).
* The convergence is slow: empirically, D_σ(x) − c(σ) oscillates with an envelope of size about 1.2/log x (± 0.035 at x = 10^15). So the data constants agree with c(σ) to two or three decimals, while c(σ) itself is known to ten. This is the expected behaviour, not a deficiency of the computation: the oscillation is expected to come from the low-lying zeros of L(s, ρ) on the critical line — the leading oscillation has period 2π/γ1 in log x.
* Since Σ_σ π_(1/2)(x; σ) = π_(1/2)(x) − R, the left-hand sides for the different σ are not independent; their sum is |G|·R for every x, and Σ_σ c(σ) = |G|·R. For a quadratic field the two data constants are mirror images of each other around R.

## 2.4 Fields

| Directory | Field / group | Status |
|---|---|---|
| [`mod4_Dirichlet_char/`](mod4_Dirichlet_char/) | **Q**(i), G = {σ1, σ−1} | complete (x ≤ 10^15) |
| `mod5_Dirichlet_char/` | **Q**(ζ5), G = {σ1, σ2, σ3, σ4} | to be added |
| `mod8_Dirichlet_char/` | **Q**(ζ8), G = {σ1, σ3, σ5, σ7} | to be added |

## Reference

[1] M. Aoki and S. Koyama, Chebyshev's bias against splitting and principal primes in global fields, J. Number Theory 245 (2023), 233–262.
