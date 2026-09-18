# The characters mod 5

The Dirichlet characters modulo 5 are:

| a mod 5 | 0 | 1 | 2 | 3 | 4 | conductor | |
|---|---|---|---|---|---|---|---|
| χ₁(a) | 0 | 1 | 1 | 1 | 1 | 1 | trivial, not primitive |
| χ₂(a) | 0 | 1 | i | −i | −1 | 5 | primitive, order 4 |
| χ₃(a) = conj χ₂(a) | 0 | 1 | −i | i | −1 | 5 | primitive, order 4 |
| χ₄(a) = (a/5) | 0 | 1 | −1 | −1 | 1 | 5 | primitive, real (quadratic) |

For i = 2, 3, 4 one has m = ord_(s=1/2) L(s, χᵢ) = 0. Since χ₂² = χ₃² = χ₄ ≠ 1 and χ₄² = 1, we have ν(χ₂) = ν(χ₃) = 0 and ν(χ₄) = 1, so DRH (B) states

**lim_(x→∞) ∏_(p≤x) (1 − p^(−1/2) χᵢ(p))^(−1) = L(1/2, χᵢ)  (i = 2, 3),  = √2 · L(1/2, χ₄)  (i = 4).**

The conjectured limits are

| character | limit |
|---|---|
| χ₂ | L(1/2, χ₂) = 0.763747880117 + 0.216964767519 i |
| χ₃ | L(1/2, χ₃) = 0.763747880117 − 0.216964767519 i  (complex conjugate of the above) |
| χ₄ | √2 · L(1/2, χ₄) = 0.327745333052994  (L(1/2, χ₄) = 0.231750947504016) |

Note that for the complex characters χ₂, χ₃ there is **no factor √2**, while for the real character χ₄ there is. The partial products for χ₃ are the complex conjugates of those for χ₂ at every x (because χ₃(p) = conj χ₂(p)), so they carry the same information.

## Table: partial Euler products ∏_(p≤x) (1 − p^(−1/2) χᵢ(p))^(−1)

| x | χ₂ | χ₃ | χ₄ |
|---|---|---|---|
| 10^5 | 0.745059925 +0.238382275 i | 0.745059925 -0.238382275 i | 0.326547073261 |
| 10^6 | 0.791488633 +0.186397429 i | 0.791488633 -0.186397429 i | 0.330797485580 |
| 10^7 | 0.781746548 +0.216994407 i | 0.781746548 -0.216994407 i | 0.320619903626 |
| 10^8 | 0.773259603 +0.220921370 i | 0.773259603 -0.220921370 i | 0.332076753951 |
| 10^9 | 0.765529242 +0.213361130 i | 0.765529242 -0.213361130 i | 0.323971830472 |
| 10^10 | 0.768953108 +0.204042673 i | 0.768953108 -0.204042673 i | 0.321140911486 |
| 10^11 | 0.774842607 +0.244666981 i | 0.774842607 -0.244666981 i | 0.341454396346 |
| 10^12 | 0.750290467 +0.227003726 i | 0.750290467 -0.227003726 i | 0.326875113352 |
| 10^13 | 0.766680240 +0.203144322 i | 0.766680240 -0.203144322 i | 0.333938501336 |
| 10^14 | 0.757637016 +0.214064963 i | 0.757637016 -0.214064963 i | 0.315760239824 |
| 10^15 | 0.767629772 +0.214339108 i | 0.767629772 -0.214339108 i | 0.325443055428 |

(The full data, 128 checkpoints per decade from x = 100 to 10^15, is in `partial_euler_products_mod5.csv`: columns `x`, `pi_x`, `product_chi2_re`, `product_chi2_im`, `product_chi4`.)

## Figure

![partial Euler products, characters mod 5](fig_mod5.svg)

Top: the real character χ₄, with the limit √2 L(1/2, χ₄) (dashed). Middle and bottom: real and imaginary parts of the product for the complex character χ₂, with Re L(1/2, χ₂) and Im L(1/2, χ₂) (dashed). All three curves oscillate around their limits inside envelopes that narrow like 1/log x; at x = 10^15 the χ₄ product differs from its limit by 0.70 % and the χ₂ product by 0.59 % (in absolute value). If the factor √2 were wrongly applied to χ₂, the discrepancy would be about 29 % throughout, so the data distinguish the two cases ν = 0 and ν = 1 clearly.
