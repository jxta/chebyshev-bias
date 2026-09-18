"""build_mod58.py — builds Number_Fields/{1_partial_euler_products,2_bias}/{mod5,mod8}_Dirichlet_char from
mod58_1e15.csv (class sums, Euler logs, auxiliary sums) and drh_traj_1e15.csv (mod-5 Euler logs). Output under out/."""
import math, os, re, numpy as np, pandas as pd, matplotlib
SUB = {"1":"₁","2":"₂","3":"₃","4":"₄","5":"₅","6":"₆","7":"₇","8":"₈","9":"₉"}
def subs(t):
    t = re.sub(r"([χρ])([1-9])", lambda m: m.group(1)+SUB[m.group(2)], t)
    t = t.replace("χi²", "χᵢ²")
    return re.sub(r"χi\b", "χᵢ", t)
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from minisvg import figure as _figure
def figure(panels, xlo, xhi, path):
    n = _figure(panels, xlo, xhi, path)
    t = open(path).read(); open(path, "w").write(subs(t)); return n
from mpmath import mp, zeta, mpf, mpc, findroot, fabs, euler, mertens
mp.dps = 30
g, M = float(euler), float(mertens); cQ = g - M; s2 = math.sqrt(2)
def Lf(q, chi, s): return sum(chi[a]*zeta(s, mpf(a)/q) for a in chi) * mp.power(q, -s)
def L(q, chi): return complex(Lf(q, chi, mpf(1)/2))
def lowest_zero(q, chi, both=False, tmax=15.0):
    f = lambda s: Lf(q, chi, s); best = None
    ts = np.arange(-tmax if both else 0.3, tmax, 0.05)
    v = np.array([float(fabs(f(mpc(0.5, t)))) for t in ts])
    for i in range(1, len(ts)-1):
        if v[i] < v[i-1] and v[i] < v[i+1] and v[i] < 0.3:
            r = findroot(f, mpc(0.5, ts[i]), solver="muller")
            if abs(float(r.real)-0.5) < 1e-10 and float(fabs(f(r))) < 1e-15:
                gz = float(r.imag)
                if best is None or abs(gz) < abs(best): best = gz
    return best

d = pd.read_csv("mod58_1e15.csv"); drh = pd.read_csv("drh_traj_1e15.csv")
assert (d.x.values == drh.x.values).all()
x = d.x.values.astype(float); lx = np.log(x); llx = np.log(lx)
sel = (np.arange(len(x)) % 4) == 0; xs = x[sel]; m = x >= 1e3
ks = list(range(5, 16)); idx = {k: np.where(x == 10.0**k)[0][0] for k in ks}
plt.rcParams.update({"font.size": 9.5, "axes.grid": True, "grid.alpha": 0.25, "figure.dpi": 130, "axes.axisbelow": True})
def png(panels, path, fs, legend_loc="upper right"):
    fig, axs = plt.subplots(len(panels), 1, figsize=fs, sharex=True); axs = np.atleast_1d(axs)
    for ax, p in zip(axs, panels):
        for (cx, cy, col, lw, lab) in p["curves"]: ax.semilogx(cx, cy, color=col, lw=lw, label=subs(lab), ls=p.get("ls", {}).get(lab, "-"))
        for (hv, col, lab) in p.get("hlines", []): ax.axhline(hv, color=col, ls="--", lw=1.2, label=subs(lab))
        ax.set_ylabel(subs(p.get("ylabel", ""))); ax.legend(loc=legend_loc, fontsize=7.8)
        if p.get("title"): ax.set_title(subs(p["title"]), fontsize=10)
    axs[-1].set_xlabel("x"); fig.tight_layout(); fig.savefig(path, bbox_inches="tight"); plt.close(fig)
O = "out/Number_Fields"
for sub in ("1_partial_euler_products/mod5_Dirichlet_char", "1_partial_euler_products/mod8_Dirichlet_char",
            "2_bias/mod5_Dirichlet_char", "2_bias/mod8_Dirichlet_char"): os.makedirs(f"{O}/{sub}", exist_ok=True)

# ===================== characters and constants =====================
chi52 = {1: 1, 2: 1j, 3: -1j, 4: -1}; chi54 = {1: 1, 2: -1, 3: -1, 4: 1}
L52 = L(5, chi52); L54 = L(5, chi54).real
ch8 = {2: {1: 1, 3: -1, 5: 1, 7: -1}, 3: {1: 1, 3: 1, 5: -1, 7: -1}, 4: {1: 1, 3: -1, 5: -1, 7: 1}}
L8 = {j: L(8, ch8[j]).real for j in (2, 3, 4)}
z52 = lowest_zero(5, chi52, both=True); z54 = lowest_zero(5, chi54); z8 = {j: lowest_zero(8, ch8[j]) for j in (2, 3, 4)}
print("lowest zeros: chi2 mod5", z52, " chi4 mod5", z54, " mod8", z8)

# ===================== 1. partial Euler products =====================
P54 = np.exp(drh.S5.values); P52 = np.exp(drh.S5c_re.values) * np.exp(1j*drh.S5c_im.values)
E5 = f"{O}/1_partial_euler_products/mod5_Dirichlet_char"
pd.DataFrame({"x": d.x.values, "pi_x": d.pi_x.values, "product_chi2_re": [f"{v:.15g}" for v in P52.real],
              "product_chi2_im": [f"{v:.15g}" for v in P52.imag], "product_chi4": [f"{v:.15g}" for v in P54]}).to_csv(f"{E5}/partial_euler_products_mod5.csv", index=False)
pan = [dict(curves=[(xs, P54[sel], "#000", 0.85, "partial Euler product, χ4 (real)")], hlines=[(s2*L54, "#666", f"√2·L(1/2,χ4) = {s2*L54:.9f}…")], ylabel="χ4",
            title="Characters mod 5 (Q(ζ5)): partial Euler products up to x = 10^15"),
       dict(curves=[(xs, P52.real[sel], "#000", 0.85, "Re of partial Euler product, χ2")], hlines=[(L52.real, "#666", f"Re L(1/2,χ2) = {L52.real:.9f}…")], ylabel="χ2: real part"),
       dict(curves=[(xs, P52.imag[sel], "#000", 0.85, "Im of partial Euler product, χ2")], hlines=[(L52.imag, "#666", f"Im L(1/2,χ2) = {L52.imag:.9f}…")], ylabel="χ2: imaginary part")]
figure(pan, 1e3, 1e15, f"{E5}/fig_mod5.svg"); png([dict(p, curves=[(x[m], c[1] if False else None, 0, 0, "")]) for p, c in []] or
    [dict(curves=[(x[m], P54[m], "#000", 0.85, "partial Euler product, χ4 (real)")], hlines=[(s2*L54, "#666", f"√2·L(1/2,χ4) = {s2*L54:.9f}…")], ylabel="χ4", title="Characters mod 5 (Q(ζ5)): partial Euler products up to x = 10^15"),
     dict(curves=[(x[m], P52.real[m], "#000", 0.85, "Re of partial Euler product, χ2")], hlines=[(L52.real, "#666", f"Re L(1/2,χ2) = {L52.real:.9f}…")], ylabel="χ2: real part"),
     dict(curves=[(x[m], P52.imag[m], "#000", 0.85, "Im of partial Euler product, χ2")], hlines=[(L52.imag, "#666", f"Im L(1/2,χ2) = {L52.imag:.9f}…")], ylabel="χ2: imaginary part")], f"{E5}/fig_mod5.png", (7.4, 8.6))
rows = "\n".join(f"| 10^{k} | {P52[idx[k]].real:.9f} {P52[idx[k]].imag:+.9f} i | {P52[idx[k]].real:.9f} {-P52[idx[k]].imag:+.9f} i | {P54[idx[k]]:.12f} |" for k in ks)
open(f"{E5}/README.md", "w").write(subs(f"""# The characters mod 5

The Dirichlet characters modulo 5 are:

| a mod 5 | 0 | 1 | 2 | 3 | 4 | conductor | |
|---|---|---|---|---|---|---|---|
| χ1(a) | 0 | 1 | 1 | 1 | 1 | 1 | trivial, not primitive |
| χ2(a) | 0 | 1 | i | −i | −1 | 5 | primitive, order 4 |
| χ3(a) = conj χ2(a) | 0 | 1 | −i | i | −1 | 5 | primitive, order 4 |
| χ4(a) = (a/5) | 0 | 1 | −1 | −1 | 1 | 5 | primitive, real (quadratic) |

For i = 2, 3, 4 one has m = ord_(s=1/2) L(s, χi) = 0. Since χ2² = χ3² = χ4 ≠ 1 and χ4² = 1, we have ν(χ2) = ν(χ3) = 0 and ν(χ4) = 1, so DRH (B) states

**lim_(x→∞) ∏_(p≤x) (1 − p^(−1/2) χi(p))^(−1) = L(1/2, χi)  (i = 2, 3),  = √2 · L(1/2, χ4)  (i = 4).**

The conjectured limits are

| character | limit |
|---|---|
| χ2 | L(1/2, χ2) = {L52.real:.12f} + {L52.imag:.12f} i |
| χ3 | L(1/2, χ3) = {L52.real:.12f} − {L52.imag:.12f} i  (complex conjugate of the above) |
| χ4 | √2 · L(1/2, χ4) = {s2*L54:.15f}  (L(1/2, χ4) = {L54:.15f}) |

Note that for the complex characters χ2, χ3 there is **no factor √2**, while for the real character χ4 there is. The partial products for χ3 are the complex conjugates of those for χ2 at every x (because χ3(p) = conj χ2(p)), so they carry the same information.

## Table: partial Euler products ∏_(p≤x) (1 − p^(−1/2) χi(p))^(−1)

| x | χ2 | χ3 | χ4 |
|---|---|---|---|
{rows}

(The full data, 128 checkpoints per decade from x = 100 to 10^15, is in `partial_euler_products_mod5.csv`: columns `x`, `pi_x`, `product_chi2_re`, `product_chi2_im`, `product_chi4`.)

## Figure

![partial Euler products, characters mod 5](fig_mod5.svg)

Top: the real character χ4, with the limit √2 L(1/2, χ4) (dashed). Middle and bottom: real and imaginary parts of the product for the complex character χ2, with Re L(1/2, χ2) and Im L(1/2, χ2) (dashed). All three curves oscillate around their limits inside envelopes that narrow like 1/log x; at x = 10^15 the χ4 product differs from its limit by {100*abs(P54[idx[15]]/(s2*L54)-1):.2f} % and the χ2 product by {100*abs(P52[idx[15]]-L52)/abs(L52):.2f} % (in absolute value). If the factor √2 were wrongly applied to χ2, the discrepancy would be about 29 % throughout, so the data distinguish the two cases ν = 0 and ν = 1 clearly.
"""))

E8 = f"{O}/1_partial_euler_products/mod8_Dirichlet_char"
P8 = {2: np.exp(d.S_m4.values), 3: np.exp(d.S_m8.values), 4: np.exp(d.S_8.values)}; T8 = {j: s2*L8[j] for j in (2,3,4)}
pd.DataFrame({"x": d.x.values, "pi_x": d.pi_x.values, **{f"product_chi{j}": [f"{v:.15g}" for v in P8[j]] for j in (2,3,4)}}).to_csv(f"{E8}/partial_euler_products_mod8.csv", index=False)
lab8 = {2: "χ2 (conductor 4)", 3: "χ3 (conductor 8)", 4: "χ4 (conductor 8)"}
pan = [dict(curves=[(xs, P8[j][sel], "#000", 0.85, f"partial Euler product, {lab8[j]}")], hlines=[(T8[j], "#666", f"√2·L(1/2,χ{j}) = {T8[j]:.9f}…")], ylabel=f"χ{j}",
            title=("Characters mod 8 (Q(ζ8)): partial Euler products up to x = 10^15" if j == 2 else None)) for j in (2,3,4)]
figure(pan, 1e3, 1e15, f"{E8}/fig_mod8.svg")
png([dict(curves=[(x[m], P8[j][m], "#000", 0.85, f"partial Euler product, {lab8[j]}")], hlines=[(T8[j], "#666", f"√2·L(1/2,χ{j}) = {T8[j]:.9f}…")], ylabel=f"χ{j}",
          title=("Characters mod 8 (Q(ζ8)): partial Euler products up to x = 10^15" if j == 2 else None)) for j in (2,3,4)], f"{E8}/fig_mod8.png", (7.4, 8.6))
rows = "\n".join(f"| 10^{k} | {P8[2][idx[k]]:.12f} | {P8[3][idx[k]]:.12f} | {P8[4][idx[k]]:.12f} |" for k in ks)
open(f"{E8}/README.md", "w").write(subs(f"""# The characters mod 8

The Dirichlet characters modulo 8 are:

| a mod 8 | 1 | 3 | 5 | 7 | conductor | |
|---|---|---|---|---|---|---|
| χ1(a) | 1 | 1 | 1 | 1 | 1 | trivial, not primitive |
| χ2(a) | 1 | −1 | 1 | −1 | 4 | not primitive mod 8 (it is the character mod 4 of `mod4_Dirichlet_char/`) |
| χ3(a) | 1 | 1 | −1 | −1 | 8 | primitive; χ3(p) = (−2/p), the character of **Q**(√−2) |
| χ4(a) | 1 | −1 | −1 | 1 | 8 | primitive; χ4(p) = (2/p), the character of **Q**(√2) |

(χi(a) = 0 for even a.) For i = 2, 3, 4 one has m = ord_(s=1/2) L(s, χi) = 0 and χi² = 1, so ν(χi) = 1 and DRH (B) states

**lim_(x→∞) ∏_(p≤x) (1 − p^(−1/2) χi(p))^(−1) = √2 · L(1/2, χi)  (i = 2, 3, 4).**

| character | L(1/2, χi) | limit √2 · L(1/2, χi) |
|---|---|---|
| χ2 | {L8[2]:.15f} | {T8[2]:.15f} |
| χ3 | {L8[3]:.15f} | {T8[3]:.15f} |
| χ4 | {L8[4]:.15f} | {T8[4]:.15f} |

Since χ2(p) coincides with the character mod 4 for every prime p, its partial Euler products are exactly the numbers of `mod4_Dirichlet_char/`.

## Table: partial Euler products ∏_(p≤x) (1 − p^(−1/2) χi(p))^(−1)

| x | χ2 | χ3 | χ4 |
|---|---|---|---|
{rows}

(The full data, 128 checkpoints per decade from x = 100 to 10^15, is in `partial_euler_products_mod8.csv`: columns `x`, `pi_x`, `product_chi2`, `product_chi3`, `product_chi4`.)

## Figure

![partial Euler products, characters mod 8](fig_mod8.svg)

Each panel shows the partial Euler product of one character with its conjectured limit √2 L(1/2, χi) (dashed). The curves oscillate around the limits inside envelopes that narrow like 1/log x; at x = 10^15 the products differ from their limits by {100*abs(P8[2][idx[15]]/T8[2]-1):.2f} % (χ2), {100*abs(P8[3][idx[15]]/T8[3]-1):.2f} % (χ3) and {100*abs(P8[4][idx[15]]/T8[4]-1):.2f} % (χ4).
"""))

# ===================== 2. bias =====================
def bias_page(name, field, Rq, ram, Mfun, classes, LHS, Csig, consts_md, theorem_md, Gdesc, fig_titles, csvname, extra_md=""):
    B = f"{O}/2_bias/{name}"
    DC = {a: LHS[a] - Mfun[a]*llx for a in classes}
    cols = {"x": d.x.values, "pi_x": d.pi_x.values}
    for a in classes: cols[f"lhs_sigma{a}"] = [f"{v:.15g}" for v in LHS[a]]; cols[f"data_constant_sigma{a}"] = [f"{v:.15g}" for v in DC[a]]
    pd.DataFrame(cols).to_csv(f"{B}/{csvname}", index=False)
    colr = dict(zip(classes, ["#000", "#c00000", "#0050c0", "#777"]))
    # fig 1: trajectories + laws
    curves = []
    for a in classes:
        law = Mfun[a]*llx + Csig[a]
        curves.append((xs, LHS[a][sel], colr[a], 0.85, f"σ{a}: π½(x) − 4π½(x;σ{a})"))
        curves.append((xs, law[sel], colr[a], 1.5, f"M(σ{a}) log log x + c(σ{a})  (theory)"))
    figure([dict(curves=curves, ylabel="π½(x) − 4π½(x;σ)", title=fig_titles[0])], 1e3, 1e15, f"{B}/fig_trajectories_{name.split('_')[0]}.svg")
    fig, ax = plt.subplots(figsize=(7.4, 5.2))
    for a in classes:
        ax.semilogx(x[m], LHS[a][m], color=colr[a], lw=0.75, label=f"σ{a}: π½(x) − 4π½(x;σ{a})")
        ax.semilogx(x[m], (Mfun[a]*llx + Csig[a])[m], color=colr[a], lw=1.5, ls="--", label=f"M(σ{a}) log log x + c(σ{a})  (theory)")
    ax.set_xlabel("x"); ax.set_ylabel("π½(x) − 4π½(x;σ)"); ax.set_title(fig_titles[0], fontsize=10); ax.legend(fontsize=7.4, loc="center left")
    fig.tight_layout(); fig.savefig(f"{B}/fig_trajectories_{name.split('_')[0]}.png", bbox_inches="tight"); plt.close(fig)
    # fig 2: data constants
    figure([dict(curves=[(xs, DC[a][sel], colr[a], 0.85, f"σ{a}: data constant") for a in classes],
                 hlines=[(Csig[a], colr[a], f"c(σ{a}) = {Csig[a]:.7f}…") for a in classes], ylabel="data constant", title=fig_titles[1])],
           1e3, 1e15, f"{B}/fig_constants_{name.split('_')[0]}.svg")
    png([dict(curves=[(x[m], DC[a][m], colr[a], 0.8, f"σ{a}: data constant") for a in classes],
              hlines=[(Csig[a], colr[a], f"c(σ{a}) = {Csig[a]:.7f}…") for a in classes], ylabel="data constant", title=fig_titles[1])],
        f"{B}/fig_constants_{name.split('_')[0]}.png", (7.4, 5.0))
    hdr = "| x | " + " | ".join(f"LHS σ{a} | data const. σ{a} (→ {Csig[a]:.7f})" for a in classes) + " |\n|---|" + "---|"*(2*len(classes))
    rows = "\n".join(f"| 10^{k} | " + " | ".join(f"{LHS[a][idx[k]]:.10f} | {DC[a][idx[k]]:.10f}" for a in classes) + " |" for k in ks)
    mm = x >= 1e5
    check = ", ".join(f"σ{a}: {np.mean(DC[a][mm]-Csig[a]):+.4f}" for a in classes)
    open(f"{B}/README.md", "w").write(subs(f"""# L = **{field}**, G = Gal(L/**Q**) = {{{Gdesc}}}

{theorem_md}

## Constants

{consts_md}

## Table

Left-hand sides π_(1/2)(x) − 4π_(1/2)(x; σa) and data constants D_σ(x) = left-hand side − M(σa) log log x.

{hdr}
{rows}

(The full data, 128 checkpoints per decade from x = 100 to 10^15, is in `{csvname}`. Since the class sums are of size 10^6 at x = 10^15, about 10 decimals of the left-hand sides are significant; 10 are printed.)

## Figure 1: the trajectories and the law M(σ) log log x + c(σ)

![trajectories]({f"fig_trajectories_{name.split('_')[0]}.svg"})

The thin curves are π_(1/2)(x) − 4π_(1/2)(x; σa) for the four classes, 10^3 ≤ x ≤ 10^15; the thick curves are M(σa) log log x + c(σa) with the theoretical constants above (no fitting). The sum of the four left-hand sides is 4R for every x.

## Figure 2: the data constants

![data constants]({f"fig_constants_{name.split('_')[0]}.svg"})

The curves are the data constants D_σ(x) for the four classes; the dashed lines are the theoretical constants c(σa). Each curve oscillates around its own constant inside an envelope that narrows like 1/log x.{extra_md}
"""))
    return Csig

# ---- Q(zeta_5) ----
R5 = 1/math.sqrt(5)
c54 = -0.1 + float(d.K3_5.iloc[-1]); c52 = 0.5*float(d.H4.iloc[-1]) + complex(float(d.K3_2re.iloc[-1]), float(d.K3_2im.iloc[-1]))
f5 = {2: np.log(L52) - c52, 3: np.log(L52.conjugate()) - c52.conjugate(), 4: math.log(L54) - c54}
M5 = {a: chi54[a]/2 for a in (1,2,3,4)}
def c5(a):
    s = 0
    for j, chi in ((2, chi52), (3, {k: v.conjugate() for k, v in chi52.items()}), (4, chi54)): s += chi[a].conjugate()*f5[j]
    return (M5[a]*g + R5 - M5[a]*(math.log(2)+cQ) - s).real
C5 = {a: c5(a) for a in (1,2,3,4)}
pihalf5 = d.P5_1.values + d.P5_2.values + d.P5_3.values + d.P5_4.values + R5
LHS5 = {a: pihalf5 - 4*d[f"P5_{a}"].values for a in (1,2,3,4)}
consts5 = f"""| quantity | value |
|---|---|
| L(1/2, χ2) = conj L(1/2, χ3) | {L52.real:.12f} + {L52.imag:.12f} i |
| L(1/2, χ4) | {L54:.15f} |
| m(χ2) = m(χ3) = m(χ4) = ord_(s=1/2) L(s, χ) | 0 |
| lowest non-trivial zero of L(s, χ2): s = 1/2 + i γ1 (for χ3 the sign of γ1 is reversed) | γ1 = {z52:.10f} |
| lowest non-trivial zero of L(s, χ4): s = 1/2 ± i γ1 | γ1 = {z54:.10f} |
| c_Q = γ − M_Mertens | {cQ:.13f} |
| c(χ4) = −1/10 + Σ_(k≥3) Σ_p χ4(p^k)/(k p^(k/2)) | {c54:.10f} |
| c(χ2) = ½ Σ_p χ4(p)/p + Σ_(k≥3) Σ_p χ2(p^k)/(k p^(k/2)) = conj c(χ3) | {c52.real:.10f} {c52.imag:+.10f} i |
| c(σ1) | **{C5[1]:.10f}** |
| c(σ2) | **{C5[2]:.10f}** |
| c(σ3) | **{C5[3]:.10f}** |
| c(σ4) | **{C5[4]:.10f}** |

(Σ_σ c(σ) = 4R = 4/√5 = {4*R5:.10f}.)"""
thm5 = f"""σa is the automorphism ζ5 ↦ ζ5^a, and a prime p ≠ 5 has Frob_p = σa exactly when p ≡ a (mod 5); the prime 5 is ramified, so R = Σ_(p|D_L) p^(−1/2) = 1/√5. The non-trivial irreducible characters of G are one-dimensional: ρ2 with ρ2(σ1) = 1, ρ2(σ2) = i, ρ2(σ3) = −i, ρ2(σ4) = −1; ρ3 = conj ρ2; and ρ4 = ρ2² with ρ4(σ1) = ρ4(σ4) = 1, ρ4(σ2) = ρ4(σ3) = −1, which is real. Since only ρ4 is real, M(σ) = ½ ρ4(σ): **M(σ1) = M(σ4) = 1/2, M(σ2) = M(σ3) = −1/2**, and m(σ) = 0 for all σ. If we identify ρ2, ρ3, ρ4 with the Dirichlet characters χ2, χ3, χ4 modulo 5 of `1_partial_euler_products/mod5_Dirichlet_char/`, then the Artin L-functions for ρ2, ρ3, ρ4 coincide with the Dirichlet L-functions for χ2, χ3, χ4.

The theorem reads, for σ ∈ G,

**π_(1/2)(x) − 4 π_(1/2)(x; σ) = M(σ) log log x + c(σ) + o(1),**

with c(σ) = M(σ) γ + R − M(σ)(log 2 + c_Q) − Σ_(χ = χ2, χ3, χ4) χ̄(σ)( log L(1/2, χ) − c(χ) ), c_Q = − Σ_p (1/p + log(1 − 1/p)), c(χ4) = −1/10 + Σ_(k≥3) Σ_p χ4(p^k)/(k p^(k/2)) and, for the complex characters, c(χ2) = ½ Σ_p χ4(p)/p + Σ_(k≥3) Σ_p χ2(p^k)/(k p^(k/2)), c(χ3) = conj c(χ2). The coefficient of log L(1/2, χ) − c(χ) is the complex conjugate χ̄(σ) = χ(σ)^(−1); for the real character χ4 it is χ4(σ)."""
extra5 = "\n\nThe pair σ2, σ3 (primes ≡ 2 and ≡ 3 mod 5) is not symmetric: the two constants differ by 2 Im( log L(1/2, χ2) − c(χ2) ), which is where the complex characters enter."
bias_page("mod5_Dirichlet_char", "Q(ζ5)", R5, 5, M5, (1,2,3,4), LHS5, C5, consts5, thm5, "σ1, σ2, σ3, σ4",
          ("Q(ζ5): Chebyshev's bias π½(x) − 4π½(x;σ) and the law M(σ) log log x + c(σ)", "Q(ζ5): data constants vs theoretical constants, x ≤ 10^15"), "bias_mod5.csv", extra5)

# ---- Q(zeta_8) ----
R8 = 1/s2
cchi = {j: -0.25 + float(d[{2: "K3_m4", 3: "K3_m8", 4: "K3_8"}[j]].iloc[-1]) for j in (2,3,4)}
M8 = {a: 0.5*sum(ch8[j][a] for j in (2,3,4)) for a in (1,3,5,7)}
C8 = {a: M8[a]*g + R8 - M8[a]*(math.log(2)+cQ) - sum(ch8[j][a]*(math.log(L8[j]) - cchi[j]) for j in (2,3,4)) for a in (1,3,5,7)}
pihalf8 = d.P8_1.values + d.P8_3.values + d.P8_5.values + d.P8_7.values + R8
LHS8 = {a: pihalf8 - 4*d[f"P8_{a}"].values for a in (1,3,5,7)}
consts8 = f"""| quantity | value |
|---|---|
| L(1/2, χ2), L(1/2, χ3), L(1/2, χ4) | {L8[2]:.12f}, {L8[3]:.12f}, {L8[4]:.12f} |
| m(χ2) = m(χ3) = m(χ4) = ord_(s=1/2) L(s, χ) | 0 |
| lowest non-trivial zero s = 1/2 ± i γ1 of L(s, χ2), L(s, χ3), L(s, χ4) | γ1 = {z8[2]:.10f}, {z8[3]:.10f}, {z8[4]:.10f} |
| c_Q = γ − M_Mertens | {cQ:.13f} |
| c(χ2), c(χ3), c(χ4) = −1/4 + Σ_(k≥3) Σ_p χ(p^k)/(k p^(k/2)) | {cchi[2]:.10f}, {cchi[3]:.10f}, {cchi[4]:.10f} |
| c(σ1) | **{C8[1]:.10f}** |
| c(σ3) | **{C8[3]:.10f}** |
| c(σ5) | **{C8[5]:.10f}** |
| c(σ7) | **{C8[7]:.10f}** |

(Σ_σ c(σ) = 4R = 2√2 = {4*R8:.10f}.)"""
thm8 = f"""σa is the automorphism ζ8 ↦ ζ8^a, and an odd prime p has Frob_p = σa exactly when p ≡ a (mod 8); the prime 2 is ramified, so R = Σ_(p|D_L) p^(−1/2) = 1/√2. The non-trivial irreducible characters of G are one-dimensional and real: ρ2, ρ3, ρ4 with (ρ2, ρ3, ρ4)(σ1) = (1, 1, 1), (σ3) = (−1, 1, −1), (σ5) = (1, −1, −1), (σ7) = (−1, −1, 1). Hence M(σ) = ½ (ρ2(σ) + ρ3(σ) + ρ4(σ)): **M(σ1) = 3/2, M(σ3) = M(σ5) = M(σ7) = −1/2**, and m(σ) = 0 for all σ. If we identify ρ2, ρ3, ρ4 with the Dirichlet characters χ2, χ3, χ4 modulo 8 of `1_partial_euler_products/mod8_Dirichlet_char/`, then the Artin L-functions for ρ2, ρ3, ρ4 coincide with the Dirichlet L-functions for χ2, χ3, χ4.

The theorem reads, for σ ∈ G,

**π_(1/2)(x) − 4 π_(1/2)(x; σ) = M(σ) log log x + c(σ) + o(1),**

with c(σ) = M(σ) γ + R − M(σ)(log 2 + c_Q) − Σ_(χ = χ2, χ3, χ4) χ(σ)( log L(1/2, χ) − c(χ) ), c_Q = − Σ_p (1/p + log(1 − 1/p)) and c(χ) = −1/4 + Σ_(k≥3) Σ_p χ(p^k)/(k p^(k/2)) for each of the three characters."""
bias_page("mod8_Dirichlet_char", "Q(ζ8)", R8, 2, M8, (1,3,5,7), LHS8, C8, consts8, thm8, "σ1, σ3, σ5, σ7",
          ("Q(ζ8): Chebyshev's bias π½(x) − 4π½(x;σ) and the law M(σ) log log x + c(σ)", "Q(ζ8): data constants vs theoretical constants, x ≤ 10^15"), "bias_mod8.csv")
print("c(σ) mod5:", {a: round(C5[a], 7) for a in C5}); print("c(σ) mod8:", {a: round(C8[a], 7) for a in C8})
