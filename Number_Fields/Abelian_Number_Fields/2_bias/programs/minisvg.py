"""minisvg.py — compact SVG figures (log-x axis) without matplotlib bloat.
panel: dict(curves=[(x, y, color, width, label)], hlines=[(y, color, label)], ylabel, title, ylim=None)
"""
import math, numpy as np

W, H = 760, 300          # panel size (px); margins below
ML, MR, MT, MB = 62, 16, 26, 36

def _fmt(v):
    s = f"{v:.1f}".rstrip("0").rstrip(".")
    return s if s not in ("", "-0") else "0"

def _ticks_y(lo, hi):
    span = hi - lo
    step = 10 ** math.floor(math.log10(span / 4))
    for m in (1, 2, 2.5, 5, 10):
        if span / (step * m) <= 6: step *= m; break
    t0 = math.ceil(lo / step) * step
    ticks = []
    while t0 <= hi + 1e-12: ticks.append(round(t0, 12)); t0 += step
    return ticks

def panel(p, xlo, xhi, oy=0, xaxis=True):
    """returns svg fragment for one panel placed at vertical offset oy"""
    ys = np.concatenate([np.asarray(c[1], float) for c in p["curves"]] + [np.array([h[0]]) for h in p.get("hlines", [])])
    lo, hi = (p["ylim"] if p.get("ylim") else (ys.min(), ys.max()))
    pad = 0.04 * (hi - lo); lo -= pad; hi += pad
    lx0, lx1 = math.log10(xlo), math.log10(xhi)
    def X(x): return ML + (math.log10(x) - lx0) / (lx1 - lx0) * (W - ML - MR)
    def Y(y): return oy + MT + (hi - y) / (hi - lo) * (H - MT - MB)
    out = []
    out.append(f'<rect x="{ML}" y="{oy+MT}" width="{W-ML-MR}" height="{H-MT-MB}" fill="white" stroke="#000" stroke-width="0.8"/>')
    # grid + ticks
    for k in range(math.ceil(lx0), math.floor(lx1) + 1):
        xx = X(10.0**k)
        out.append(f'<line x1="{_fmt(xx)}" y1="{oy+MT}" x2="{_fmt(xx)}" y2="{oy+H-MB}" stroke="#ddd" stroke-width="0.5"/>')
        if xaxis and (k - math.ceil(lx0)) % 2 == 0:
            out.append(f'<text x="{_fmt(xx)}" y="{oy+H-MB+15}" font-size="11" text-anchor="middle">10<tspan dy="-5" font-size="8">{k}</tspan></text>')
    for t in _ticks_y(lo, hi):
        yy = Y(t)
        out.append(f'<line x1="{ML}" y1="{_fmt(yy)}" x2="{W-MR}" y2="{_fmt(yy)}" stroke="#ddd" stroke-width="0.5"/>')
        out.append(f'<text x="{ML-5}" y="{_fmt(yy+4)}" font-size="11" text-anchor="end">{t:g}</text>')
    # hlines
    for (yv, col, lab) in p.get("hlines", []):
        out.append(f'<line x1="{ML}" y1="{_fmt(Y(yv))}" x2="{W-MR}" y2="{_fmt(Y(yv))}" stroke="{col}" stroke-width="1.3" stroke-dasharray="6,4"/>')
    # curves
    for (xs, ys_, col, lw, lab) in p["curves"]:
        xs = np.asarray(xs, float); ys_ = np.asarray(ys_, float)
        m = (xs >= xlo) & (xs <= xhi)
        pts = " ".join(f"{_fmt(X(a))},{_fmt(Y(b))}" for a, b in zip(xs[m], ys_[m]))
        out.append(f'<polyline points="{pts}" fill="none" stroke="{col}" stroke-width="{lw}"/>')
    # legend (top right)
    ly = oy + MT + 14
    for (xs, ys_, col, lw, lab) in p["curves"]:
        out.append(f'<line x1="{W-MR-330}" y1="{ly-4}" x2="{W-MR-300}" y2="{ly-4}" stroke="{col}" stroke-width="{lw}"/>')
        out.append(f'<text x="{W-MR-294}" y="{ly}" font-size="11">{lab}</text>'); ly += 15
    for (yv, col, lab) in p.get("hlines", []):
        out.append(f'<line x1="{W-MR-330}" y1="{ly-4}" x2="{W-MR-300}" y2="{ly-4}" stroke="{col}" stroke-width="1.3" stroke-dasharray="6,4"/>')
        out.append(f'<text x="{W-MR-294}" y="{ly}" font-size="11">{lab}</text>'); ly += 15
    # labels
    if p.get("ylabel"):
        cy = oy + MT + (H - MT - MB) / 2
        out.append(f'<text x="14" y="{_fmt(cy)}" font-size="11" text-anchor="middle" transform="rotate(-90 14 {_fmt(cy)})">{p["ylabel"]}</text>')
    if p.get("title"):
        out.append(f'<text x="{ML + (W-ML-MR)/2}" y="{oy+16}" font-size="13" text-anchor="middle">{p["title"]}</text>')
    if xaxis:
        out.append(f'<text x="{ML + (W-ML-MR)/2}" y="{oy+H-4}" font-size="11" text-anchor="middle">x</text>')
    return "\n".join(out)

def figure(panels, xlo, xhi, path):
    n = len(panels); Ht = H * n
    body = "\n".join(panel(p, xlo, xhi, oy=i*H, xaxis=(i == n-1)) for i, p in enumerate(panels))
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{Ht}" viewBox="0 0 {W} {Ht}" '
           f'font-family="Helvetica, Arial, sans-serif">\n<rect width="{W}" height="{Ht}" fill="white"/>\n{body}\n</svg>\n')
    open(path, "w").write(svg)
    return len(svg)
