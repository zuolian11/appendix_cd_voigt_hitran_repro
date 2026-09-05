# -*- coding: utf-8 -*-
"""
plot_split_images.py — 把 CO 结果画成两张独立图片
  co_fullband_1atm_T1000K.svg     全带 1900-2300 cm^-1, P=1 atm
  co_zoom_p1_vs_p5atm_T1000K.svg  放大 2044-2058 cm^-1, P=1 vs 5 atm
条件与 simulate_co_figure.py 相同: T=1000 K, X_CO=1%(air), L=1 cm
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hitran_par import parse_par
from simulate_co_figure import spectrum


def single_svg(path, title, xlabel, ylabel, xarr, series, xlim, ylim,
               legend=None):
    """series: list[(y, color, width)]；xarr 与各 y 等长。"""
    W, H, ML, MR, MT, MB = 920, 420, 76, 20, 48, 58
    PW, PH = W - ML - MR, H - MT - MB

    def mxy(x, y):
        return (ML + (x - xlim[0]) / (xlim[1] - xlim[0]) * PW,
                MT + PH - (y - ylim[0]) / (ylim[1] - ylim[0]) * PH)

    out = [f'<line x1="{ML}" y1="{MT}" x2="{ML}" y2="{MT+PH}" stroke="#999"/>',
           f'<line x1="{ML}" y1="{MT+PH}" x2="{W-MR}" y2="{MT+PH}" stroke="#999"/>']
    for t in np.linspace(ylim[0], ylim[1], 6):
        _, py = mxy(xlim[0], t)
        out.append(f'<line x1="{ML-4}" y1="{py:.1f}" x2="{ML}" y2="{py:.1f}" '
                   f'stroke="#666"/><text x="{ML-7}" y="{py+4:.1f}" '
                   f'text-anchor="end" font-size="11">{t:.3f}</text>')
    for x in np.linspace(xlim[0], xlim[1], 7):
        px, _ = mxy(x, ylim[0])
        out.append(f'<line x1="{px:.1f}" y1="{MT+PH}" x2="{px:.1f}" '
                   f'y2="{MT+PH+4}" stroke="#666"/><text x="{px:.1f}" '
                   f'y="{MT+PH+17}" text-anchor="middle" font-size="11">'
                   f'{x:.0f}</text>')
    for ys, color, width in series:
        pts = []
        for xv, yv in zip(xarr, ys):
            px, py = mxy(xv, yv)
            pts.append(f"{px:.1f},{py:.1f}")
        out.append(f'<polyline points="{" ".join(pts)}" fill="none" '
                   f'stroke="{color}" stroke-width="{width}"/>')
    if legend:
        lx = ML + 8
        for label, (color) in legend:
            out.append(f'<line x1="{lx}" y1="{MT+16}" x2="{lx+22}" '
                       f'y2="{MT+16}" stroke="{color}" stroke-width="2.2"/>'
                       f'<text x="{lx+28}" y="{MT+20}" font-size="12" '
                       f'font-family="sans-serif">{label}</text>')
            lx += 34 + 8 * len(label)
    out.insert(0, (f'<text x="{ML}" y="20" font-size="14" '
                   f'font-family="sans-serif">{title}</text>'))
    out.append(f'<text x="{ML+PW/2}" y="{MT+PH+32}" text-anchor="middle" '
               f'font-size="13" font-family="sans-serif">{xlabel}</text>')
    out.append(f'<text x="16" y="{MT+PH/2}" text-anchor="middle" font-size="13" '
               f'font-family="sans-serif" transform="rotate(-90 16 {MT+PH/2})">'
               f'{ylabel}</text>')
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" '
           f'height="{H}"><rect width="{W}" height="{H}" fill="white"/>'
           + "".join(out) + "</svg>")
    with open(path, "w", encoding="utf-8") as f:
        f.write(svg)
    print("->", path)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    par = os.path.join(here, "hitran_data", "6a9a90ff.par")
    T, X, L = 1000.0, 0.01, 1.0
    d = parse_par(par, mol=5, nu_min=1900.0, nu_max=2300.0)
    print("lines in band:", d["parsed_n"])

    # 全带 P=1
    nuL = np.arange(1900.0, 2300.0 + 1e-9, 0.005)
    aL = spectrum(d, nuL, T, 1.0, X, L)
    single_svg(os.path.join(here, "co_fullband_1atm_T1000K.svg"),
               f"CO (0,1) rovibrational band - full band, P=1 atm, T={T:.0f} K, "
               f"X_CO={X:g} in air, L={L:g} cm (real HITRAN, {d['parsed_n']} lines)",
               "Wavenumber [cm-1]", "Absorbance [-ln(I/I0)]",
               nuL, [(aL, "#1f4e79", 1.6)], (1900.0, 2300.0),
               (0.0, float(aL.max()) * 1.15))

    # 缩放 P=1 vs 5
    nuR = np.arange(2044.0, 2058.0 + 1e-9, 0.001)
    aR1 = spectrum(d, nuR, T, 1.0, X, L)
    aR5 = spectrum(d, nuR, T, 5.0, X, L)
    ymax = float(max(aR1.max(), aR5.max())) * 1.15
    single_svg(os.path.join(here, "co_zoom_p1_vs_p5atm_T1000K.svg"),
               f"CO zoom 2044-2058 cm-1 - P=1 vs P=5 atm (T={T:.0f} K, "
               f"X_CO={X:g} in air, L={L:g} cm)",
               "Wavenumber [cm-1]", "Absorbance [-ln(I/I0)]",
               nuR, [(aR1, "#000000", 1.4), (aR5, "#0000bb", 1.4)],
               (2044.0, 2058.0), (0.0, ymax),
               legend=[("P = 1 atm", "#000000"), ("P = 5 atm", "#0000bb")])
    print("peaks:", "full", round(float(aL.max()), 4),
          "| zoom 1atm", round(float(aR1.max()), 4),
          "| zoom 5atm", round(float(aR5.max()), 4))


if __name__ == "__main__":
    main()
