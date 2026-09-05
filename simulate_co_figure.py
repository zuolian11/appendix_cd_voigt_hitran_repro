# -*- coding: utf-8 -*-
"""
simulate_co_figure.py — 复现速成课"Example spectrum"图（CO 基频带）
左图: 全谱带 1900-2300 cm^-1,  T=1000K, P=1atm, X_CO=0.01(air), L=1 cm
右图: 缩放 2044-2058 cm^-1,   P=1atm(黑) vs P=5atm(蓝)
数据: 真实 HITRAN CO .par（用户按 PPT 教程下载的全库文件）
"""
import os, sys, argparse
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hitran_par import parse_par
from voigt_humlicek import voigt_humlicek
from simulate_hitran_co import (linestrength_scale, per_atm,
                                doppler_fwhm, M_CO, N_SELF, M_AIR, T0, HC_K)

CLASSIC = "#000000"
BLUE = "#0000aa"


def spectrum(d, nu, T, P, X, L):
    """用已解析的线表 d 计算频率网格 nu 上的吸光度 alpha。"""
    X_air = 1.0 - X
    nu0 = d["nu0"] + P * X_air * d["delta"] * (T0 / T) ** M_AIR
    S_atm = per_atm(linestrength_scale(d["S296"], d["Epp"], d["nu0"], T), T)
    gam_air = d["gam_air"] * (T0 / T) ** d["n_air"]
    gam_self = d["gam_self"] * (T0 / T) ** N_SELF
    dnuC = P * (X_air * 2.0 * gam_air + X * 2.0 * gam_self)
    dnuD = doppler_fwhm(d["nu0"], T, M_CO)
    a = np.sqrt(np.log(2.0)) * dnuC / np.maximum(dnuD, 1e-12)
    alpha = np.zeros_like(nu)
    # 翼部截断（Voigt 在 |X|>~35 时已可忽略），每条线只写它真正影响到的网格点
    CUT = 35.0
    for j in range(len(d["nu0"])):
        alphaD = dnuD[j] / (2.0 * np.sqrt(np.log(2.0)))
        if alphaD <= 0:
            continue
        half = CUT * alphaD
        lo, hi = np.searchsorted(nu, [nu0[j] - half, nu0[j] + half])
        if hi <= lo:
            continue
        Xg = (nu[lo:hi] - nu0[j]) / alphaD
        phi = voigt_humlicek(Xg, a[j]) / (alphaD * np.sqrt(np.pi))
        alpha[lo:hi] += S_atm[j] * P * X * phi * L
    return alpha


def panel(ax_x, ax_y, w, h, title, xlabel, ylabel, series, xarr, xlim, ylim, xticks):
    """返回该面板的 SVG 字符串（手写，无 matplotlib）。
    series: list of (y_array, color, linewidth)"""
    def mapxy(x, y):
        px = ax_x + (x - xlim[0]) / (xlim[1] - xlim[0]) * w
        py = ax_y + h - (y - ylim[0]) / (ylim[1] - ylim[0]) * h
        return px, py

    out = [f'<line x1="{ax_x}" y1="{ax_y}" x2="{ax_x}" y2="{ax_y+h}" stroke="#999"/>',
           f'<line x1="{ax_x}" y1="{ax_y+h}" x2="{ax_x+w}" y2="{ax_y+h}" stroke="#999"/>']
    # y ticks
    for t in np.linspace(ylim[0], ylim[1], 6):
        _, py = mapxy(xlim[0], t)
        out.append(f'<line x1="{ax_x-4}" y1="{py:.1f}" x2="{ax_x}" y2="{py:.1f}" stroke="#666"/>'
                   f'<text x="{ax_x-6}" y="{py+4:.1f}" text-anchor="end" font-size="11">{t:.2f}</text>')
    for x in xticks:
        px, _ = mapxy(x, ylim[0])
        out.append(f'<line x1="{px:.1f}" y1="{ax_y+h}" x2="{px:.1f}" y2="{ax_y+h+4}" stroke="#666"/>'
                   f'<text x="{px:.1f}" y="{ax_y+h+16}" text-anchor="middle" font-size="11">{x:.0f}</text>')
    for (ys, color, width) in series:
        pts = []
        for xv, yv in zip(xarr, ys):
            xp, yp = mapxy(xv, yv)
            pts.append(f"{xp:.1f},{yp:.1f}")
        out.append(f'<polyline points="{" ".join(pts)}" fill="none" stroke="{color}" '
                   f'stroke-width="{width}"/>')
    out.append(f'<text x="{ax_x}" y="{ax_y-8}" font-size="13" '
               f'font-family="sans-serif">{title}</text>')
    out.append(f'<text x="{ax_x+w/2}" y="{ax_y+h+26}" text-anchor="middle" '
               f'font-size="12" font-family="sans-serif">{xlabel}</text>')
    out.append(f'<text x="{ax_x-40}" y="{ax_y+h/2}" text-anchor="middle" font-size="12" '
               f'font-family="sans-serif" transform="rotate(-90 {ax_x-40} {ax_y+h/2})">{ylabel}</text>')
    return "".join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("par")
    ap.add_argument("--T", type=float, default=1000.0, dest="T")
    ap.add_argument("--X", type=float, default=0.01)
    ap.add_argument("--L", type=float, default=1.0)
    ap.add_argument("--nu-min", type=float, default=1900.0, dest="nu_min")
    ap.add_argument("--nu-max", type=float, default=2300.0, dest="nu_max")
    ap.add_argument("--zoom-a", type=float, default=2044.0, dest="za")
    ap.add_argument("--zoom-b", type=float, default=2058.0, dest="zb")
    ap.add_argument("--step-band", type=float, default=0.005, dest="sband",
                    help="full-band grid step [cm-1] (dense, ~10+ pts per line)")
    ap.add_argument("--step-zoom", type=float, default=0.001, dest="szoom",
                    help="zoom grid step [cm-1] (very dense)")
    args = ap.parse_args()

    T, X, L = args.T, args.X, args.L
    d = parse_par(args.par, mol=5, nu_min=args.nu_min, nu_max=args.nu_max)
    print(f"lines in {args.nu_min:.0f}-{args.nu_max:.0f}: {d['parsed_n']}")

    # ---- 左图: 全谱带 @ P=1 atm ----
    nuL = np.arange(args.nu_min, args.nu_max + 1e-9, args.sband)
    aL = spectrum(d, nuL, T, 1.0, X, L)
    print(f"left: peak alpha={aL.max():.4f} @ {nuL[np.argmax(aL)]:.3f} cm^-1")

    # ---- 右图: 缩放窗口, P=1 vs P=5 ----
    nuR = np.arange(args.za, args.zb + 1e-9, args.szoom)
    aR1 = spectrum(d, nuR, T, 1.0, X, L)
    aR5 = spectrum(d, nuR, T, 5.0, X, L)
    print(f"right: peak(1atm)={aR1.max():.4f} @ {nuR[np.argmax(aR1)]:.3f}; "
          f"peak(5atm)={aR5.max():.4f} @ {nuR[np.argmax(aR5)]:.3f}")

    here = os.path.dirname(os.path.abspath(__file__))
    tag = f"T{T:.0f}K_X{X:g}_L{L:g}"
    # CSV
    for name, nu, a in (("full", nuL, aL), ("zoom_p1", nuR, aR1), ("zoom_p5", nuR, aR5)):
        with open(os.path.join(here, f"co_{tag}_{name}.csv"), "w", encoding="utf-8") as f:
            f.write("nu_cm-1,alpha\n")
            for i in range(len(nu)):
                f.write(f"{nu[i]:.3f},{a[i]:.6e}\n")

    # SVG 双面板
    W, H = 980, 420
    left = panel(60, 50, 420, 300, "Full band (P=1 atm)", "Frequency [cm-1]",
                 "Absorbance [-ln(I/I0)]", [(aL, CLASSIC, 1.0)], nuL,
                 (args.nu_min, args.nu_max), (0, float(aL.max())*1.12),
                 np.linspace(args.nu_min, args.nu_max, 5))
    ymaxR = float(max(aR1.max(), aR5.max())) * 1.12
    right = panel(520, 50, 420, 300, "Zoom (P=1 vs P=5 atm)", "Frequency [cm-1]",
                  "Absorbance [-ln(I/I0)]",
                  [(aR1, CLASSIC, 1.2), (aR5, BLUE, 1.2)], nuR,
                  (args.za, args.zb), (0, ymaxR), np.linspace(args.za, args.zb, 5))

    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}">'
           f'<rect x="0" y="0" width="{W}" height="{H}" fill="white"/>'
           f'<text x="60" y="22" font-size="15" font-family="sans-serif">'
           f'CO (v&#8243;,v&#8242;)=(0,1) rovibrational band, REAL HITRAN data ({d["parsed_n"]} lines): '
           f'T={T:.0f}K, P=1 atm, X_CO={X:g} in air, L={L:g} cm</text>'
           f'<line x1="540" y1="60" x2="560" y2="60" stroke="{CLASSIC}" stroke-width="2"/>'
           f'<text x="565" y="64" font-size="12" font-family="sans-serif">P = 1 atm</text>'
           f'<line x1="640" y1="60" x2="660" y2="60" stroke="{BLUE}" stroke-width="2"/>'
           f'<text x="665" y="64" font-size="12" font-family="sans-serif">P = 5 atm</text>'
           + left + right + '</svg>')
    svg_path = os.path.join(here, f"co_figure_{tag}.svg")
    with open(svg_path, "w", encoding="utf-8") as f:
        f.write(svg)
    print(f"\nSVG -> {svg_path}")
    for n in ("full", "zoom_p1", "zoom_p5"):
        print(f"CSV -> co_{tag}_{n}.csv")


if __name__ == "__main__":
    main()
