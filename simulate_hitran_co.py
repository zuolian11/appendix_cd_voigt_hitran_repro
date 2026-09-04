# -*- coding: utf-8 -*-
"""
simulate_hitran_co.py — 用真实 HITRAN 数据做 CO 基频带吸收谱模拟
=================================================================
流程（对应速成课第 24 页"Spectrum simulation"）：
  解析 .par -> 窗口内逐线 -> 线强温度标定 -> Doppler/碰撞宽 -> Voigt
  -> alpha_j = S(T) P X phi_j(ν) L 逐线叠加 -> 输出谱 + CSV/SVG

条件（可改）：T=1000 K, P=1 atm, X_CO=1% (空气平衡), L=10 cm
窗口默认取 CO 基频带 2050-2250 cm^-1（可 --window 修改）。
"""
import os, sys, argparse
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hitran_par import parse_par
from voigt_humlicek import voigt_humlicek

HC_K = 1.438776877            # hc/k [cm*K]
T0 = 296.0
CO_OMEGA = 2169.8136          # CO 振动基频 [cm^-1]（RRHO 配分比用，v=1<->0 近似）
M_CO = 28.010                 # g/mol
N_SELF = 0.75                 # 自展宽温度指数（库不提供，速成课默认）
M_AIR  = 0.96                 # 空气压力位移温度指数（速成课默认）


def q_ratio_rrho(T1, T2, omega=CO_OMEGA):
    """RRHO: 线性分子 Q_rot ∝ T，振动单模谐振子 => Q(T1)/Q(T2)。"""
    rot = T1 / T2
    f1 = 1.0 - np.exp(-HC_K * omega / T1)
    f2 = 1.0 - np.exp(-HC_K * omega / T2)
    return rot * (f2 / f1)


def linestrength_scale(S296, Epp, nu0, T, T0=T0):
    """把 S(296K) 按式 D.6（配分函数比用 RRHO 近似）标定到 T。"""
    qr = q_ratio_rrho(T0, T)
    bol = np.exp(-HC_K * Epp * (1.0 / T - 1.0 / T0))
    stim = (1.0 - np.exp(-HC_K * nu0 / T)) / (1.0 - np.exp(-HC_K * nu0 / T0))
    return S296 * qr * bol * stim


def per_atm(s_mol, T):
    """S̃[cm^-1/(molecule cm^-2)] -> S[cm^-2 atm^-1] = S̃ * 7.34e21 / T。"""
    return s_mol * 7.34e21 / T


def doppler_fwhm(nu0, T, M):
    return 7.1623e-7 * nu0 * np.sqrt(T / M)


def main():
    ap = argparse.ArgumentParser(description="Real-HITRAN CO spectrum simulation")
    ap.add_argument("par", help="path to HITRAN .par file (CO)")
    ap.add_argument("--nu-min", type=float, default=2050.0, dest="nu_min")
    ap.add_argument("--nu-max", type=float, default=2250.0, dest="nu_max")
    ap.add_argument("--T", type=float, default=1000.0, dest="T")
    ap.add_argument("--P", type=float, default=1.0)
    ap.add_argument("--X", type=float, default=0.01, dest="X",
                    help="CO mole fraction (rest is air)")
    ap.add_argument("--L", type=float, default=10.0)
    ap.add_argument("--step", type=float, default=0.005,
                    help="wavenumber grid step [cm^-1]")
    ap.add_argument("--show-lines", action="store_true",
                    help="print per-line contributions summary")
    args = ap.parse_args()

    T, P, X, L = args.T, args.P, args.X, args.L
    X_air = 1.0 - X

    d = parse_par(args.par, mol=5, nu_min=args.nu_min, nu_max=args.nu_max)
    n = d["parsed_n"]
    print(f"lines in window {args.nu_min:.0f}-{args.nu_max:.0f}: {n}")
    if n == 0:
        return

    # ---- 逐线参数（矢量化） ----
    nu0 = d["nu0"] + P * X_air * d["delta"] * (T0 / T) ** M_AIR   # 压力位移
    S_atm = per_atm(linestrength_scale(d["S296"], d["Epp"], d["nu0"], T), T)

    gam_air_T = d["gam_air"] * (T0 / T) ** d["n_air"]
    gam_self_T = d["gam_self"] * (T0 / T) ** N_SELF
    dnuC = P * (X_air * 2.0 * gam_air_T + X * 2.0 * gam_self_T)     # FWHM
    dnuD = doppler_fwhm(d["nu0"], T, M_CO)
    a = np.sqrt(np.log(2.0)) * dnuC / dnuD

    # ---- 频率网格 ----
    nu = np.arange(args.nu_min, args.nu_max + args.step * 0.5, args.step)
    alpha = np.zeros_like(nu)

    # 逐线叠加（每条线向量化在整个网格上）
    for j in range(n):
        if dnuD[j] <= 0 or dnuC[j] < 0:
            continue
        alphaD = dnuD[j] / (2.0 * np.sqrt(np.log(2.0)))
        Xg = (nu - nu0[j]) / alphaD
        phi = voigt_humlicek(Xg, a[j]) / (alphaD * np.sqrt(np.pi))
        alpha += S_atm[j] * P * X * phi * L

    trans = np.exp(-alpha)

    # ---- 摘要 ----
    jmax = int(np.argmax(alpha))
    print(f"peak alpha = {alpha[jmax]:.4f} at {nu[jmax]:.3f} cm^-1")
    print(f"min transmittance = {trans.min():.4f}")
    A_int = np.trapezoid(alpha, nu)
    S_tot = float(np.sum(S_atm)) * P * X * L
    print(f"integrated alpha over window = {A_int:.4e}")
    print(f"sum_j S_j P X L (all lines)  = {S_tot:.4e}")
    print(f"(ratio grid/integral ~ window truncation: {A_int / max(S_tot,1e-300):.3f})")

    if args.show_lines:
        order = np.argsort(d["S296"])[::-1][:10]
        print("\ntop-10 strongest lines in window (296 K):")
        print(f"{'nu0':>12} {'S(1000K)/atm':>13} {'Epp':>10} {'gam_air':>8} "
              f"{'n_air':>6} {'iso':>3}")
        for j in order:
            print(f"{d['nu0'][j]:12.4f} {S_atm[j]:13.3e} {d['Epp'][j]:10.3f} "
                  f"{d['gam_air'][j]:8.4f} {d['n_air'][j]:6.2f} {d['iso'][j]:3d}")

    # ---- 输出文件 ----
    here = os.path.dirname(os.path.abspath(__file__))
    base = f"co_T{T:.0f}K_P{P:g}atm_X{X:g}_L{L:g}"
    csv_path = os.path.join(here, f"co_fundamental_spectrum_{base}.csv")
    svg_path = os.path.join(here, f"co_fundamental_spectrum_{base}.svg")
    with open(csv_path, "w", encoding="utf-8") as f:
        f.write("nu_cm-1,alpha,transmittance\n")
        for i in range(len(nu)):
            f.write(f"{nu[i]:.3f},{alpha[i]:.6e},{trans[i]:.6e}\n")
    make_svg(svg_path, nu, alpha,
             f"CO fundamental band (HITRAN real data, {n} lines)  "
             f"T={T:.0f}K P={P:g}atm X={X:g} L={L:g}cm")
    print(f"\nCSV -> {csv_path}")
    print(f"SVG -> {svg_path}")


def make_svg(path, nu, alpha, title):
    W, H, ML, MR, MT, MB = 980, 420, 70, 20, 44, 60
    PW, PH = W - ML - MR, H - MT - MB
    ymax = float(alpha.max()) * 1.1

    def polyl(yvals, color, width=1.6):
        xs = ML + (nu - nu.min()) / (nu.max() - nu.min()) * PW
        ys = MT + PH - (np.asarray(yvals) / ymax) * PH
        pts = " ".join(f"{x:.1f},{y:.1f}" for x, y in zip(xs, ys))
        return (f'<polyline points="{pts}" fill="none" stroke="{color}" '
                f'stroke-width="{width}"/>')

    # 简单 Y 刻度
    ticks = "".join(
        f'<line x1="{ML-4}" y1="{MT+PH-(t/ymax)*PH:.1f}" x2="{ML}" '
        f'y2="{MT+PH-(t/ymax)*PH:.1f}" stroke="#666"/>'
        f'<text x="{ML-8}" y="{MT+PH-(t/ymax)*PH+4:.1f}" text-anchor="end" '
        f'font-size="11">{t:.1f}</text>'
        for t in np.linspace(0, ymax, 6))
    xticks = "".join(
        f'<line x1="{ML+(x-nu.min())/(nu.max()-nu.min())*PW:.1f}" y1="{MT+PH}" '
        f'x2="{ML+(x-nu.min())/(nu.max()-nu.min())*PW:.1f}" y2="{MT+PH+4}" '
        f'stroke="#666"/><text x="{ML+(x-nu.min())/(nu.max()-nu.min())*PW:.1f}" '
        f'y="{MT+PH+18}" text-anchor="middle" font-size="11">{x:.0f}</text>'
        for x in np.linspace(nu.min(), nu.max(), 7))

    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}">
<rect x="0" y="0" width="{W}" height="{H}" fill="white"/>
<text x="{ML}" y="18" font-size="14" font-family="sans-serif">{title}</text>
<line x1="{ML}" y1="{MT}" x2="{ML}" y2="{MT+PH}" stroke="#999"/>
<line x1="{ML}" y1="{MT+PH}" x2="{W-MR}" y2="{MT+PH}" stroke="#999"/>
{ticks}{xticks}
{polyl(alpha, "#1f4e79", 2.2)}
<text x="{ML}" y="{MT+PH+34}" font-size="12" font-family="sans-serif">Wavenumber, cm-1</text>
<text x="16" y="{MT+PH/2}" font-size="12" font-family="sans-serif" transform="rotate(-90 16 {MT+PH/2})" text-anchor="middle">alpha(nu)</text>
</svg>"""
    with open(path, "w", encoding="utf-8") as f:
        f.write(svg)


if __name__ == "__main__":
    main()
