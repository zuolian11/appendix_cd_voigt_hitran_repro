# -*- coding: utf-8 -*-
"""
liuze
"""

import numpy as np
import os
from voigt_humlicek import voigt_humlicek

# ----------------------------------------------------------------------
# 0. 物理常量
# ----------------------------------------------------------------------
HC_K  = 1.438776877  # hc/k [cm*K]
KELVIN_CONV = 7.34e21  # n/P = 1013250/(k[erg/K]) [molecules cm^-3 atm^-1 * K]

# H2O 三个基频（式 D.2 末尾给出，cm^-1）：用于 RRHO 振动配分函数比
H2O_NU_I = (3657.05, 1594.75, 3755.93)

T0 = 296.0          # HITRAN 参考温度 [K]
T  = 1000.0         # 算例温度 [K]
P  = 1.0            # 总压 [atm]
X_H2O = 0.10        # H2O 摩尔分数（空气中）
X_AIR = 1.0 - X_H2O # 其余为空气
L  = 10.0           # 光程 [cm]
M_H2O = 18.015      # 分子量 [g/mol]（书中用 18 也可, 影响 <0.1 %）

# 书内 Table D.3 (HITRAN2012): 每条线的 (nu0, S(T0), gamma_air, gamma_self,
#                               E'', n_air)
#   S(T0) 单位: cm^-1/(molecule cm^-2); gamma 为 HWHM 系数 [cm^-1 atm^-1]@296K
LINES = [
    #  Line  nu0[cm^-1]    S(296K)      gamma_air   gamma_self   E''[cm^-1]  n_air
    dict(nu0=7185.596571, S=2.00e-22, ga=0.0342, gs=0.371, Epp=1045.0583, n=0.62),
    dict(nu0=7185.596909, S=5.98e-22, ga=0.0421, gs=0.195, Epp=1045.0577, n=0.62),
]

N_SELF = 0.75   # 书内假设的自加宽温度指数


# ----------------------------------------------------------------------
# 1. 线强温度标定
# ----------------------------------------------------------------------
def q_ratio_rrho(T1, T2):
    """RRHO 近似下 Q(T1)/Q(T2)（H2O: 非对称陀螺 Qrot ~ T^{3/2} × Π 谐振子）
    只差一个与 T 无关的常数, 比值可直接写为："""
    rot = (T1 / T2) ** 1.5
    vib = 1.0
    for nu_i in H2O_NU_I:
        f1 = 1.0 - np.exp(-HC_K * nu_i / T1)
        f2 = 1.0 - np.exp(-HC_K * nu_i / T2)
        vib *= f2 / f1
    return rot * vib


def linestrength_per_molecule(S_ref, Epp, nu0, T, T0=296.0):
    """式 D.2/D.6：把 HITRAN 的 S(T0)[cm^-1/(molecule cm^-2)] 标定到 T。
    S(T) = S(T0) * [Q(T0)/Q(T)] * exp[-(hc/k)E''(1/T - 1/T0)]
           * {1-exp(-hc nu0/kT)}/{1-exp(-hc nu0/kT0)}          (式 D.6)
    其中 Q(T0)/Q(T) 用 RRHO（对 H2O 即式 D.2, 与 HITRAN Q 差 <2 %@296-1500 K）
    """
    qr = q_ratio_rrho(T0, T)              # Q(T0)/Q(T)
    bol = np.exp(-HC_K * Epp * (1.0 / T - 1.0 / T0))
    stim = (1.0 - np.exp(-HC_K * nu0 / T)) / (1.0 - np.exp(-HC_K * nu0 / T0))
    return S_ref * qr * bol * stim


def s_per_atm(s_per_molecule, Temp):
    """式 D.3 单位换算 [cm^-1/(molecule cm^-2)] -> [cm^-2 atm^-1]：
    S[cm^-2/atm] = S[cm^-1/(molecule cm^-2)] * 7.34e21 / T
    （7.34e21/T 即 1 atm 下的分子数密度 n[cm^-3] 除以 P[atm]，见式 7.66-7.68）"""
    return s_per_molecule * KELVIN_CONV / Temp


# ----------------------------------------------------------------------
# 2. 展宽
# ----------------------------------------------------------------------
def doppler_fwhm(nu0, T, M):
    """式 8.25: dnuD = 7.1623e-7 * nu0 * sqrt(T/M)  [cm^-1]"""
    return 7.1623e-7 * nu0 * np.sqrt(T / M)


def gamma_T(gamma_ref, T, n, T0=296.0):
    """式 8.21: gamma(T) = gamma(T0) * (T0/T)^n（HWHM 系数）"""
    return gamma_ref * (T0 / T) ** n


def collisional_fwhm(ga, gs, n, T, P, X_h2o, X_air, T0=296.0):
    """式 8.19: dnuC = P * sum_k X_k * 2*gamma_k(T)
    自加宽用温度指数 0.75, 空气加宽用库给 n_air。返回 FWHM [cm^-1]"""
    g_air  = gamma_T(ga, T, n,    T0)
    g_self = gamma_T(gs, T, N_SELF, T0)
    return P * (X_air * 2.0 * g_air + X_h2o * 2.0 * g_self)


# ----------------------------------------------------------------------
# 3. 组装
# ----------------------------------------------------------------------
def main():
    print("=" * 78)
    print("Appendix D.1 reproduction - H2O doublet @ 1392.67 nm "
          "(T=1000 K, P=1 atm, X=10%, L=10 cm)")
    print("=" * 78)

    rows = []
    for k, ln in enumerate(LINES, start=1):
        S1000_m = linestrength_per_molecule(ln["S"], ln["Epp"], ln["nu0"], T)
        S1000_a = s_per_atm(S1000_m, T)
        dnuD = doppler_fwhm(ln["nu0"], T, M_H2O)
        dnuC = collisional_fwhm(ln["ga"], ln["gs"], ln["n"], T, P,
                                X_H2O, X_AIR)
        a = np.sqrt(np.log(2.0)) * dnuC / dnuD   # Voigt a 参数
        rows.append((k, ln["nu0"], S1000_m, S1000_a, dnuD, dnuC, a))
        print(f"\nLine {k}: nu0 = {ln['nu0']:.6f} cm^-1")

        print(f"  [D.1.1] S(1000K) = {S1000_m:.3e} cm^-1/(molecule cm^-2)"
              f"   (book: {1.02e-21 if k==1 else 3.05e-21:.2e})")
        print(f"          -> S[cm^-2/atm]@1000K = {S1000_a:.3e}"
              f"   (book prints: {'7.487e-3' if k==1 else '2.237e-3'})")
        print(f"  [D.1.2] Doppler  FWHM = {dnuD:.4f} cm^-1"
              f"   (book: 0.0384)")
        print(f"          Collis.  FWHM = {dnuC:.4f} cm^-1"
              f"   (book: {'0.0587' if k==1 else '0.0513'})")
        print(f"          Voigt a param  = {a:.4f}")

    # ---- 吸光度谱 ----
    nu = np.arange(7185.20, 7186.00 + 1e-9, 0.001)
    alpha = np.zeros_like(nu)
    series = []
    for k, ln in enumerate(LINES, start=1):
        S1000_m = linestrength_per_molecule(ln["S"], ln["Epp"], ln["nu0"], T)
        S1000_a = s_per_atm(S1000_m, T)
        dnuD = doppler_fwhm(ln["nu0"], T, M_H2O)
        dnuC = collisional_fwhm(ln["ga"], ln["gs"], ln["n"], T, P, X_H2O, X_AIR)
        a = np.sqrt(np.log(2.0)) * dnuC / dnuD
        alpha_D = dnuD / (2.0 * np.sqrt(np.log(2.0)))   # 1/e 半宽
        X = (nu - ln["nu0"]) / alpha_D
        phi = voigt_humlicek(X, a) / (alpha_D * np.sqrt(np.pi))  # ∫φ dnu=1
        alpha_j = S1000_a * P * X_H2O * phi * L
        series.append(alpha_j)
        alpha += alpha_j

    trans = np.exp(-alpha)

    imax = int(np.argmax(alpha))
    i1 = int(np.argmax(series[0]))
    i2 = int(np.argmax(series[1]))
    print("\n[D.1.3] peak values on grid:")
    print(f"   line-1 peak alpha = {series[0][i1]:.4f}  at {nu[i1]:.4f} cm^-1")
    print(f"   line-2 peak alpha = {series[1][i2]:.4f}  at {nu[i2]:.4f} cm^-1")
    print(f"   total    peak alpha = {alpha[imax]:.4f}  at {nu[imax]:.4f} cm^-1")
    print(f"   min transmittance   = {trans.min():.4f}")

    # 积分吸光度（数值） vs S*P*X*L
    A_int = np.trapezoid(alpha, nu)
    S_tot = sum(s_per_atm(
        linestrength_per_molecule(ln["S"], ln["Epp"], ln["nu0"], T),
        T) for ln in LINES)
    print(f"   integrated alpha (window 7185.2-7186.0) = {A_int:.4e}")
    print(f"   expected S1*P*X*L + S2*P*X*L             = {S_tot*P*X_H2O*L:.4e}")

    # ---- 写 CSV ----
    here = os.path.dirname(os.path.abspath(__file__))
    csv_path = os.path.join(here, "h2o_doublet_absorbance.csv")
    with open(csv_path, "w", encoding="utf-8") as f:
        f.write("nu_cm-1,alpha_line1,alpha_line2,alpha_total,transmittance\n")
        for i in range(len(nu)):
            f.write(f"{nu[i]:.3f},{series[0][i]:.6e},{series[1][i]:.6e},"
                    f"{alpha[i]:.6e},{trans[i]:.6e}\n")
    print(f"\nCSV written: {csv_path}")

    # ---- 写 SVG 图（无 matplotlib 依赖, 浏览器可直接打开） ----
    make_svg(os.path.join(here, "h2o_doublet_spectrum.svg"), nu, series, alpha)


def make_svg(path, nu, series, alpha):
    W, H, ML, MR, MT, MB = 920, 400, 70, 20, 40, 60
    PW, PH = W - ML - MR, H - MT - MB

    def to_xy(nu_arr, y_arr, ymax):
        x = ML + (nu_arr - nu.min()) / (nu.max() - nu.min()) * PW
        y = MT + PH - (y_arr / ymax) * PH
        return x, y

    def polyline(nu_arr, y_arr, ymax, color, width=1.6):
        xs, ys = to_xy(nu_arr, y_arr, ymax)
        pts = " ".join(f"{x:.1f},{y:.1f}" for x, y in zip(xs, ys))
        return (f'<polyline points="{pts}" fill="none" stroke="{color}" '
                f'stroke-width="{width}"/>')

    ymax = float(alpha.max()) * 1.12
    ticks = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30]
    tick_svg = ""
    for t in ticks:
        y = MT + PH - (t / ymax) * PH
        tick_svg += (f'<line x1="{ML-4}" y1="{y:.1f}" x2="{ML}" y2="{y:.1f}" '
                     f'stroke="#666"/><text x="{ML-8}" y="{y+4:.1f}" '
                     f'text-anchor="end" font-size="11">{t:.2f}</text>')
    xlab = 7185.2
    xlab_svg = ""
    while xlab <= 7186.0001:
        x = ML + (xlab - nu.min()) / (nu.max() - nu.min()) * PW
        xlab_svg += (f'<line x1="{x:.1f}" y1="{MT+PH}" x2="{x:.1f}" '
                     f'y2="{MT+PH+4}" stroke="#666"/>'
                     f'<text x="{x:.1f}" y="{MT+PH+18}" text-anchor="middle" '
                     f'font-size="11">{xlab:.2f}</text>')
        xlab += 0.2

    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}">
<rect x="0" y="0" width="{W}" height="{H}" fill="white"/>
<title>H2O doublet absorbance @ 1000 K, 1 atm, 10% H2O, L=10 cm</title>
<text x="{ML}" y="18" font-size="14" font-family="sans-serif">
H2O doublet @ ~7185.6 cm-1 (1392.67 nm): absorbance alpha(nu)  T=1000K P=1atm X=10% L=10cm</text>
<g stroke="#ddd"><line x1="{ML}" y1="{MT}" x2="{ML}" y2="{MT+PH}"/>
<line x1="{ML}" y1="{MT+PH}" x2="{W-MR}" y2="{MT+PH}"/></g>
{tick_svg}{xlab_svg}
{polyline(nu, series[0], ymax, "#4c72b0")}
{polyline(nu, series[1], ymax, "#dd8452")}
{polyline(nu, alpha, ymax, "#000000", 2.4)}
<line x1="{ML+10}" y1="{MT+10}" x2="{ML+34}" y2="{MT+10}" stroke="#4c72b0"/>
<text x="{ML+40}" y="{MT+14}" font-size="12" font-family="sans-serif">Line 1</text>
<line x1="{ML+90}" y1="{MT+10}" x2="{ML+114}" y2="{MT+10}" stroke="#dd8452"/>
<text x="{ML+120}" y="{MT+14}" font-size="12" font-family="sans-serif">Line 2</text>
<line x1="{ML+180}" y1="{MT+10}" x2="{ML+204}" y2="{MT+10}" stroke="#000000" stroke-width="2.4"/>
<text x="{ML+210}" y="{MT+14}" font-size="12" font-family="sans-serif">Sum</text>
<text x="{ML}" y="{MT+PH+32}" font-size="12" font-family="sans-serif">Wavenumber, cm-1</text>
<text x="16" y="{MT+PH/2}" font-size="12" font-family="sans-serif" transform="rotate(-90 16 {MT+PH/2})" text-anchor="middle">alpha(nu)</text>
</svg>"""
    with open(path, "w", encoding="utf-8") as f:
        f.write(svg)
    print(f"SVG written: {path}")


if __name__ == "__main__":
    main()
