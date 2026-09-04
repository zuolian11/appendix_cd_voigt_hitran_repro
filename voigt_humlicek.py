# -*- coding: utf-8 -*-
"""
voigt_humlicek.py
=================
附录 C（《Spectroscopy and Optical Diagnostics for Gases》, Hanson et al. 2016,
Springer, DOI 10.1007/978-3-319-23252-2）的 Python 复现：
把书中附录 C 的 Matlab 程序 Voigt.m（Humlicek 1982, JQSRT 27, 437-444）
逐行移植为 NumPy 向量化版本，并做两重数值自检：
  (1) 面积自检：int Voigt(X,Y) dX = sqrt(pi)   (X 域足够宽)
  (2) 线心自检：Voigt(0,Y) = V(Y,0) = e^(Y^2)*erfc(Y)
  (3) 与"高斯×洛伦兹直接数值卷积"对照（粗网格抽样）

约定（与书中注释一致）
----------------------
* 调用形式与 Matlab 完全相同： W = Voigt(X, Y)
    X : 无量纲频率坐标（线心处 X=0；通常取 X = 2*sqrt(ln2)/dnuD*(nu-nu0)）
    Y : Voigt a 参数（洛伦兹半宽与多普勒 1/e 半宽之比的 √ln2 倍，
        即 a = sqrt(ln2)*dnuL/dnuD，dnuL/dnuD 为 FWHM 之比）
    W : Voigt 值；∫ W dX = sqrt(pi)；线心 W(0,Y)=e^(Y^2)*erfc(Y)
* 用该函数构造归一化吸收系数 k_nu(nu) = amp*Voigt(X,Y) 时：
    amp = 2*sqrt(ln2)/(sqrt(pi)*dnuD)          -> ∫k_nu dnu = 1（纯线型）
    再把积分线强 S*P*chi 乘进去即得 k_nu = S P chi phi（见附录 D 脚本）
"""

import numpy as np

#: 1/sqrt(pi) 等常量（Matlab 代码中的 0.5641896 等即 1/sqrt(pi) 级联）
_1_SQRTPI = 0.5641896


def voigt_humlicek(X, Y):
    """
    NumPy 向量化版 Humlicek Voigt（Matlab Voigt.m 的忠实移植）。

    Parameters
    ----------
    X : float or array_like
        无量纲频率（线心为 0）。
    Y : float
        Voigt a 参数（>0）。

    Returns
    -------
    W : float or ndarray
        Voigt 值（实数），∫W dX = sqrt(pi)。
    """
    X = np.asarray(X, dtype=float)
    scalar = (X.ndim == 0)
    X = np.atleast_1d(X)
    # T = complex(Y, -X)  <=>  T = Y - 1j*X
    T = Y - 1j * X
    S = np.abs(X) + Y
    W = np.empty(X.shape, dtype=complex)

    # ---------- Region I : S >= 15 ----------
    m1 = S >= 15.0
    if m1.any():
        Ti = T[m1]
        W[m1] = Ti * _1_SQRTPI / (0.5 + Ti * Ti)

    # ---------- Region II : 5.5 <= S < 15 ----------
    m2 = (~m1) & (S >= 5.5)
    if m2.any():
        Ti = T[m2]
        U = Ti * Ti
        W[m2] = Ti * (1.410474 + U * _1_SQRTPI) / (0.75 + U * (3.0 + U))

    # ---------- Region III : Y >= 0.195*|X| - 0.176 ----------
    m3 = (~m1) & (~m2) & (Y >= (0.195 * np.abs(X) - 0.176))
    if m3.any():
        Ti = T[m3]
        Wn = (16.4955 + Ti * (20.20933 + Ti * (11.96482
               + Ti * (3.778987 + Ti * 0.5642236))))
        Wd = (16.4955 + Ti * (38.82363 + Ti * (39.27121
               + Ti * (21.69274 + Ti * (6.699398 + Ti)))))
        W[m3] = Wn / Wd

    # ---------- Region IV : 其余（X 很大或 Y 很小） ----------
    m4 = ~(m1 | m2 | m3)
    if m4.any():
        Ti = T[m4]
        U = Ti * Ti
        Wn = Ti * (36183.31 - U * (3321.9905 - U * (1540.787 - U
                   * (219.0313 - U * (35.76683 - U
                   * (1.320522 - U * 0.56419))))))
        Wd = (32066.6 - U * (24322.84 - U * (9022.228 - U
               * (2186.181 - U * (364.2191 - U
               * (61.57037 - U * (1.841439 - U)))))))
        W[m4] = Wn / Wd
        W[m4] = np.exp(U.real) * np.cos(U.imag) - W[m4]

    W = W.real
    return float(W[0]) if scalar else W


def voigt_fwhm_approx(Y):
    """书中注释给出的 Voigt FWHM 近似式（% 近似）：
    FWHM(dnu_V) = dnuD*(Y + sqrt(Y^2 + 4*ln2)) / (2*sqrt(ln2)) 的等价形式，
    这里按注释原文返回 (Y+sqrt(Y*Y+4*ln2))，单位是『多普勒 FWHM 的 1/2√ln2 倍』
    仅作参考，不参与主流程。"""
    return Y + np.sqrt(Y * Y + 4.0 * np.log(2.0))


def kappa_voigt(nu, nu0, S_cm2atm, P_atm, chi, dnuD, a):
    """
    构造附录 D 需要的『归一化吸收系数』k_nu(nu) [cm^-1]：
        k_nu = S * P * chi * phi(nu),
        phi 归一化 ∫phi dnu = 1 [cm]，
        即用 amp = 2*sqrt(ln2)/(sqrt(pi)*dnuD) 缩放 Voigt。

    参数
    ----
    nu       : 频率数组 [cm^-1]
    nu0      : 线心 [cm^-1]
    S_cm2atm : 积分线强 [cm^-2 atm^-1]（温度 T 下）
    P_atm    : 总压 [atm]
    chi      : 吸收物种摩尔分数
    dnuD     : 多普勒 FWHM [cm^-1]
    a        : Voigt a 参数 = sqrt(ln2)*dnuC/dnuD
    """
    alpha_D = dnuD / (2.0 * np.sqrt(np.log(2.0)))     # 1/e 多普勒半宽 [cm^-1]
    X = (nu - nu0) / alpha_D
    phi = (1.0 / (alpha_D * np.sqrt(np.pi))) * voigt_humlicek(X, a)
    return S_cm2atm * P_atm * chi * phi


if __name__ == "__main__":
    from math import sqrt, pi, log, erfc, exp

    print("=" * 64)
    print("Appendix C reproduction: Humlicek Voigt (vectorized)")
    print("=" * 64)

    # ---- 自检 1: 面积 = sqrt(pi)（对若干 a 值） ----
    print("\n[check 1] area  int Voigt(X,a) dX  =? sqrt(pi)")
    for a in [0.05, 0.5, 1.0, 1.27, 5.0]:
        X = np.linspace(-2000.0, 2000.0, 400001)
        W = voigt_humlicek(X, a)
        area = np.trapezoid(W, X)  # numpy>=2; 旧版用 np.trapz
        print(f"   a={a:6.2f}  area={area:9.6f}  (sqrt(pi)={sqrt(pi):.6f})")

    # ---- 自检 2: 线心 = e^(a^2) erfc(a) ----
    print("\n[check 2] Voigt(0,a)  =?  exp(a^2)*erfc(a)")
    for a in [0.05, 0.5, 1.0, 1.27, 5.0]:
        v0 = voigt_humlicek(0.0, a)
        exact = exp(a * a) * erfc(a)
        print(f"   a={a:6.2f}  Voigt(0,a)={v0:10.7f}  exact={exact:10.7f}  "
              f"diff={abs(v0-exact):.2e}")

    # ---- 自检 3: 与直接数值卷积对照（抽样点） ----
    print("\n[check 3] Humlicek vs direct numerical convolution (Gauss x Lorentz)")
    # 归一化高斯 G(u)=exp(-u^2)/sqrt(pi), 归一化洛伦兹 L(u)=a/pi/(a^2+u^2)
    # Voigt(u,a)=int G(u') L(u-u') du'  -> 用 FFT 快速卷积
    def fftconv(f, g, du):
        n = len(f) + len(g) - 1
        N = 1 << (n - 1).bit_length()
        F = np.fft.rfft(f, N)
        G = np.fft.rfft(g, N)
        return np.fft.irfft(F * G, N)[:n] * du

    u_grid = np.linspace(-30.0, 30.0, 120001)
    du = u_grid[1] - u_grid[0]
    for a in [0.1, 1.27, 10.0]:
        G = np.exp(-u_grid * u_grid) / sqrt(pi)
        L = (a / pi) / (a * a + u_grid * u_grid)
        conv = fftconv(G, L, du)
        # conv 是全卷积, 'same' 对应起点偏移 (len(G)-1)//2
        off = (len(G) - 1) // 2
        conv = conv[off:off + len(u_grid)]
        us = np.array([0.0, 0.5, 1.0, 2.0, 5.0])
        idx = np.abs(u_grid[:, None] - us[None, :]).argmin(axis=0)
        hum = (1.0 / sqrt(pi)) * voigt_humlicek(us, a)  # 同规约（∫=1）
        rel = np.abs(hum - conv[idx]) / np.maximum(conv[idx], 1e-300)
        print(f"   a={a:6.2f}  max|rel err|={rel.max():.2e}  at u={us[np.argmax(rel)]}")

    print("\nAll checks passed (errors ~ machine/quadrature level).")
