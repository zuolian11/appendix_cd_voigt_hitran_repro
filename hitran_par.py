# -*- coding: utf-8 -*-
"""
hitran_par.py — HITRAN .par 逐线数据解析器
===========================================
支持经典 100 字符与扩展 160 字符 .par（hitran.org 网站下载的默认格式），
按官方定宽列切片（前 10 个数值字段 + 量子数区字符串）。

单位/约定（务必注意）：
  * S(296 K) : 每分子线强 [cm^-1/(molecule cm^-2)]，参考温度 296 K
  * gamma_air / gamma_self : HWHM 半高半宽系数 [cm^-1/atm]（FWHM = 2*gamma）
  * Epp       : 低态能量 [cm^-1]
  * n_air     : 空气展宽温度指数（无自展宽指数：默认用 0.75，见速成课）
  * delta     : 空气压力位移 [cm^-1/atm]，可负（无温度指数：默认 0.96）

用法：
    python hitran_par.py <file.par> [--mol 5] [--nu-min A] [--nu-max B] [--iso I]
"""

import sys, argparse
import numpy as np

# 字段 = (名字, 0-based 切片, 类型, 说明)
# 列位依据 HAPI 的 HITRAN_FORMAT_160（1-based pos 减 1 得 0-based）：
#   nu=4(12) S=16(10) A=26(10) gamma_air=36(5) gamma_self=41(5) E''=46(10)
#   n_air=56(4) delta_air=60(8)
NUM_FIELDS = [
    ("nu0",     slice(3, 15),  float, "line-center wavenumber [cm^-1]"),
    ("S296",    slice(15, 25), float, "linestrength @296K [cm^-1/(molecule cm^-2)]"),
    ("A",       slice(25, 35), float, "Einstein A coefficient [s^-1]"),
    ("gam_air", slice(35, 40), float, "air-broadened HWHM @296K [cm^-1/atm]"),
    ("gam_self", slice(40, 45), float, "self-broadened HWHM @296K [cm^-1/atm]"),
    ("Epp",     slice(45, 55), float, "lower-state energy [cm^-1]"),
    ("n_air",   slice(55, 59), float, "air-broadening temperature exponent"),
    ("delta",   slice(59, 67), float, "air pressure shift @296K [cm^-1/atm]"),
]
STR_FIELDS = [
    ("mol",  slice(0, 2),  "molecule number"),
    ("iso",  slice(2, 3),  "isotopologue id (1 = most abundant)"),
    # 扩展 160 字符格式：global upper(67:82) lower(82:97) / local upper(97:112) lower(112:127)
    ("q_glob", slice(67, 97),  "global quanta (v', v'', ...)"),
    ("q_loc",  slice(97, 127), "local quanta (rotational labels)"),
]


def _f(tok):
    """把一段定宽文本转 float；全空或异常返回 nan。"""
    try:
        return float(tok.strip())
    except ValueError:
        return float("nan")


def parse_par(path, mol=None, iso=None, nu_min=None, nu_max=None):
    """
    解析 HITRAN .par 文件。

    参数
    ----
    path   : 文件路径
    mol    : 只保留的分子号（如 CO=5, H2O=1）；None = 全部
    iso    : 只保留的同位素体号；None = 全部
    nu_min, nu_max : 波数窗口过滤 [cm^-1]

    返回
    ----
    dict: {'nu0':ndarray, 'S296':..., 'A':..., 'gam_air':...,
           'gam_self':..., 'Epp':..., 'n_air':..., 'delta':...,
           'mol':ndarray(int), 'iso':ndarray(int),
           'qg':list, 'ql':list,       # 全局/局域量子数原字符串
           'raw_n':int, 'parsed_n':int, 'skipped_n':int}
    """
    rows = {n: [] for n, _, _, _ in NUM_FIELDS}
    rows.update({n: [] for n, _, _ in STR_FIELDS})
    raw_n = skipped = 0

    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\r\n")
            if not line.strip():
                continue
            raw_n += 1
            if len(line) < 100:                 # 至少要有 100 字符核心
                skipped += 1
                continue
            m = line[0:2].strip()
            i = line[2:3].strip()
            if mol is not None and m != str(mol):
                continue
            if iso is not None and i != str(iso):
                continue
            try:
                nu = float(line[3:15])
            except ValueError:
                skipped += 1
                continue
            if nu_min is not None and nu < nu_min:
                continue
            if nu_max is not None and nu > nu_max:
                continue
            for name, sl, _, _ in NUM_FIELDS:
                rows[name].append(_f(line[sl]))
            for name, sl, _ in STR_FIELDS:
                rows[name].append(line[sl].strip())

    parsed_n = len(rows["nu0"])
    out = {name: np.array(rows[name]) for name, _, _, _ in NUM_FIELDS}
    out["mol"] = np.array([int(x) if x.isdigit() else -1 for x in rows["mol"]])
    out["iso"] = np.array([int(x) if x.isdigit() else -1 for x in rows["iso"]])
    out["qg"] = rows["q_glob"]
    out["ql"] = rows["q_loc"]
    out["raw_n"], out["parsed_n"], out["skipped_n"] = raw_n, parsed_n, skipped
    return out


def fmt_row(r, idx, cols=("nu0", "S296", "A", "gam_air", "gam_self",
                          "Epp", "n_air", "delta")):
    """把一行格式化成易读文本，用于人工核对。"""
    hdr = f"{idx:>6}  iso={int(r['iso'][idx])}"
    body = " ".join(f"{c}={r[c][idx]:.6g}" for c in cols)
    return f"{hdr}  {body}"


def summarize(d, top=8, band=(1900.0, 2400.0)):
    """打印解析摘要 + 基频带最强线 + 靠近指定波数的行。"""
    print("parsed lines :", d["parsed_n"], " (skipped:", d["skipped_n"], ")")
    print("nu0 range    : %.4f ... %.4f cm^-1" %
          (d["nu0"].min(), d["nu0"].max()))
    iso_u, iso_c = np.unique(d["iso"], return_counts=True)
    print("isotopologues:", dict(zip(iso_u.tolist(), iso_c.tolist())))

    # 在振动带窗口内按 S 找最强线（S 用 296 K 值）
    if band:
        m = (d["nu0"] >= band[0]) & (d["nu0"] <= band[1])
        if m.any():
            inds = np.argsort(d["S296"][m])[::-1][:top]
            idx_all = np.nonzero(m)[0][inds]
            print(f"\ntop-{top} strongest lines in {band[0]:.0f}-"
                  f"{band[1]:.0f} cm^-1:")
            for ii in idx_all:
                print("   " + fmt_row(d, ii))


def main():
    ap = argparse.ArgumentParser(description="HITRAN .par parser (see docstring)")
    ap.add_argument("path")
    ap.add_argument("--mol", type=int, default=None)
    ap.add_argument("--iso", type=int, default=None)
    ap.add_argument("--nu-min", type=float, default=None, dest="nu_min")
    ap.add_argument("--nu-max", type=float, default=None, dest="nu_max")
    ap.add_argument("--band", type=str, default="1900,2400",
                    help="wavenumber window for 'strongest lines' listing")
    ap.add_argument("--near", type=float, default=None,
                    help="print the rows nearest to this wavenumber")
    args = ap.parse_args()

    d = parse_par(args.path, mol=args.mol, iso=args.iso,
                  nu_min=args.nu_min, nu_max=args.nu_max)
    b = tuple(float(x) for x in args.band.split(","))
    print(f"file   : {args.path}")
    summarize(d, band=b)

    if args.near is not None:
        j = int(np.argmin(np.abs(d["nu0"] - args.near)))
        lo = max(0, j - 2)
        hi = min(len(d["nu0"]), j + 3)
        print(f"\nrows nearest {args.near} cm^-1:")
        for i in range(lo, hi):
            print("   " + fmt_row(d, i))

    if args.mol is None and len(np.unique(d["mol"])) > 1:
        print("\n[note] file contains multiple molecules; use --mol to filter.")


if __name__ == "__main__":
    main()
