# -*- coding: utf-8 -*-
"""
official_hapi_check.py — 官方 hapi 取数/解析 与 本地 .par 一致性校验
====================================================================
在几个代表性窗口里，用官方 hapi.fetch_by_ids 抓 CO 谱线，
与本地下载的 .par（同源同版本数据）逐条对比：
  nu0 绝对差、S 相对差、行数 —— 全零/全等则说明官方解析=本地解析。

用法: python official_hapi_check.py <co.par>
"""
import os, sys, ssl

ssl._create_default_https_context = ssl._create_unverified_context
P = r"C:\Users\hong1\workplace\.pylibs"
if P not in sys.path:
    sys.path.insert(0, P)
import numpy as np
import hapi

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hitran_par import parse_par

DB = r"C:\Users\hong1\workplace\.tools\hapi_db"
CO_ISO = [26, 27, 28, 29, 30, 31]   # CO 主同位素体(官方全局 iso id)
WINDOWS = [(1950.0, 1955.0), (2100.0, 2105.0), (2172.0, 2174.0),
           (2255.0, 2262.0)]


def compare_window(local, i, w):
    wlo, whi = w
    os.makedirs(DB, exist_ok=True)
    tab = f"co_w{i}"
    hapi.fetch_by_ids(tab, CO_ISO, wlo, whi)
    nuO = np.asarray(hapi.getColumn(tab, "nu"), float)
    swO = np.asarray(hapi.getColumn(tab, "sw"), float)
    o = np.argsort(nuO)
    nuO, swO = nuO[o], swO[o]

    # 本地 .par 同窗口（iso 已含 1..6，即全局 26..31）
    m = (local["nu0"] >= wlo) & (local["nu0"] <= whi)
    nuL = np.sort(local["nu0"][m])
    # 按 nu 对齐 S
    swL = local["S296"][m][np.argsort(local["nu0"][m])]

    n = min(len(nuO), len(nuL))
    if len(nuO) != len(nuL):
        print(f"[{wlo:.1f}-{whi:.1f}] COUNT differs: official={len(nuO)} "
              f"local={len(nuL)}")
        return False
    dnu = np.abs(nuO[:n] - nuL[:n]).max()
    dsw = np.abs(swO[:n] / swL[:n] - 1.0).max() if n else 0.0
    ok = dnu < 1e-9 and dsw < 1e-6
    print(f"[{wlo:.1f}-{whi:.1f}] lines={n}  max|dnu|={dnu:.2e}  "
          f"max|dS/S|={dsw:.2e}  -> {'OK' if ok else 'MISMATCH'}")
    return ok


def main():
    par = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "hitran_data", "6a9a90ff.par")
    print("local .par:", par)
    local = parse_par(par, mol=5)
    print("local parsed lines:", local["parsed_n"])

    hapi.db_begin(DB)
    results = [compare_window(local, i, w) for i, w in enumerate(WINDOWS)]
    hapi.db_commit()
    print("\nALL WINDOWS MATCH:", all(results))


if __name__ == "__main__":
    main()
