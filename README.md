# Real-HITRAN CO Spectrum Simulation

Reproduction of the *"CO (v″, v′) = (0,1) rovibrational band, T = 1000 K, P = 1 atm,
χ_CO = 0.01 in air, L = 1 cm"* example spectrum using **real HITRAN line-by-line data**
(downloaded from [hitran.org](https://www.hitran.org) following the Standard "Crash
course on absorption spectroscopy" tutorial) — instead of the two hand-picked lines in
the Appendix D of Hanson et al. 2016.

## Files

| File | Role |
|---|---|
| `hitran_par.py` | HITRAN `.par` **parser** (fixed-width 100/160-char columns → arrays) |
| `voigt_humlicek.py` | **Voigt lineshape** algorithm (Humlíček 1982, 4-region rational approx) + self-checks |
| `simulate_hitran_co.py` | Real-data CO spectrum simulator (temperature scaling → widths → Voigt → sum) |
| `simulate_co_figure.py` | Two-panel figure: full band + zoom (P = 1 vs 5 atm) |
| `co_figure_T1000K_X0.01_L1.svg` | **Result figure** (open in browser) |
| `co_T1000K_X0.01_L1_full.csv` | Grid: ν, α(ν) — full band (1900–2300 cm⁻¹) |
| `co_T1000K_X0.01_L1_zoom_p1.csv` | Grid: ν, α(ν) — zoom (P = 1 atm) |
| `co_T1000K_X0.01_L1_zoom_p5.csv` | Grid: ν, α(ν) — zoom (P = 5 atm) |

> The raw `.par` dataset is stored locally in `hitran_data/` (kept out of the repo to keep
> it clean/download-size small). To run, point the scripts at your own CO `.par`.

## Quick start

```bash
python hitran_par.py <co.par> --band 2000,2300        # parse + sanity checks
python simulate_co_figure.py <co.par>                 # regenerate the two-panel figure
```

Requires Python ≥ 3.8 + NumPy.

## Result (real HITRAN data, 1258 lines in window)

- **Full band** (left panel): peak absorbance **α ≈ 0.20 @ 2206 cm⁻¹** (R-branch), with the
  characteristic P-branch (~2100 cm⁻¹) / R-branch (~2200 cm⁻¹) double-lobe structure.
- **Zoom 2044–2058 cm⁻¹** (right panel): P = 1 atm → α ≈ 0.146, P = 5 atm → α ≈ 0.151
  (@ 2055.4 cm⁻¹); the 5 atm lines are broader and slightly taller (pressure broadening
  ∝ P scales with total absorption).
- **Validation**: integrated absorbance `∫α dν ≈ Σ S·P·X·L` (energy conservation, ratio ≈ 1.000)
  and the anchor line `ν₀ = 2147.08 cm⁻¹, E″ = 0` == CO 1-0 **R(0)**.

## Pipeline implemented

`parse .par → windowing → linestrength temperature scaling (partition/Boltzmann/stimulated-emission)
→ per-atm unit conversion (×7.34e21/T) → Doppler + collisional FWHM (mixture-weighted, HWHM→2γ)
→ Voigt (Humlíček) → α(ν) = Σ_j S_j(T)·P·X·φ_j(ν)·L → transmittance`

Key pitfalls handled: `.par` fixed-width columns (100 vs 160 chars; `n_air` 4-char, `δ` 8-char
fields), HWHM vs FWHM, missing `n_self`/shift temperature exponents (defaults 0.75 / 0.96),
and grid step ≪ line width.
