# Real-HITRAN CO Spectrum Simulation

Reproduction of the *"CO (v″, v′) = (0,1) rovibrational band, T = 1000 K, P = 1 atm,
χ_CO = 0.01 in air, L = 1 cm"* example spectrum with **real HITRAN line-by-line data**.

## Data handling (route C)

- **Fetching/parsing follows the official HITRANonline library `hapi`
  (`https://hitran.org/hapi/`, `pip install hapi`; official repo hitranonline/hapi)**:
  `official_hapi_check.py` uses `hapi.fetch_by_ids` + `hapi.getColumn` (named columns:
  nu, sw, gamma_air, …) on several windows across the CO band.
- **Validation**: official `hapi` output is compared line-by-line with the local `.par`
  dataset (downloaded from hitran.org following the "Crash course on absorption
  spectroscopy" tutorial) → **identical** (max |Δν| = 0, max |ΔS/S| = 0 over all windows).
- The local `.par` (same official source, full band) is then used for the band-wide
  spectrum simulation below — equivalent data, more convenient for the full range.

## Files

| File | Role |
|---|---|
| `hitran_par.py` | HITRAN `.par` fixed-width parser (official 100/160-char columns, named arrays) |
| `voigt_humlicek.py` | Voigt lineshape (Humlíček 1982) + self-checks |
| `simulate_co_figure.py` | Two-panel figure: full band + zoom (P = 1 vs 5 atm), dense grids |
| `simulate_hitran_co.py` | Single-window spectrum simulator (CLI) |
| `official_hapi_check.py` | Official `hapi` fetch + named-column parsing vs local `.par` consistency check |
| `co_figure_T1000K_X0.01_L1.svg` | **Result figure** (open in browser) |
| `co_T1000K_X0.01_L1_full.csv` | Full band (1900–2300 cm⁻¹), step 0.005 cm⁻¹ (80,000 pts) |
| `co_T1000K_X0.01_L1_zoom_p1.csv` | Zoom 2044–2058 cm⁻¹, step 0.001 cm⁻¹, P = 1 atm |
| `co_T1000K_X0.01_L1_zoom_p5.csv` | Zoom 2044–2058 cm⁻¹, step 0.001 cm⁻¹, P = 5 atm |

## Result (real HITRAN data, 1258 lines in window; dense grids)

| Panel | Result |
|---|---|
| Full band, P = 1 atm (step 0.005 cm⁻¹) | peak **α ≈ 0.219 @ 2193.4 cm⁻¹** (R branch); P/R double-lobe structure |
| Zoom 2044–2058 cm⁻¹ (step 0.001 cm⁻¹) | P = 1 atm → peak 0.146; P = 5 atm → peak 0.150 (@2055.4 cm⁻¹); 5 atm lines broader & slightly taller |

> Grid-resolution note: with the previous coarse step (0.05 cm⁻¹) the full-band peak
> read 0.201 @ 2206 cm⁻¹ — the true resolved peak is **0.219 @ 2193.4 cm⁻¹**. Grid step
> must be ≪ line width (lines here ~0.05 cm⁻¹ FWHM → steps of 0.005/0.001 used).

Validation: integrated absorbance `∫α dν ≈ Σ S·P·X·L` (≈1.000); anchor line
`ν₀ = 2147.08 cm⁻¹, E″ = 0` = CO 1-0 **R(0)**; official-hapi ↔ local-`.par` data agreement.

## Pipeline implemented

`parse .par → windowing → linestrength temperature scaling (partition/Boltzmann/stimulated-emission)
→ per-atm unit conversion (×7.34e21/T) → Doppler + collisional FWHM (mixture-weighted, HWHM→2γ)
→ Voigt (Humlíček) → α(ν) = Σ_j S_j(T)·P·X·φ_j(ν)·L → transmittance`

Key pitfalls handled: fixed-width columns (100 vs 160 chars; `n_air` 4-char, `δ` 8-char),
HWHM vs FWHM, missing `n_self` / shift temperature exponents (defaults 0.75 / 0.96),
grid step ≪ line width.
