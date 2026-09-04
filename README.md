# Appendix C & D Reproduction — Voigt profile & HITRAN H₂O spectrum example

Python reproduction of **Appendix C** (Humlíček Voigt fitting program) and **Appendix D.1** (H₂O doublet absorbance example) of:

> R. K. Hanson, R. M. Spearrin, C. S. Goldenstein, *Spectroscopy and Optical Diagnostics for Gases*, Springer International Publishing Switzerland, 2016.
> ISBN 978-3-319-23251-5 · DOI [10.1007/978-3-319-23252-2](https://doi.org/10.1007/978-3-319-23252-2)

## Files

| File | Description |
|---|---|
| `voigt_humlicek.py` | Faithful NumPy port of the book's Appendix C `Voigt.m` (Humlíček 1982) + 3 self-checks (area = √π, line-center `e^(a²)erfc(a)`, cross-check vs FFT convolution) |
| `reproduce_appendixD.py` | Full Appendix D.1 example: HITRAN linestrength temperature scaling → unit conversion → Doppler/collisional widths → Voigt → absorbance spectrum |
| `h2o_doublet_absorbance.csv` | Grid data: α₁, α₂, α_total, transmittance vs wavenumber |
| `h2o_doublet_spectrum.svg` | Absorbance spectrum figure (browser-viewable; generated without matplotlib) |
| `附录C与D_复现与讲解.md` | Detailed walkthrough in Chinese: physics, formulas, knowledge & skill checklist, exercises |

## Quick start

```bash
python voigt_humlicek.py        # self-checks of the Humlíček implementation
python reproduce_appendixD.py   # reproduces the H₂O doublet example (prints + CSV/SVG)
```

Requires only Python ≥ 3.8 with NumPy.

## Example being reproduced (Appendix D.1)

H₂O doublet near 7185.6 cm⁻¹ (≈1392.67 nm), HITRAN2012 parameters (Table D.3):

| Line | ν₀ [cm⁻¹] | S(296 K) [cm⁻¹/(molecule·cm⁻²)] | γ_air | γ_self | E″ [cm⁻¹] | n_air |
|---|---|---|---|---|---|---|
| 1 | 7185.596571 | 2.00e-22 | 0.0342 | 0.371 | 1045.0583 | 0.62 |
| 2 | 7185.596909 | 5.98e-22 | 0.0421 | 0.195 | 1045.0577 | 0.62 |

Conditions: **T = 1000 K, P = 1 atm, 10 % H₂O in air, L = 10 cm**.

### Reproduction vs book values

| Quantity | This repo | Book |
|---|---|---|
| S₁(1000 K), per molecule | 1.025e-21 | 1.02e-21 |
| S₂(1000 K), per molecule | 3.066e-21 | 3.05e-21 |
| Doppler FWHM | 0.0383 cm⁻¹ | 0.0384 cm⁻¹ |
| Collisional FWHM (Line 1 / 2) | 0.0587 / 0.0513 cm⁻¹ | 0.0587 / 0.0513 cm⁻¹ |
| S₁ → cm⁻² atm⁻¹ (1000 K) | 7.53e-3 | 7.487e-3 |
| S₂ → cm⁻² atm⁻¹ (1000 K) | 2.25e-2 | printed 2.237e-3 ⚠️ |
| Total peak absorbance | ≈ 0.29 | Fig. D.3 ≈ 0.3 |

⚠️ The book prints the per-atm linestrength of Line 2 as `2.237e-3`; combined with its own
`S₂(1000 K) = 3.05e-21` and the conversion factor `7.34e21/T` this appears to be an
order-of-magnitude typo — the self-consistent value is ≈ `2.24e-2 cm⁻² atm⁻¹`, which also
agrees with the total peak absorbance ≈ 0.3 in the book's Fig. D.3. See the walkthrough
document for the full argument.

## Pipeline implemented

`look up HITRAN params → scale S(T) → convert units → Doppler & collisional widths →
Voigt (Humlíček) → per-line αⱼ(ν) = Sⱼ(T)·P·X·φⱼ(ν)·L → sum → transmittance`

The linestrength scaling uses the RRHO approximation of Eq. D.2 (within ~2 % of the
HITRAN-partition-function result D.6 for 296–1500 K, per the book's Fig. D.1).

*Educational notes (in Chinese) are included in `附录C与D_复现与讲解.md`.*
