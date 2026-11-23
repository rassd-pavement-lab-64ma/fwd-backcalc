# PaveLab FWD — Backcalculation Studio

Multilayer elastic FWD backcalculation web app (Burmister/Schiffman forward model +
bounded least-squares multi-start inversion, numpy/scipy, Vercel serverless).

Reference methodology: FHWA-HRT-16-011, Chapter 3.

## Validation vs EVERCALC (FHWA-HRT-16-011 Ch.3, WSDOT SR-395 basins)

| Section   | This engine (MPa)      | EVERCALC (MPa)       | RMS % (this / EVERCALC) |
|-----------|------------------------|----------------------|-------------------------|
| MP208.00  | 2,998 / 120 / 117      | 2,882 / 124 / 103    | 0.70 / 0.70             |
| MP207.85  | 3,542 / 118 / 94       | 3,640 / 117 / 83     | 3.34 / 3.09             |
| MP209.00  | 1,945 / 35 / 98        | 1,951 / 34 / 90      | 12.54 / 13.40           |
| MP208.50* | 13,773 / 165 / 282     | 1,703 / 117 / 248    | 2.65 / 3.39             |

*Thin-HMA section (119 mm): backcalculation is ill-conditioned; the engine reaches a
lower-RMS optimum at the HMA upper bound while EVERCALC reports a lower-modulus local
solution. Use per-layer bounds or fix the HMA modulus (lab/temperature-adjusted) for
thin sections. Boussinesq half-space check: < 0.02 % error.

## Stack
- `api/engine.py` — forward model + inversion (pure numpy/scipy)
- `api/backcalc.py` — Vercel serverless handler (`pyproject.toml` `[tool.vercel] entrypoint`)
- `public/` — static enterprise UI (canvas basin chart, measured vs calculated)
