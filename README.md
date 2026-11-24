# PaveLab FWD — Backcalculation Studio

Multilayer elastic FWD backcalculation web app (Burmister/Schiffman forward model +
bounded least-squares multi-start inversion, numpy/scipy, Vercel serverless).

Reference methodology: FHWA-HRT-16-011, Chapter 3.

## Validated scope

Validated on exactly four WSDOT SR-395 deflection basins from FHWA-HRT-16-011
Chapter 3 (MP207.85, MP208.00, MP208.50, MP209.00), compared against the EVERCALC
results published in that chapter only. No comparison against BAKFAA or any other
backcalculation program has been performed, and parity beyond these four basins is
not claimed.

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

## Engineering-review flags

A free layer whose backcalculated modulus lands within 2% of its `E_min`/`E_max`
bound is flagged in the API (`at_bound` per layer, `review_required` overall) and in
the UI as needing engineering review: the optimum wanted to pass the bound, so treat
the value as a limit, not a measurement. The response echoes each layer's bounds,
seed, fixed state, thickness and Poisson ratio next to the result.

## Data handling

The browser sends the entered basin and layer data to the server endpoint
`POST /api/backcalc` for computation; the app stores no inputs.

## Stack
- `api/engine.py` — forward model + inversion (pure numpy/scipy)
- `api/backcalc.py` — Vercel serverless handler (`pyproject.toml` `[tool.vercel] entrypoint`)
- `public/` — static enterprise UI (canvas basin chart, measured vs calculated)
