# First-pass halo fit — machinery output, NOT a measurement

> 🚨 **SUPERSEDED the same day** by
> `results/validation_and_corrected_fit_2026-09-04.md`. Every template in
> this fit omitted the pixel solid-angle factor (cos b, a 2x latitude-
> dependent shape error across the ROI), found by the validation step
> hours after this was written. The numbers below are kept as the record
> of what the uncorrected pipeline produced; do not use them.

**Date:** 2026-09-04. **Run:** `runs/fff2de5bc7e0` (w009–w791, UltraClean,
zenith < 100°). **Model:** 19 background components (ring-based diffuse, one
f_l per bin each) + one NFW-ρ² halo template at Totani's Via Lactea II
parameters, amplitude free in sign (his §2.3 convention). Full per-bin output
with background amplitudes: `results/fits/halo_fit_2026-09-04.json`.

## The numbers

| bin | E_c (GeV) | halo counts | observed | fraction | TS | mean\|σ\| halo cells (with / without halo) |
|---|---|---|---|---|---|---|
| 0 | 1.51 | 67,355 | 1,860,861 | 3.6% | 392.7 | 6.88 / 6.94 |
| 1 | 2.55 | 52,949 | 967,387 | 5.5% | 485.0 | 5.09 / 5.21 |
| 2 | 4.31 | 34,949 | 462,354 | 7.6% | 446.2 | 3.41 / 3.53 |
| 3 | 7.28 | 18,720 | 209,493 | 8.9% | 313.6 | 2.39 / 2.59 |
| 4 | 12.29 | 9,995 | 95,415 | 10.5% | 228.6 | 1.82 / 1.96 |
| 5 | 20.76 | 5,433 | 43,195 | 12.6% | 172.5 | 1.46 / 1.59 |
| 6 | 35.06 | 2,557 | 20,497 | 12.5% | 86.1 | 1.06 / 1.12 |
| 7 | 59.22 | 1,189 | 8,892 | 13.4% | 42.4 | 1.04 / 1.10 |
| 8 | 100.02 | 423 | 3,645 | 11.6% | 13.1 | 0.89 / 0.89 |

All 9 bins converged; profile-likelihood 1σ intervals are in the JSON (typical
half-width 5–20% of the amplitude).

## What can and cannot be read from this

**The halo amplitude is positive in every bin.** The fit wants a
centrally-brightening, halo-shaped component at the 4–13% level of the
halo-region counts, rising with energy to a ~12–13% plateau.

🚨 **The TS values are formally enormous (up to 485 ≈ "22σ") and MUST NOT be
quoted as a significance.** The background-only fit already leaves mean
per-cell residuals of ~7σ at bin 0 — the background model is inadequate at a
level far beyond Poisson, so ANY additional smooth template will absorb some
residual with a formally huge TS. A significance is meaningful only under a
correct null model, and ours is demonstrably not. This is exactly the trap
Totani's own systematics sections (§3.4) exist for, and exactly the
plausible-looking wrong number rule 1 warns about.

**The spectral shape does NOT show a 20 GeV peak.** A crude SED proxy
(halo counts × E_c, valid up to slowly-varying exposure for log-uniform bins)
peaks at **bin 2, ~4 GeV**, declining gently above — more reminiscent of a
GCE-like few-GeV spectrum than of Totani's 20 GeV feature. Read NOTHING into
this yet: the exposure-weighting convention is unvalidated against any
published number, the background inadequacy could shape the recovered spectrum
arbitrarily, and the halo template only absorbs a SMALL part of the residual
structure (mean|σ| improves by ~2–8% per bin) — most of what the background
fails to describe is NOT halo-shaped.

## Before any of this can be believed

1. Validate the pipeline against Totani's own background spectra (his Fig. 2/3)
   — the exposure/units conventions have never been checked against a published
   number.
2. His systematics program: quadrant splits (§3.4.1), the PS-masked variant
   (§3.4.3), bubble-template variants (§3.4.2).
3. Understand the residual structure the halo does NOT absorb — until the
   background describes the data at the few-σ level in the halo region, no
   component added on top can be interpreted.
