# Oscar's two questions of 2026-09-18, answered on the stack

`scripts/bubble_lowE_check.py` → `results/fits/bubble_lowE_check_2026-09-18.json`.
315 s, peak RSS 7.1 GB. Structured (`macias`) baseline, cell space.

## Context: what Oscar said

The structured map is **Figure 12 of Ackermann et al. 2017 (arXiv:1704.03910)**,
not the 2014 bubbles paper as every note here had said (corrected in the code
and notes 2026-09-18). Verified against the paper: §5.2 derives it from the
**1–10 GeV** residuals *after subtracting gas-correlated emission and point
sources* (their GALPROP 5-ring gas, 3FGL sources) by a three-component
spectral decomposition — soft E^−2.4, medium E^−1.9 "bubble-like", hard
MSP-like — with a **1.5σ significance cut** on the medium component; Oscar's
inpainting fills the point-source holes. So it should be at its *best* at low
energy, and above 10 GeV its morphology is an extrapolation (motivated by the
2014 finding that the bubble spectrum is uniform and hard). He wondered
whether the poorer fit at 1.5–2.6 GeV "might be related to some other
templates", and proposed a SkyFACT-style controlled flexibility of the bubble
template, tested on simulated skies with and without an injected halo.

## (a) Yes — it is Loop I

Δ(−2lnL), background-only, structured − flat (negative = structured better),
per energy bin:

| treatment | 1.5 | 2.6 | 4.3 | 7.3 | 12.3 | 20.8 | 35 | 59 | 100 |
|---|---|---|---|---|---|---|---|---|---|
| baseline (PS modelled) | **+243** | **+45** | −54 | −100 | +48 | −10 | −54 | −26 | +7 |
| no dust residuals | +109 | −8 | −99 | −126 | +33 | −20 | −60 | −31 | +4 |
| **no Loop I** | **−215** | **−258** | **−216** | **−180** | +8 | −25 | −57 | −27 | +6 |
| no dust, no Loop I | −334 | −317 | −269 | −214 | −12 | −37 | −65 | −33 | +3 |
| PS masked (95%) | +93 | +33 | +1 | −49 | +67 | −6 | −63 | −26 | +8 |
| PS masked (1° fixed) | +135 | +28 | −23 | −44 | +35 | −19 | −38 | −14 | +10 |

🚨 **With the two Loop I shells removed, the structured map is the better
description at every energy below 12 GeV, by a wide margin.** The +243 at
1.5 GeV — the number that has been carried since 09-07 as "the one place the
structured map fits worse" — is not a property of the map. It is what happens
when the map and our two geometric Loop I shells are both in the fit: a
1–10 GeV residual template built with a 1.5σ significance cut plausibly
contains Loop I / North Polar Spur emission (same sky, similar soft-ish
spectrum at these energies), and our shells then fight it. The flat mask,
having no structure there, coexists with the shells. Oscar's instinct that it
was another template was right, and the template is Loop I.

Point-source masking cuts the deficit to +93/+135 but does not remove it;
dropping the dust residuals takes it to +109. Neither is the main effect.

Above 12 GeV nothing changes in any row: the structured map is modestly
better from 20 GeV up regardless. So the 35+ GeV extrapolation Oscar flagged
is not where the trouble is, on this evidence.

**Consequence for the halo.** Loop I has been a recurring suspect (the
|b| = 45–60° residual was attributed to it and withdrawn on 09-06;
`no_loopI` moves the 1.5 GeV halo 106,612 → 130,842 and TS 217 → 383 but
leaves 20.8 GeV at 4,612). The right test now is **structured bubbles with no
Loop I shells as a baseline candidate**, scored the usual way (halo, TS,
quadrants, injection) — not run here because it is a baseline change, which
is Kev's call. It is one cell-space fit.

## (b) Regional bubble normalisations: unbiased, but not stable

⚠ **This is not SkyFACT, and Kev said so (2026-09-18).** SkyFACT (Storm,
Weniger & Calore 2017, `docs/Storm_SKYFACT_1705.04065v2.pdf`, §2) gives every
template a *per-pixel* spatial modulation τ_p and a per-bin spectral
modulation σ_b, fitted under a penalised Poisson likelihood with a
maximum-entropy regulariser (effective variance 1/√λ per parameter) and
explicit smoothing terms — of order 10⁵ nuisance parameters, convex, solved
with L-BFGS-B, errors from the Fisher matrix. What follows is a
hard-partitioned, unregularised split of one template into at most eight
regional normalisations on the 10° cell grid. It answers a narrower question
— *does giving the bubble template piecewise freedom move the halo or bias
the estimator?* — and says nothing about what a regularised pixel-level
modulation would do. The SkyFACT question is Oscar's and Kev's to take up;
this section is not an answer to it.

The structured template given 1, 2 (N/S), 4 (quadrants) and 8 (quadrants ×
|b| ≤ 30° / > 30°) independent normalisations per bin:

| bubble norms | halo @ 20.8 | 68% | TS | quadrant Σ/joint | injection slope (min / bin 5) |
|---|---|---|---|---|---|
| 1 | 4,612 | [4,026, 5,237] | 23.5 | 0.57, one negative | 1.00 / 1.00 |
| 2 (N/S) | 4,443 | [3,844, 5,065] | 22.1 | 0.59, one negative | 1.00 / 1.00 |
| 4 (quadrants) | 4,387 | [3,794, 5,002] | 21.4 | 0.60, one negative | 1.00 / 1.00 |
| **8 (quadrants × lat)** | **6,637** | [5,901, 7,449] | 36.9 | 0.63, one negative | 1.00 / 1.00 |

Halo per bin at 8 norms: 145,380 / 85,654 / 42,744 / 18,822 / 13,045 /
6,637 / 1,068 / 822 / 166 — every bin up 40–300% on the baseline.

Three things:

1. **The injection gate passes at every rung** (slope 1.00 at 20.8 GeV, worst
   1.04 at 1.5 GeV). Splitting the bubble template into eight regional pieces
   does *not* let it absorb an injected spherical halo. The estimator stays
   unbiased. This is the half of Oscar's simulated-sky test that asks "does
   flexibility eat a real signal?", answered on the real stack: at this level
   of freedom, no.
2. **But the answer is not stable.** 1 → 4 norms moves the halo by −5%;
   the eighth rung moves it by **+44%**, outside the 68% interval, with TS
   rising to 37. An unbiased estimator giving a different answer under a
   different bubble model is a *model* systematic, not an estimator one — and
   it is the other half of Oscar's test, the one that needs simulated skies
   with different bubble morphologies, because the real sky cannot say which
   rung is "right". (Δ(−2lnL) per rung was not recorded; it should be, so the
   extra parameters can be judged.) A regularised method would replace the
   arbitrary choice of rung with a hyper-parameter — which is the point of
   SkyFACT, and why this split is not a substitute for it.
3. **The quadrant test never passes.** Σ/joint 0.57 → 0.63, the l < 0, b > 0
   quadrant negative at every rung. Regional bubble freedom is not what the
   quadrants disagree about either. After today, that failure has survived:
   the bubble template choice, N/S, quadrant and 8-way bubble splits, the
   eROSITA outer template, and every 09-07 variant. **It is the sharpest thing
   in this analysis and nothing yet touches it.**

## On the SkyFACT proposal

Not assessed here. The regional split above shows that the halo is sensitive
to how much freedom the bubble template has — a model sensitivity, with the
estimator unbiased — which is a reason the question matters, not an answer to
it. Oscar's simulated-sky test (different bubble morphologies, with and
without an injected halo) is the right instrument for a regularised
pixel-level modulation, and it is his and Kev's to design. The step this
note does support is the Loop I one above, which is a baseline decision.

Analysis and text: Claude Fable 5.1 (`claude-fable-5-1[1m]`), under Kev's
direction.
