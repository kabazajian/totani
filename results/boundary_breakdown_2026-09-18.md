# Where inside the eROSITA edge the boundary excess sits — 2026-09-18

`scripts/boundary_breakdown.py` → `results/fits/boundary_breakdown_2026-09-18.json`.
Open item 15 of `STATE_2026-09-17`, run at Kev's request 2026-09-18. Both
bubble templates; 233 s; **peak RSS 12.2 GB** (the largest job in this repo).

## The construction

The inside region of `boundary_test_2026-09-08` — inside Oscar's eROSITA
boundary, outside the flat Fermi-bubble footprint, |b| ≥ 10°, per hemisphere —
is cut three ways: by great-circle distance from the Fermi-bubble edge
(0–5°, 5–10°, 10–20°, >20°), by sign of longitude, and by latitude band
(10–30°, 30–45°, 45–60°). Each piece is matched in θ_GC against the **same full
outside region** (`erosita.matched_theta_weights_regions`), so all pieces
refer to one common reference and are comparable; the statistic is also given
per matched square degree. The `all` piece is the whole inside and
**reproduces the 09-08 result to 0.0 relative** on both templates — the
generalised weights are the old weights.

The fitted halo's own predicted difference is ≤ 0.4 counts in every piece at
1.5 GeV and 0.0 at 20.8 GeV: the matching is blind to a spherical component in
every sub-region, not only in the whole.

Sub-regions near the Fermi edge lose most of their area to the θ < 34.5° hole
(edge 0–5°: 435 of 781 deg² dropped in the north). Their error bars are
correspondingly larger; that is the honest size of the question.

## 20.8 GeV, structured baseline (`macias`), N + S combined

lobe = inside − θ-matched outside, minus the fitted halo's predicted
difference. Flat-mask (`ackermann`) result alongside.

| piece | structured | | flat | |
|---|---|---|---|---|
| **all** | **+404 ± 200** | 2.0σ | +806 ± 200 | 4.0σ |
| edge 0–5° | **−120 ± 81** | −1.5σ | **+350 ± 78** | **4.5σ** |
| edge 5–10° | +86 ± 76 | 1.1σ | +189 ± 75 | 2.5σ |
| **edge 10–20°** | **+332 ± 100** | **3.3σ** | **+322 ± 101** | **3.2σ** |
| edge > 20° | +55 ± 56 | 1.0σ | +76 ± 56 | 1.3σ |
| **l < 0** | **+415 ± 123** | **3.4σ** | +674 ± 122 | 5.5σ |
| l > 0 | −38 ± 121 | −0.3σ | +198 ± 121 | 1.6σ |
| \|b\| 10–30° | +295 ± 130 | 2.3σ | +351 ± 131 | 2.7σ |
| \|b\| 30–45° | +220 ± 100 | 2.2σ | +393 ± 99 | 4.0σ |
| \|b\| 45–60° | **−1 ± 90** | 0.0σ | +229 ± 89 | 2.6σ |

Per hemisphere (structured): l < 0 is N +257 ± 102, S +158 ± 68; edge 10–20°
is N +222 ± 78, S +110 ± 62; |b| 10–30° is N +295 ± 87, S 0 ± 97.

### What this says

1. **It does not hug the Fermi-bubble edge.** Within 5° of the edge the
   structured baseline leaves *nothing* (−120 ± 81; −129 ± 56 in the north, a
   slight over-subtraction). The flat mask leaves +350 ± 78 there (4.5σ) —
   **that** is spillover from a hard-edged template, and the structured map
   removes it. This is the mechanism behind the 09-08 finding that the
   structured map absorbs half the excess: it absorbs the part at the edge.

2. **What survives sits 10–20° outside the Fermi edge**, +332 ± 100 (3.3σ),
   and it is **the same on both templates** (+322 flat). Whatever it is, the
   bubble template does not touch it, so it is not bubble-edge modelling.

3. **It is entirely at negative longitude**: +415 ± 123 at l < 0, −38 ± 121 at
   l > 0 — and in both hemispheres (N +257, S +158). The eROSITA lobes lean to
   negative l in the north (⟨l⟩ = −9.5°) but positive in the south (+3.8°), so
   this is not the footprint's asymmetry; it is an east–west asymmetry in the
   sky the model leaves behind.

4. **It is at |b| = 10–45°, not 45–60°** (−1 ± 90 at high latitude). So
   **open item 10 — the antisymmetric residual at |b| = 45–60° — is not this,
   at 20.8 GeV.** At 1.5 GeV that band does carry +2,140 ± 594 (3.6σ), all
   southern (S +2,832, N −692), which is the sign pattern item 10 reported;
   the two may overlap at low energy but they are not one thing at 20 GeV.

5. **Open item 7.** The quadrant whose halo amplitude is negative in every
   variant is l < 0, b > 0. The 20.8 GeV excess is at l < 0 in both hemispheres
   — so the north-west quadrant contains +257 ± 102 counts of unmodelled
   emission 10–20° outside the bubble edge. Yet adding a template for exactly
   that region (`erosita_variant_2026-09-18`) left the quadrant at −313. The
   two problems share a quadrant and are not the same problem; that is now
   established rather than suspected.

## 1.5 GeV, structured baseline

| piece | N | S | combined |
|---|---|---|---|
| all | +5,246 ± 1,052 (5.0σ) | **+14,788 ± 871 (17.0σ)** | +20,034 (14.7σ) |
| edge 0–5° | **−2,355 ± 378 (−6.2σ)** | +3,912 ± 350 (11.2σ) | +1,558 (3.0σ) |
| edge 5–10° | +989 ± 370 | +4,636 ± 324 (14.3σ) | +5,625 (11.4σ) |
| edge 10–20° | +2,895 ± 548 (5.3σ) | **+9,124 ± 419 (21.8σ)** | +12,019 (17.4σ) |
| edge > 20° | **+3,768 ± 310 (12.2σ)** | +25 ± 240 | +3,793 (9.7σ) |
| l < 0 | +6,415 ± 708 (9.1σ) | +7,807 ± 446 (17.5σ) | +14,222 (17.0σ) |
| l > 0 | −1,135 ± 555 (−2.0σ) | +8,989 ± 587 (15.3σ) | +7,855 (9.7σ) |
| \|b\| 10–30° | **+9,252 ± 619 (14.9σ)** | +4,417 ± 668 (6.6σ) | +13,669 (15.0σ) |
| \|b\| 30–45° | +988 ± 519 | **+9,720 ± 447 (21.8σ)** | +10,708 (15.6σ) |
| \|b\| 45–60° | −692 ± 460 | +2,832 ± 376 (7.5σ) | +2,140 (3.6σ) |

At low energy the excess is **pervasive and mostly southern** — 74% of it in
the south, spread across every distance band and every latitude band there.
The north is different: the structured map *over*-subtracts at the northern
edge (−6.2σ within 5°), the northern excess is at low latitude (10–30°) and
far from the bubbles (> 20°, +3,768, where the south has none), and the
northern l > 0 quarter is slightly negative. This is not one component with
one shape; it is the low-energy foreground being wrong in several places at
once, which is what the 15σ-falling-to-nothing spectrum already said.

⚠ **This corrects the cell-space fit's reading.** `erosita_variant_2026-09-18`
assigned all of the outer-region emission to a *north* template; that was the
19-template cell fit choosing where to put it, not where it is. The direct
contrast here is the measurement.

## Flat-mask spillover, quantified

edge 0–5°, flat mask: +7,721 ± 509 at 1.5 GeV, +350 ± 78 at 20.8 GeV, both
hemispheres. Structured: +1,558 ± 515 and −120 ± 81. The flat mask's edge
excess is 30% of its total at 1.5 GeV and 43% at 20.8 GeV; the structured map
removes it. That is the on-the-sky content of "the structured template is the
more physical one".

## Bottom line

- The 20.8 GeV boundary excess is **not bubble-edge spillover** (the structured
  map has already removed that), is **10–20° outside the Fermi edge**, is
  **entirely at l < 0** in both hemispheres, and is at **|b| = 10–45°**.
- Item 10 (|b| = 45–60°) is **not** this at 20 GeV; possible overlap at
  1.5 GeV only.
- Item 7 (quadrants) shares the l < 0, b > 0 quadrant with it and is
  **not** fixed by modelling it.
- The low-energy excess is pervasive and mostly southern; the north shows the
  structured map over-subtracting at the edge.
- An l < 0 / l > 0 asymmetry 10–20° outside the bubbles, flat in the bubble
  template choice, is the thing to look at next — and the eROSITA
  *brightness* map Oscar offered is the natural template to test against it.

Analysis and text: Claude Fable 5.1 (`claude-fable-5-1[1m]`), under Kev's
direction.
