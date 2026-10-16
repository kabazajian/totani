# Systematics: the halo amplitude is not robust — it is systematics-dominated

> ⚠ **SUPERSEDED 2026-09-07.** The baseline bubble template changed to Oscar
> Macias' structured map (PI decision, physical grounds). The re-baselined
> table is `systematics_2026-09-07.md`; the 20.8 GeV baseline is now 4,612
> [4,026, 5,237], TS 24, with no excess above 35 GeV. This document is kept
> as the record of the flat-mask baseline. Its validation content stands.

**Data:** `runs/fff2de5bc7e0`, disk-excluded (Totani §3.2), profile estimator
(injection slope 1.00–1.04 in every variant). Full output:
`results/fits/systematics_2026-09-04.json`.

## Halo amplitude vs background choice (counts)

| variant | bin 0 (1.5 GeV) | bin 2 (4.3) | bin 5 (20.8) | bin 8 (100) |
|---|---|---|---|---|
| baseline | 150,226 | 50,570 | 7,705 | 555 |
| no Loop I | 162,341 | 54,182 | 7,818 | 555 |
| no bubble | 242,300 | 90,906 | **15,242** | 1,609 |
| no Loop I + no bubble | 227,762 | 88,497 | 14,961 | 1,604 |
| no dust | 152,649 | 57,794 | 9,111 | 679 |
| PS masked (95%) | **3,921** | 29,484 | 7,424 | 464 |
| PS masked (1° fixed) | 28,250 | 10,779 | **2,931** | 89 |

🚨 **At 20.8 GeV — the energy Totani claims — the fitted halo ranges over a
factor of 5 (2,931 to 15,242) across defensible background choices, while the
statistical 1σ interval is ~±10%.** At 1.5 GeV the range is a factor of 60.
The systematic uncertainty dominates the statistical one by more than an
order of magnitude, so no amplitude here supports a claim in either direction.

## Quadrant consistency (sum of quadrants / baseline; 1.0 = consistent)

| variant | bin 0 | bin 2 | bin 5 |
|---|---|---|---|
| baseline | 0.55 (one negative) | 0.76 | 0.84 |
| no Loop I | 0.68 (negative) | 0.90 | 0.91 |
| no bubble | 0.93 | 1.06 | 0.93 |
| no Loop I + no bubble | 1.10 | 1.10 | 0.95 |
| no dust | −0.05 (negative) | 0.16 (negative) | 0.68 |
| PS masked (95%) | 0.87 (negative) | 0.71 | 1.05 |
| PS masked (1°) | 1.15 | 0.76 | 0.79 |

## ⚠ Correcting an over-read

An earlier note said removing the bubble template "restores quadrant
consistency," implying the bubble template was *the* contaminant and the rest
was clean. That reading is wrong, and the amplitude column is why: **removing
the bubble makes the quadrants agree AND doubles the halo** (7,705 → 15,242 at
bin 5).

The bubbles are real emission; deleting their template does not delete their
photons. The halo — the only other broad, centrally-brightening component —
absorbs them. And because the bubbles are two lobes roughly symmetric about
the GC in both hemispheres, a halo absorbing them looks *quadrant-consistent*.
So the `no_bubble` configuration is not a cleaner measurement; it is a
demonstration that **the halo template can masquerade as the Fermi bubbles.**
Quadrant consistency is a necessary condition for a halo, not a sufficient one.

## What the point-source variants show

Masking sources instead of modelling them (Totani §3.4.3) removes most of the
low-energy signal: bin 0 falls from 150,226 to 3,921 (95% radii) — a factor of
38. So the bin-0 "halo" was largely point-source-associated residual, not a
diffuse component. At 20.8 GeV the 95%-radius variant agrees with the baseline
(7,424 vs 7,705) while the deliberately-oversized 1° mask gives 2,931; the 1°
mask removes far more sky than containment requires at that energy, so the two
bracket rather than contradict.

## Bottom line

- **No evidence for a 20 GeV halo.** The SED proxy is flat 2.5–60 GeV in the
  baseline and no variant produces a 20 GeV peak.
- **No evidence against one either.** The amplitude is not determined to better
  than a factor of several by this background model.
- **The limiting factor is the background model, exactly as Totani (§3.4.4)
  and Stenhouse et al. both state.** One diffuse family cannot settle this.

That is precisely the gap this project exists to close: the categorical
diffuse-model marginalization needs the |b| ≤ 60 Galp21 delivery (6 GALPROP
realizations) to turn a spread-across-arbitrary-choices into a marginalized
uncertainty. Until then these numbers are a pipeline demonstration, not a
measurement.
