# Halo search, corrected estimator and templates — v3

> ⚠ **SUPERSEDED 2026-09-07.** The baseline bubble template changed to Oscar
> Macias' structured map (PI decision, physical grounds). The re-baselined
> table is `systematics_2026-09-07.md`; the 20.8 GeV baseline is now 4,612
> [4,026, 5,237], TS 24, with no excess above 35 GeV. This document is kept
> as the record of the flat-mask baseline. Its validation content stands.

**Supersedes** v1, v2 and the quadrant JSON (all void per
`validation/audit_2026-09-04.md`). **Data:** `runs/fff2de5bc7e0`.
**Full output:** `results/fits/halo_fit_v3_2026-09-04.json`.

## What changed since the void results

1. **Masked sources dropped from the point-source template** (audit finding 2).
   240 of 2,651 ROI sources sat inside extended-source mask circles carrying
   29–45% of predicted point-source flux while their photons were masked out
   of the data. Fixing it moved the fitted point-source normalization from
   **0.65–0.96 to 0.91–1.17** — now straddling unity, a materially better
   absolute validation than before.
2. **Profile MLE replaces block alternation** (audit finding 1). Acceptance
   gate is an injection test, and it passes: **slope 0.998–1.012** baseline,
   0.975–1.129 per quadrant, against 0.355–0.636 for the old estimator.
3. **The disk is excluded**, which is Totani's actual procedure (§3.2) and was
   wrong in every earlier fit here.

## Baseline (disk-excluded)

| bin | E (GeV) | halo counts | 1σ | fraction | TS | mean\|σ\| |
|---|---|---|---|---|---|---|
| 0 | 1.51 | 150,226 | [145.6k, 153.1k] | 8.1% | 390 | 4.34 |
| 1 | 2.55 | 91,952 | [90.6k, 96.7k] | 9.5% | 317 | 3.20 |
| 2 | 4.31 | 50,570 | [48.7k, 51.9k] | 10.9% | 201 | 2.26 |
| 3 | 7.28 | 26,187 | [25.2k, 27.4k] | 12.5% | 142 | 1.79 |
| 4 | 12.29 | 14,660 | [12.9k, 15.7k] | 15.4% | 61 | 1.31 |
| 5 | 20.76 | 7,705 | [7.3k, 8.3k] | 17.8% | 69 | 1.09 |
| 6 | 35.06 | 3,843 | [3.1k, 4.2k] | 18.7% | 21 | 0.93 |
| 7 | 59.22 | 1,494 | [1.3k, 1.7k] | 16.8% | 19 | 0.90 |
| 8 | 100.02 | 555 | [416, 696] | 15.2% | 5 | 0.78 |

The SED proxy (counts × E) is **flat to within ~30% from 2.5 to 60 GeV** and
falls at both ends — no 20 GeV peak, and no few-GeV peak either now that the
disk is out. The halo fraction rises monotonically with energy to ~19%.

> ⚠ **See `systematics_2026-09-04.md` (same day).** The reading below --
> that the quadrant failure points at background contamination -- stands,
> but the follow-up attribution to the bubble template specifically was an
> OVER-READ: removing the bubble makes the quadrants agree *and doubles the
> halo*, i.e. the halo absorbs the bubbles' real photons. Quadrant
> consistency is necessary, not sufficient.

## 🚨 The quadrant split does NOT support a halo interpretation

Totani's own §3.4.1 systematic is the sharpest test available here, and it
fails badly. A genuine, roughly spherical halo must appear in all four
quadrants at comparable strength. Instead:

| bin | baseline | l+b+ | l−b+ | l+b− | l−b− |
|---|---|---|---|---|---|
| 0 | 150,226 | 23,172 | **−16,616** | 45,223 | 31,451 |
| 1 | 91,952 | 16,228 | **−13,846** | 28,247 | 19,584 |
| 5 | 7,705 | 1,567 | 730 | 2,372 | 1,836 |

- **One quadrant (l−b+) goes NEGATIVE** at low energy — the fit prefers to
  *remove* halo-shaped emission there while adding it elsewhere. A real halo
  cannot do that.
- **Quadrant spread is a factor of ~3** at every energy where statistics
  allow, far outside the quoted 1σ intervals.
- **The quadrants sum to 0.46–0.84 of the baseline**, not ~1. The baseline
  amplitude is not a sum of consistent parts; it is partly absorbing a
  large-scale gradient the background model cannot describe.

The natural reading is that the "halo" component is largely soaking up
residual asymmetric structure — plausibly Loop I / bubble mismodeling, which
is exactly north–south and east–west asymmetric — rather than an isotropic
halo. That is consistent with the residuals: mean |σ| is still 4.3 per cell at
bin 0, so the background is demonstrably inadequate at low energy where the
quadrant disagreement is worst. Note the disagreement shrinks with energy
(quadrant sum ratio 0.55 → 0.84 from bin 0 to bin 5) as residuals approach
Poisson.

## Status

The pipeline now passes every check it has: absolute normalization (point
sources ≈ 1, isotropic 0.9–1.5× IGRB), estimator unbiasedness (injection slope
≈ 1), and internal consistency. **The physics conclusion is negative so far:
no evidence for a coherent 20 GeV halo, and positive evidence that the fitted
amplitude is contaminated by asymmetric background residuals.** Before this is
a result rather than a status report: the bubble/Loop I systematic variants
(§3.4.2), the PS-masked variant (§3.4.3), and ideally a second diffuse family.
