# Structured bubbles + no Loop I as baseline: tested, and NOT recommended — 2026-09-18

Kev's decision #3 of 2026-09-18 ("structured + no Loop I as baseline") was
made on the framing in `bubble_lowE_check_2026-09-18.md` — that the structured
map's poor fit at 1.5–2.6 GeV "is Loop I" — and the candidate was run in full:
`systematics.py --no-loopI` → `fits/systematics_2026-09-18.json`;
`boundary_test.py --no-loopI` → `fits/boundary_test_2026-09-18_noLoopI.json`;
`boundary_breakdown.py --no-loopI` → `fits/boundary_breakdown_2026-09-18_noLoopI.json`
(reproduces its parent to 0.0); `loopi_amplitudes.py` →
`fits/loopi_amplitudes_2026-09-18.json`. **The last of these shows the
framing was incomplete, and the analysis's recommendation is to keep Loop I.
The candidate is recorded here and not adopted; the 09-07 baseline stands
unless Kev rules otherwise.**

## What the candidate gives

| | structured + Loop I (09-07 baseline) | structured, no Loop I |
|---|---|---|
| halo @ 20.8 GeV | 4,612 [4,026, 5,237], TS 23.5 | **4,612** [4,127, 5,332], TS 26.8 |
| halo @ 1.5 GeV | 106,612, TS 217 | **130,842** (+23%), TS 383 |
| halo @ 4.3 GeV | 31,305 | 35,519 (+13%) |
| above 35 GeV | none (TS ≤ 0.8) | none (TS ≤ 1.2) |
| quadrants @ 20.8 | [783, −313, 1,446, 707], Σ 0.57 | [1,085, **−566**, 1,281, 707], Σ 0.54 |
| quadrants @ 1.5 | one negative | [36,514, **−85,074**, 33,917, 28,198] |
| injection slopes | 1.00–1.01 | 1.00–1.09 |
| envelope @ 20.8 (physical variants) | 4,326–15,242 (3.5×) | 4,474–14,961 (3.3×) |
| boundary excess @ 20.8 (structured) | +404 ± 200 (2.0σ) | **+688 ± 201 (3.4σ)** |
| boundary excess @ 1.5 | +20,034 (14.7σ) | **+41,007 (29.9σ)** |

The 20.8 GeV number does not move. Everything else moves the wrong way.

## Why: the shells are not double-counted by the structured map

Fitted Loop I counts over the |b| ≥ 10 cells, background-only
(`loopi_amplitudes_2026-09-18.json`):

| E (GeV) | 1.5 | 2.6 | 4.3 | 7.3 | 12.3 | 20.8 | 35 |
|---|---|---|---|---|---|---|---|
| beside the **flat** mask | 139,782 | 80,934 | 33,975 | 14,184 | 7,508 | 2,860 | 1,177 |
| beside the **structured** map | 128,227 | 73,108 | 29,907 | 12,171 | 6,441 | 2,431 | 988 |
| ratio | 0.92 | 0.90 | 0.88 | 0.86 | 0.86 | 0.85 | 0.84 |

If the structured map (a 1–10 GeV residual template) contained the Loop I /
North Polar Spur emission, the shells would collapse beside it. They drop
8–15%. And what the two shell parameters buy in −2lnL, beside the structured
map: **−2,986 at 1.5 GeV, −1,586, −664, −224, −75, −23 at 20.8 GeV** (beside
the flat mask: −3,444 … −38). They model ~10⁵ counts at 1.5 GeV that no other
template in the stack carries.

So the correct reading of `bubble_lowE_check` part (a) is narrower than
"it is Loop I": *with the shells present*, flat + shells describes 1.5–2.6 GeV
better than structured + shells by 243 and 45; *without them*, structured
beats flat by 215–258. That is an interaction — the structured map and the
geometric shells disagree about the shape in the region where they overlap
(the northern sky at l < 0, which is Loop I territory: centre (l, b) ≈ (−31°,
+18°), radius 58°) — and it is worth 243 units at 1.5 GeV. Dropping the shells
costs 2,986 units at the same energy. The interaction is a **template-shape
mismatch in the overlap**, not a redundancy.

## What dropping them does to the sky tests

The +404 → +688 at 20.8 GeV and 20,034 → 41,007 at 1.5 GeV in the boundary
test is the emission the shells had been modelling inside the northern
eROSITA lobe reappearing as "unmodelled". In the breakdown it lands where
Loop I is: l < 0 (+533 ± 123 vs +100 at l > 0), and |b| 45–60° comes back
(+150 ± 90, was −1). The quadrant that goes negative is l < 0, b > 0 — the
Loop I quadrant — and it goes *more* negative without the shells (−313 →
−566). All three sky tests say the same thing: the shells were doing work.

## Recommendation

**Keep structured + Loop I as the baseline.** Treat the +243 as what it is: a
disagreement between two templates about the low-energy shape of the northern
l < 0 sky, which is also where the boundary excess sits and where the quadrant
test fails. Two cheap next steps address that directly, neither of which is
"drop a template":

1. ~~Swap the two geometric shells for the Pohl `LoopI_Normalized.fits`.~~
   **Checked 2026-09-18 evening, at Kev's instruction, and withdrawn:** that
   file is a 40°×40° cutout (|l|, |b| ≤ 20°, 0.1°), nonzero everywhere, and a
   least-squares fit to our two Wolleben shells gives R² = 0.9994 with equal
   weights on both shells — it *is* the same geometric model. No empirical
   Loop I map exists on this machine. Candidates would be the Haslam 408 MHz
   all-sky map (Remazeilles+ 2015 destriped/desourced, LAMBDA) or the Planck
   synchrotron map, sampled to the 0.125° grid; radio synchrotron traces the
   electrons rather than their gamma-ray IC, so this is a modelling choice,
   not a drop-in. Kev's call.
2. **Oscar's eROSITA brightness template** (next week) in the same region.

And the SkyFACT-style question Oscar raised is, on this evidence, exactly a
question about this overlap region — which is an argument for it, not
against.

Analysis and text: Claude Fable 5.1 (`claude-fable-5-1[1m]`), under Kev's
direction. Peak RSS: systematics 7 GB, boundary test 8.3 GB, breakdown
12.4 GB.
