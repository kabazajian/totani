# Does the eROSITA-bounded residual move the halo? — 2026-09-18

`scripts/erosita_variant.py` → `results/fits/erosita_variant_2026-09-18.json`.
Kev's decision (2026-09-18) on open item 14 of `STATE_2026-09-17`: measure it.
Run on the structured (`macias`) baseline; 463 s; peak RSS 6.4 GB.

## What was added

One template: the region **inside Oscar's eROSITA boundary and outside the
Ackermann flat Fermi-bubble footprint** — Zhang et al.'s "outer region", the
region where `boundary_test_2026-09-08` found +404 ± 200 counts at 20.8 GeV
that the baseline model does not account for — as a flat-intensity shape with
a free normalisation per bin, through the same pipeline as every other diffuse
template (`stack.cell_template_from_shape`: exposure × solid angle,
extended-source mask, cells). Also as two templates split north/south, since
the Fermi bubbles themselves want separate N/S normalisations at ~3σ
(`bubble_NS`, adopted by Kev the same day).

`cell_template_from_shape` is validated in this run: fed the flat footprint, it
reproduces `systematics_2026-09-07`'s `bubble_flat` variant exactly — 7,705
[7,306, 8,290], TS 69.1, quadrant sum 0.84.

## The answer: the halo barely moves at 20.8 GeV

| variant | halo @ 20.8 GeV | 68% | TS | vs baseline |
|---|---|---|---|---|
| **baseline (structured)** | **4,612** | [4,026, 5,237] | 23.5 | — |
| + eROSITA outer | 4,556 | [3,935, 5,178] | 19.9 | **−1.2%** |
| + eROSITA outer, N/S | 4,556 | [3,914, 5,126] | 22.1 | **−1.2%** |
| bubble N/S | 4,443 | [3,844, 5,065] | 22.1 | −3.7% |
| bubble N/S + eROSITA N/S | 4,443 | [3,812, 5,032] | 20.9 | −3.7% |
| flat mask | 7,705 | [7,306, 8,290] | 69.1 | +67% |
| flat mask + eROSITA outer | 7,818 | [7,254, 8,364] | 52.8 | +70% |

✅ **The eROSITA residual is not a halo systematic at 20.8 GeV.** A −56 count
shift against a ±600 statistical interval. It does not need to enter the
envelope (4,326–15,242); if quoted, it sits inside the baseline's own error
bar. The flat-mask variant is unaffected too (+113).

At **4–7 GeV** the single-normalisation template does move the halo, by −19%
and −24% (31,305 → 25,285; 13,639 → 10,366). The N/S version does not (−4%).
The single-norm template is the worse description there (see Δ(−2lnL) below),
so the −19/−24% is a degeneracy artefact of one flat normalisation over a
region whose emission is neither flat nor N/S symmetric (see the breakdown),
not a sky statement.

Injection slopes 0.97–1.11 (the 1.11 is bin 0 of the doubly-split variant, a
mild degeneracy flag; everything at 20.8 GeV is 1.00).

## What the template absorbs, and where — the real finding

Background-only fit, counts over the |b| ≥ 10 cells:

| E (GeV) | outer (one norm) | outer **north** | outer **south** | Δ(−2lnL) one norm | Δ(−2lnL) N/S |
|---|---|---|---|---|---|
| 1.51 | 14,046 | **34,107** | **0** | −59 | **−571** |
| 2.55 | 8,180 | 18,164 | 0 | −39 | −316 |
| 4.31 | 6,296 | 10,246 | 81 | −48 | −211 |
| 7.28 | 3,688 | 4,591 | 580 | −37 | −91 |
| 12.29 | 575 | 1,995 | 0 | −2 | −42 |
| **20.76** | **543** | **708** | **62** | −4 | **−11** |
| 35.06 | 136 | 120 | 0 | 0 | −1 |
| 59.22 | 118 | 73 | 47 | −1 | −1 |

⚠ **The fit puts it all in the north — and that is a degeneracy statement,
not a sky statement.** Given the freedom, the cell-space fit assigns
essentially all of the outer-region emission to the north template and none
to the south at every energy, and prefers the N/S split by Δ(−2lnL) = −571 at
1.5 GeV, −11 at 20.8 GeV (2 extra parameters). But the direct measurement —
the matched-θ contrast in `boundary_breakdown_2026-09-18`, which does not go
through 10° cells or a 19-template fit — finds the 1.5 GeV excess **74% in the
south** (S +14,788 ± 871, N +5,246 ± 1,052) and the 20.8 GeV excess in both
hemispheres (N +230 ± 152, S +173 ± 130). So the fit's N/S assignment is about
which template shape the 10° cell fit can use: in the south the same emission
is absorbed by `bubble_S` (52,715 → 53,707 at 1.5 GeV when the outer templates
are added) and the rest of the stack, and the north template is left holding
what nothing else fits. Read the hemisphere split from the breakdown, not from
here. The 543–708 counts at 20.8 GeV are consistent with the boundary test's
+404 ± 200 (same region, different estimator).

With the **flat** mask the outer template absorbs 1,518 counts at 20.8 GeV
(vs 543 on the structured baseline), the same 2:1 ratio the boundary test
found (806 vs 404). But flat + outer template does **not** reproduce the
structured result — the halo stays at 7,818, not 4,612. So what the structured
map does to the halo is not "it covers the outer region"; it is the structure
*inside* the Fermi-bubble footprint.

## The quadrant test does not improve

| variant | bin 5 quadrants [l+b+, l−b+, l+b−, l−b−] | Σ/joint |
|---|---|---|
| baseline | [783, **−313**, 1,446, 707] | 0.57 |
| + eROSITA outer | [768, −313, 1,446, 585] | 0.55 |
| + eROSITA N/S | [783, −313, 1,446, 870] | 0.61 |
| bubble N/S | [768, −313, 1,446, 707] | 0.59 |
| bubble N/S + eROSITA N/S | [813, −313, 1,446, 694] | 0.59 |

The negative quadrant (l < 0, b > 0: north, toward negative longitude) is
**−313 in every variant**, and the sum stays at 0.55–0.61 of the joint. Adding
a template for the eROSITA outer region — even one that lives in the north —
does nothing to it. **The quadrant failure (open item 7) and the boundary
excess are not the same problem at 20.8 GeV.** Whatever the l−b+ quadrant is
disagreeing about, it is not emission bounded by the X-ray edge. (At 1.5 GeV
both north quadrants are negative and the added template makes them more so;
the low-energy quadrant failure is a different, larger animal.)

## Bottom line

- **Item 14 closed:** the eROSITA residual moves the 20.8 GeV halo by −1.2%,
  inside the statistical interval. Not a systematic. Quotable as such.
- **Item 9:** N/S bubble normalisations move the halo −3.7% (4,443) and are
  preferred at Δ(−2lnL) = −9 at 20.8 GeV, −125 at 1.5 GeV. Kev adopted them.
  Re-baselining to 4,443 is a separate decision (see STATE).
- **Where it is** is answered by `boundary_breakdown_2026-09-18`, not by the
  template fit: at 20.8 GeV, 10–20° outside the Fermi-bubble edge, at
  negative longitude, in both hemispheres, at |b| = 10–45°.
- **Item 7 is independent of all of this** and remains the sharpest open
  problem.

Analysis and text: Claude Fable 5.1 (`claude-fable-5-1[1m]`), under Kev's
direction; 138 tests passing.
