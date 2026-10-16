# Inside vs outside the eROSITA edge, matched in θ_GC — 2026-09-08

`scripts/boundary_test.py`, on Oscar's boundary reconstruction delivered
2026-09-08. JSON: `results/fits/boundary_test_2026-09-08.json` plus two
variants. Run against **both** bubble templates, which is the point of it.

## What was asked and what this answers

`STATE_2026-09-07` §"What to do when they reply" specified this test: replace
the reflection construction with inside/outside-boundary regions at matched
θ_GC, at pixel level, run on both bubble templates, as *"the independent,
on-the-sky test of the baseline choice"*. That is what this is.

## The construction

Two regions — inside and outside the eROSITA footprint — both restricted to
|b| ≥ 10°, both with the **same fixed** Fermi-bubble footprint (the Ackermann
flat 0/1 mask) excluded, so that what varies between the two runs is only which
bubble template is *fitted*. Zhang et al.'s "outer region" is exactly this:
outside the Fermi Bubbles, inside the eROSITA Bubbles.

The two regions are different shapes, so unlike the reflection test they are
matched in θ_GC only in **distribution**: in each 0.25° annulus both sides are
normalised and given the same annulus weight. A spherical component then
contributes equally to both sums whatever its radial profile.

🚨 **The matching is response-weighted, not geometric.** The statistic sums
counts, and a spherical component contributes `f(θ) × exposure × Ω`, so the
annulus normalisation is by exposure × solid angle. Matching on solid angle
alone — the first version of this, on 2026-09-08 — left the fitted NFW halo
predicting **329 more counts outside than inside at 20.8 GeV against a signal
of 960**, a 34% leak of exactly the component the test must be blind to. It was
caught by watching `halo_predicted_diff`, which is now **+0.02 counts** and is
written to the JSON every run for that reason. The weights are therefore per
energy bin, because the exposure is.

## Geometry, and the limit it puts on what this can say

| | north | south |
|---|---|---|
| θ_GC with support | 34.5–75.5° | 34.2–73.5° |
| inside area | 2,142 deg² | 1,783 deg² |
| outside area | 946 deg² | 1,182 deg² |
| inside area **dropped**, no outside at that θ | 569 deg² | 510 deg² |

⚠ **There is no support inside θ_GC ≈ 34.5°** — every pixel of the halo region
closer than that is inside the footprint, so there is nothing to compare it
with. That is not a cut anyone chose; it is the geometry. It removes 1,079 deg²
of the inside region and, more to the point, removes the range where the halo
has most of its leverage: `band_forecast_2026-09-06.json` puts 54–61% of it in
the 10–30° vs 30–60° contrast. **This test bounds a lobe's morphology. It is a
weak constraint on the halo amplitude.** Our support agrees with Oscar's
independent 0.25° quadrature (34.62° / 34.38°) to better than one annulus.

## Result — structured baseline (`macias`), both hemispheres combined

lobe = inside − θ-matched outside, minus the fitted halo's own predicted
difference. Counts on the inside region's own scale.

| E (GeV) | lobe | σ | halo inside |
|---|---|---|---|
| 1.51 | +20,034 ± 1,365 | **+14.7** | 21,613 |
| 2.55 | +10,838 ± 978 | **+11.1** | 11,493 |
| 4.31 | +5,930 ± 672 | **+8.8** | 6,341 |
| 7.28 | +3,161 ± 448 | **+7.1** | 2,760 |
| 12.29 | +976 ± 300 | +3.3 | 1,759 |
| **20.76** | **+404 ± 200** | **+2.0** | 933 |
| 35.06 | +120 ± 137 | +0.9 | 54 |
| 59.22 | +146 ± 89 | +1.6 | 38 |
| 100.02 | −15 ± 57 | −0.3 | 14 |

**There is emission bounded by the eROSITA edge that neither bubble template
captures.** It is overwhelmingly a low-energy effect: 14.7σ at 1.5 GeV, falling
monotonically, 2.0σ at 20.8 GeV, and consistent with zero above 35 GeV. A
spherical component can leak at most 0.02 counts into the 20.8 GeV number
against a 200-count error bar, so this is not the halo bleeding through.

## The baseline choice, tested on the sky

This is what the test was built for. At every energy the structured template
**absorbs a substantial part of the excess that the flat mask leaves behind**:

| E (GeV) | flat (`ackermann`) | structured (`macias`) | reduction |
|---|---|---|---|
| 1.51 | +24,821 (18.2σ) | +20,034 (14.7σ) | 19% |
| 4.31 | +7,988 (11.9σ) | +5,930 (8.8σ) | 26% |
| 7.28 | +4,450 (9.9σ) | +3,161 (7.1σ) | 29% |
| 12.29 | +1,710 (5.7σ) | +976 (3.3σ) | 43% |
| **20.76** | **+806 (4.0σ)** | **+404 (2.0σ)** | **50%** |
| 35.06 | +342 (2.5σ) | +120 (0.9σ) | 65% |

✅ **The 2026-09-07 baseline decision is supported on the sky, not only by the
fit statistic.** Kev adopted the structured map on physical grounds; this is an
independent measurement, in a region chosen by an X-ray boundary rather than by
anything in the gamma-ray fit, and the structured map is better there at every
energy — by half at 20.8 GeV. That said, **neither template removes the excess**:
at 20.8 GeV the structured baseline still leaves 2.0σ.

## Systematic variants at 20.8 GeV (combined)

| variant | flat | structured |
|---|---|---|
| baseline (`min` annulus weight, bubbles cut) | +806 ± 200 (4.0σ) | +404 ± 200 (2.0σ) |
| `--theta-weight harmonic` (variance-optimal) | +813 ± 203 (4.0σ) | +362 ± 204 (1.8σ) |
| `--exclude-bubble none` | +777 ± 252 (3.1σ) | +582 ± 252 (2.3σ) |

The annulus weighting is irrelevant, as expected. Keeping the Fermi bubbles in
raises the structured result and lowers the flat one, which is the same
statement as the table above by another route.

## The limit on a lobe, and how it compares

95% upper limit on lobe / halo counts inside, at 20.8 GeV, structured baseline:
**0.86** combined (north 1.07, south 1.00). The reflection test's best was
*lobe < 61% of the halo in R3* (`patch_test_2026-09-06.md`).

⚠ **These two numbers are not competing measurements of the same thing** and
the reflection result is not "beaten" or "confirmed" by this. The reflection
test bounds a lobe inside one 20°×15° box; this bounds a lobe *shaped like the
eROSITA footprint* over 3,900 deg². This one also has a positive central value
rather than a null, so its UL is inflated by a real excess. Above 35 GeV the
ratio is meaningless — the fitted halo inside falls to tens of counts, so
UL/halo runs to 7–9 while the lobe itself is consistent with zero.

⚠ **The eight patches of the reflection test are 100% inside the eROSITA
footprint**, measured 2026-09-08. The two tests are geometrically disjoint in
what they compare; neither result can be re-derived from the other.

## What this does not settle

The excess is strongest where the background model is least trusted, and it
sits alongside two open problems that look related and may be the same thing:

- **Open item #10**, the low-energy antisymmetric residual at |b| = 45–60°
  (−6.4/+5.8/+4.5/−3.8σ at 1.5 GeV, dying by 12 GeV). Same energies, same
  latitudes, same order of magnitude.
- **Open item #8**, the structured template fitting *worse* at 1.5–2.6 GeV
  (Δ(−2lnL) = +243 at bin 0). Oscar has still not answered the energy-range
  question, and the low-energy end of this table is the same region of the
  data.

A monotonic fall from 14.7σ at 1.5 GeV to nothing above 35 GeV is the shape of a
foreground-modelling residual, not of a 20 GeV signal. **The honest reading is
that this is primarily a statement about template completeness, and only
secondarily a limit on a lobe.** Whether the 2.0σ that survives at 20.8 GeV is
the tail of the low-energy problem or something separate is not answered here,
and a spatial breakdown of where inside the footprint the excess sits is the
next thing to do.

## Provenance of the boundary

Verified 2026-09-08 before use, not taken on trust:

- The FITS `SRCSHA` matches, **byte for byte**, the JPEG embedded on page 28 of
  our own `docs/arXiv_2408.06312v2.pdf`. Oscar's "reconstructed from the
  original publication rather than tracing by eye" is literally true, and
  `erosita.py` checks the hash on every load.
- The solid cyan curve is the eROSITA edge (outer); the Fermi edge is the
  orange dashed curve inside it, drawn as PDF vector rather than raster.
- An **independent** re-extraction — own ellipse fit to the Mollweide frame, own
  projection inverse, own cyan mask — lands a median 0.51° from his boundary
  (p95 1.6°, about the printed stroke half-width). Longitude-flipped control:
  6.5°.
- Handedness: the lobes lean opposite ways, ⟨l⟩ = −9.46° north and +3.75° south.
  `erosita.check_handedness` asserts both on every load, and
  `tests/test_erosita.py` shows an l-flipped map fails it.
- Grid phase is exact — max |x − round(x)| = 0 — so the 0.125° transfer moves no
  flux and is not a resampling. The loader raises if that ever stops holding.
- `LATTEST = F`: no LAT data was used in building it.

## Reproducing

```bash
V=/Users/aba/work/gce-fisher/.venv/bin/python
$V scripts/boundary_test.py                        # both templates, ~2.5 min, 8.2 GB peak
$V scripts/boundary_test.py --theta-weight harmonic
$V scripts/boundary_test.py --exclude-bubble none
$V -m pytest tests/test_erosita.py -q              # 21 tests, ~1 s, no run dir needed
```

The delivery in `external/erosita_validation/` is **gitignored** (89 MB, the
FITS alone 80 MB), so a fresh clone cannot run this; `erosita.boundary_path`
raises with the provenance rather than falling back to anything.

## Refactor carried out for this

`patch_sums` moved out of `scripts/patch_test.py` into
`totani.patches.region_sums`, shared by both tests so they cannot drift in how
a template is built; regions are weight maps, and the boxes are the 0/1 case.
Verified by re-running `patch_test.py` and diffing against the committed
`patch_test_2026-09-06.json`: 667 numeric leaves, worst relative difference
3.4×10⁻¹² (float summation order), and the comparison detects a planted 0.1%
change. **No committed number moved.**
