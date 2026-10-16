# Loop I from the Haslam 408 MHz map — 2026-09-18

Kev's decision (evening, 2026-09-18): an empirical alternative to the two
geometric Wolleben shells, as a variant. `scripts/loopi_haslam_variant.py` →
`results/fits/loopi_haslam_variant_2026-09-18.json`; 280 s; peak RSS 6.7 GB.
Structured (`macias`) baseline throughout.

## The template

Remazeilles et al. 2015 destriped/desourced Haslam 408 MHz (LAMBDA; sha256
in `external/haslam/README.md`), HEALPix NSIDE 512, sampled nearest-pixel
onto the 0.125° grid by a hand-written RING pixelisation (`totani.healpix`;
no healpy here; exact round trip over every pixel asserted in
`tests/test_haslam.py`). Two shapes:

- **`haslam_nps`** — the map inside the disc subtended by Wolleben shell 2
  (centre (l, b) = (332°, 37°), d = 95 pc, r_out = 82 pc → 59.7° radius),
  zero outside. ⚠ Not "the shells' footprint": shell 1 has the Sun *inside*
  its outer sphere, so its path length is non-zero along every sight line
  and its footprint is the whole sky (found by the tests). Shell 2's disc is
  the geometric model's own statement of where Loop I's bright northern
  structure is.
- **`haslam_raw`** — the map over the whole ROI: a synchrotron template
  everywhere, competing with the IC rings, not only at Loop I.

Both with one free normalisation per bin (the shells have two). 408 MHz
brightness is electrons × B-field; the gamma rays are electrons × ISRF —
related, not identical, which is why the LAT team and Totani used geometry.

## Results

−2lnL, background-only, relative to the two shells (negative = better):

| Loop I model | params | 1.5 GeV | 2.6 | 4.3 | 7.3 | 12.3 | 20.8 |
|---|---|---|---|---|---|---|---|
| two shells (baseline) | 2 | 0 | 0 | 0 | 0 | 0 | 0 |
| `haslam_nps` | 1 | **+167** | +130 | −10 | −18 | +23 | +5 |
| `haslam_raw` | 1 | **−283** | −113 | −47 | +14 | −5 | −7 |
| shells + `haslam_nps` | 3 | −110 | −58 | −60 | −37 | 0 | −1 |
| none | 0 | +2,986 | +1,586 | +664 | +224 | +75 | +23 |

Halo, quadrants, and the structured-vs-flat bubble comparison:

| Loop I model | halo @ 20.8 | 68% | TS | quadrant Σ/joint | halo @ 1.5 | struct − flat @ 1.5 / 2.6 |
|---|---|---|---|---|---|---|
| two shells | **4,612** | [4,026, 5,237] | 23.5 | 0.57, one negative | 106,612 | +243 / +45 |
| `haslam_nps` | 4,275 (−7%) | [3,698, 5,020] | 21.9 | 0.57, one negative | 99,343 | +140 / −31 |
| `haslam_raw` | **3,656 (−21%)** | [2,814, 4,466] | 15.3 | **0.25**, one negative | 41,191 | −52 / −141 |
| shells + nps | 4,556 (−1%) | [3,954, 5,180] | 23.2 | 0.55, one negative | 106,612 | +222 / +31 |
| none | 4,612 | [4,127, 5,332] | 26.8 | 0.54, one negative | 130,842 | −215 / −258 |

Injection slopes 0.97–1.04 everywhere; 1.00 at 20.8 GeV in every row.

## Reading

1. **As a Loop I model, the empirical NPS shape is not better than the
   shells.** With one parameter it loses 167 and 130 units at 1.5–2.6 GeV
   and is level above 4 GeV. Added beside the shells it carries something
   they miss (−110 at 1.5 GeV for one parameter), and the halo barely moves
   (−1.2%). So the shells are a competitive Loop I model; the empirical shape
   is a modest refinement, not a replacement.
2. **The raw synchrotron map is a different animal.** It beats the two shells
   with one parameter at every energy below 4 GeV, removes the structured
   map's low-energy disagreement with Loop I (+243 → −52), and moves the halo
   −21% — but it is not modelling Loop I; it is a whole-sky synchrotron
   template absorbing IC-like structure across the ROI, and the quadrant test
   collapses to Σ = 0.25. That −21% is a real systematic in how the diffuse
   electron emission is modelled at |b| ≥ 10°, and it belongs in the envelope
   at that level of description, not as a Loop I variant.
3. **The quadrant test fails in every row.** Loop I is not what the l < 0,
   b > 0 quadrant is disagreeing about either.
4. The structured-vs-flat interaction of the morning (+243) is a
   three-way interaction between the structured map, Loop I, and whatever
   `haslam_raw` absorbs across the ROI. It is removed only by a template that
   covers the whole sky, which points at the diffuse model rather than at
   Loop I specifically — which is what `totani_diffuse_2026-09-18` tests.

## Bottom line

- The two Wolleben shells stay. `shells + haslam_nps` is a defensible
  refinement (Δ(−2lnL) = −110 for one parameter; halo −1%).
- Loop I choice contributes 4,275–4,612 to the 20.8 GeV envelope; the raw
  synchrotron template's 3,656 is a diffuse-model statement.
- Nothing here touches the quadrant test.

Analysis and text: Claude Fable 5.1 (`claude-fable-5-1[1m]`), under Kev's
direction.
