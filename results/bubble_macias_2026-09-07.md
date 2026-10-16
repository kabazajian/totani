# Oscar's structured bubble template: what it does to the halo, and what it does not fix

**Delivered:** 2026-09-06, `external/FermiBubbles_structured_fig12_inpainted_ALLSKY_Macias.fits`
— Fig. 12 of **Ackermann+ 2017 (arXiv:1704.03910)** — *not* the 2014 bubbles
paper, as this note originally said; corrected 2026-09-18 on Oscar's
email; derived from 1–10 GeV data — (the *structured* template), point-source holes
inpainted, all-sky CAR 0.5°. **Data:** `runs/fff2de5bc7e0`, disk-excluded.
**Script:** `scripts/bubble_macias.py` (one stack build, then cell-space
variants; cells cached in `results/fits/bubble_macias_cells_2026-09-07.npz`).
**Output:** `results/fits/bubble_macias_2026-09-07.json`. Loader:
`templates.bubble_macias_on_target`, on the flux-conserving `regrid_car` path
— the flat-mask loader refuses it by design. Stack switch `bubble="macias"`;
the default stays `"ackermann"` so every earlier result reproduces.

## Tier 1 (passes)

Full ROI coverage, no NaN, no holes (the `INPAINT` HDU flags filled pixels:
759 inside the inner 40°, 5,907 all-sky). Handedness consistent with the
Ackermann mask — pixel correlation 0.375 as-is vs 0.329 l-flipped, same sign of
l-asymmetry in both hemispheres. 94% of the flat footprint is non-zero in
Oscar's map; **51% of his map's flux lies outside that footprint** (33% of ROI
pixels non-zero vs 19%). N/S flux share 0.444/0.556, essentially the flat
mask's area share 0.449/0.551. Flux through the regrid conserved to < 3%
(`tests/test_bubble_macias.py`).

## The halo amplitude, per bin

| bin | E (GeV) | Ackermann flat | Oscar structured | ratio | TS flat → Oscar |
|---|---|---|---|---|---|
| 0 | 1.51 | 150,226 | 106,612 | 0.71 | 390 → 217 |
| 1 | 2.55 | 91,952 | 56,683 | 0.62 | 317 → 130 |
| 2 | 4.31 | 50,570 | 31,305 | 0.62 | 201 → 70 |
| 3 | 7.28 | 26,187 | 13,639 | 0.52 | 142 → 39 |
| 4 | 12.29 | 14,660 | 8,697 | 0.59 | 61 → 21 |
| **5** | **20.76** | **7,705** [7,306, 8,290] | **4,612** [4,026, 5,237] | **0.60** | **69 → 24** |
| 6 | 35.06 | 3,843 | 240 [−11, 29]×35 | 0.06 | 20.5 → 0.2 |
| 7 | 59.22 | 1,494 | 232 | 0.16 | 19.3 → 1.0 |
| 8 | 100.02 | 555 | 71 | 0.13 | 5.2 → 0.1 |

**Oscar's map cuts the 20.8 GeV halo by 40%**, far outside the ±6% statistical
intervals, and **removes it entirely above 35 GeV** (TS ≤ 1). The SED proxy
(counts × E) goes from "flat 2.5–60 GeV" to: declining from 1.5 GeV, flat at
~100 from 7 to 21 GeV, then a cliff to 8 ± 20 at 35 GeV. **Still no 20 GeV
peak** — bins 3, 4, 5 are 99, 107, 96 with overlapping intervals — but now a
cutoff above 21 GeV. Injection slopes 0.997–1.013.

## Is it a better background, or a more halo-shaped one? Both, by energy

Background-only −2lnL, same data, same 19 parameters (Oscar − Ackermann):

| bin | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | sum |
|---|---|---|---|---|---|---|---|---|---|---|
| Δ(−2lnL) | **+243** | +45 | −54 | **−100** | +48 | −10 | **−54** | −26 | +7 | +98 |

Neither template is uniformly better. Oscar's is a *worse* description at
1.5–2.6 GeV (dominated by bin 0), *better* at 4–7 GeV and 35–60 GeV, and
marginally better at 20.8 GeV (−10). Where it is better at high energy it also
removes the halo completely — **above 35 GeV the "halo" was bubble structure
the flat mask could not represent.** At 20.8 GeV the fit improves slightly and
the halo drops 40% with TS 69 → 24.

The map is much more halo-like than the flat mask: correlation with the NFW-ρ²
template in the halo cells **0.52 → 0.74**; halo leverage against the full
basis 0.266 → 0.211 (full ROI), 0.183 → 0.164 (band). Bubble fitted fraction
at 20.8 GeV 9.0% → 12.5%.

## What it does not fix

**Quadrants get worse.** At 20.8 GeV the four quadrant amplitudes go from
[1,567, 730, 2,372, 1,836] (sum/baseline 0.84) to **[753, −313, 1,446, 707]**
(0.56, one negative). At 1.5 and 4.3 GeV the quadrants sum to **−0.18 and
−0.10 of the baseline** — the opposite sign. With Oscar's map the low-energy
"halo" is a common amplitude fitted to quadrants that disagree even in sign.

**The N/S problem is untouched.** Splitting his map into independent
north/south normalizations improves the background-only fit by Δ(−2lnL) = −9.0
for one parameter, the same as splitting the flat mask (−9.7); the matched-θ
on/off asymmetry at 20.8 GeV with his single-normalization map is +2.4σ/+1.8σ
in the south (flat mask: +1.5/+1.2), and only the split brings it to ≤ 1.2σ.
His N/S flux balance equals the flat mask's area balance, so this was
expected: the within-lobe structure does not carry the N/S brightness
difference the data want.

Mean |σ| per cell with the halo included is slightly *worse* with his map at
bins 0–5 (4.34 → 4.58 at 1.5 GeV; 1.09 → 1.16 at 20.8 GeV).

## Consequences

1. **The 20.8 GeV amplitude is bracketed 4,612–15,242 by the bubble treatment
   alone** (3.3×) — no bubble, flat, structured. Combined with the existing
   variants the envelope is unchanged at 5.2× (2,931 from the 1° PS mask to
   15,242), but the low end now has a second, independent, physically
   motivated member. The bubble template is the single most consequential
   diffuse choice in this ROI after the point-source treatment.
2. **The Galp21 leverage forecast is unaffected** (`scripts/band_forecast_bubble.py`):
   self-consistent Galp21 leverage in the band 0.050–0.089 (1.79×) with
   Oscar's map vs 0.048–0.088 (1.82×) with the flat mask; Pohl/Galp21 ratio
   1.8–3.3× vs 2.1–3.8×. The Youyou note's conclusions stand; its Pohl
   baseline row moves 0.183 → 0.164.
3. **Which template is "right" is energy-dependent and undecided.** The honest
   procedure is to carry both as members of the systematics envelope until the
   eROSITA boundary allows the on/off test that would discriminate them, and
   to give either one independent N/S normalizations, which the data demand
   at ~3σ regardless.
4. The high-energy result is clean and new: with a structured bubble template,
   **there is no halo-shaped excess above 35 GeV** (TS ≤ 1 in three bins).
