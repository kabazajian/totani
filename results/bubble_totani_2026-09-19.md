# Totani's bubble template, built as he did — 2026-09-19

Kev (2026-09-19): "Let's build that bubble template as he did." `totani.bubble_totani`
+ `scripts/bubble_totani.py` → `fits/bubble_totani_2026-09-19.{fits,json}`,
figure `bubble_totani_2026-09-19.png`. 256 s; peak RSS 7.3 GB.

## The recipe (his §3.1–3.2) and what was done

1. **Construction fit, disk included.** At 1.5 and 4.3 GeV, the sharp-edged
   flat mask (Ackermann+ 2014) is fitted over all 144 cells together with his
   diffuse model (GALPROP SLZ6R30T150C2 gas + ICS), the two Loop I shells,
   isotropic, 4FGL-DR4 point sources, and **the GC GeV excess as a known
   component** — NFW-ρ^2.5 with his NFW parameters (his §2.3), a new
   line-of-sight integral `nfw_rho_power_on_target`. Fitted amplitudes at
   4.3 GeV: gas 1.35, ICS 1.32, point sources 1.08 (relative to the model),
   GCE 7.4×10⁻³⁵.
2. **Bubble image** = flat template + residual, per pixel: counts −
   Σ_{l≠bubble} f_l T_l, divided by exposure × solid angle, smoothed with a
   σ = 1° Gaussian on the sky (σ/cos b along l on the plate-carrée grid;
   masked pixels neither contribute nor dilute). Point sources are deposited
   at their pixel without PSF — at 4.3 GeV the 68% PSF is 0.2° against 1°
   smoothing — so bright sources leave ring/centre residuals, as his
   figure 1 also shows.
3. **Split** the 4.3 GeV image into positive and negative parts: two
   energy-independent templates that replace any bubble template. His sign
   rule (§3.2): positive template f ≥ 0, negative template free. Our fitter
   keeps background amplitudes ≥ 0, so each bin's background-only fit is
   done with the negative template and with its inverse and the better
   likelihood kept; **the inverted template wins in every bin below 100
   GeV** — he anticipated exactly that ("the inverted negative template may
   be favored… if the discrepancy between GALPROP and the data decreases in
   higher energy bins").

Not reproduced: his iteration of the flat boundary by eye, his MCMC, his
PSF-convolved sources.

## The images

The 4.3 GeV image shows the two lobes (south strong and sharp-edged, north
present), and 90% of the flat-mask footprint is positive. Both images show
his headline feature: **a band of negative residual in the north at
b ≈ 10–20°, strongest at 1.5 GeV** — "regions of particularly strong flux in
the gas component map of GALPROP… likely due to discrepancies between the
GALPROP model and the data" (his §3.1). That band is the negative template.
Positive covers 57% of the grid, negative 41%.

## The halo search with his templates (NFW-ρ², |b| ≥ 10°)

| configuration | halo @ 20.8 | 68% | TS | quadrants | Σ/joint |
|---|---|---|---|---|---|
| **his diffuse model + pos/neg (his configuration end to end)** | **6,412** | [5,904, 6,898] | **172** | [1,100, 1,951, 1,243, 1,673] | 0.93 |
| his diffuse model + positive only | 6,637 | | 187 | | 0.87 |
| A2020 diffuse + his pos/neg | 8,437 | [7,959, 8,935] | 60 | [618, 3,188, 1,573, 1,768] | 0.85 |
| A2020 diffuse + positive only | 8,661 | | 77 | | 0.90 |

Halo per bin, his configuration (counts): −58,152 / −6,298 / 7,826 / 12,002 /
9,815 / **6,412** / 3,176 / 1,447 / 579. E × counts: −88, −16, 34, 87, 121,
**133**, 111, 86, 58 (×10³) — a peak at 20 GeV, negative below 3 GeV.
Injection slopes 0.99–1.00 in every bin.

**The consistency check** (his §3.2): the residual templates absorb
everything at 4.3 GeV except the GC excess, whose counts were subtracted in
the construction. In his ρ^2.5 fit the GCE recovers its construction value
at 4.3 GeV; in the ρ² fits the halo template takes that role. The GCE's
|b| ≥ 10° counts from the construction fit are **4,914** at 4.3 GeV; the ρ²
halo there is **7,826 [6,002, 9,339]** — the same order, 1.6× (the flatter
profile picks up some residual as well). At 1.5 GeV the GCE's 9,256 is
swamped by the −58,152 the negative low-energy residual produces. Not zero,
and not expected to be.

## Reading

- **His result reproduces in his configuration, strongly**: 6,412 ± 500 at
  20.8 GeV, TS 172, quadrants consistent to 7%, the 20 GeV peak, the
  negative low-energy bins. Every ingredient is his (diffuse model, bubble
  construction, Loop I geometry, template count) except the fitter and the
  flat boundary.
- **It is a stronger detection than with the flat mask alone** (5,006,
  TS 88) because the negative template, inverted, is absorbing the GALPROP
  low-energy over-prediction where it sits (north, b = 10–20°), which the
  flat mask could not.
- **Under our diffuse model his bubble construction gives 8,437, TS 60** —
  larger than A2020 + flat (7,705) and far from A2020 + structured (4,612).
- What this analysis has shown across the day: the 20 GeV halo amplitude
  is set by the bubble template — flat 5,006 / his residual 6,412 / Oscar's
  structured 1,575 under his diffuse model; 7,705 / 8,437 / 4,612 under
  A2020. His residual templates are, by construction, "whatever the model
  leaves at 4.3 GeV, positive and negative"; the halo search then measures
  the *spectral* departure from 4.3 GeV. Oscar's map is an independent
  morphology derived at 1–10 GeV from a different residual. They disagree
  by a factor 4 in what they leave for a halo at 20 GeV, and that
  disagreement is now the whole question.

Analysis and text: Claude Fable 5.1 (`claude-fable-5-1[1m]`), under Kev's
direction; 165 tests passing.
