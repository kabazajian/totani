# The halo search with Totani's own diffuse model — 2026-09-18

Kev (evening, 2026-09-18): "I want to use the same diffuse model as Totani."
`scripts/totani_diffuse.py` → `results/fits/totani_diffuse_2026-09-18.json`;
271 s; peak RSS 6.3 GB. **Second run**: the first (same evening) double-counted
ionized bremsstrahlung (`bremss_healpix` already contains it; 5% of brems,
~1% of gas) and lacked the ring variant; numbers below are from the corrected run
and differ from the first in the third digit.

## Which model, and where it came from

Totani (arXiv:2507.07209 §2.3) does **not** use the Ackermann+ 2017 "Sample
Model". He uses GALPROP v54 with the parameter set **SLZ6R30T150C2** from the
128-model grid of **Ackermann et al. 2012, ApJ 750, 3** (his ref. [47]) —
Lorimer pulsar source distribution, z_h = 6 kpc, R_h = 30 kpc, T_S = 150 K,
E(B−V) cut 2 mag — "the parameter set used for the ICS component of
gll_iem_v07". The 2017 paper's treatment (five Galactocentric gas rings,
three ICS components) is his §3.4.4 *variant* "gas rings & 3 ICS"; his five
parameter variants (SNR sources, z_h = 10, R_h = 20, T_S = 10⁵ K, cut 5 mag)
are neighbours in the same grid.

**The pre-computed GALPROP output of all 128 models is public**:
https://galprop.stanford.edu/PaperIISuppMaterial/ ships the galdef files
(`GALDEF.tar.gz`) and one archive per model
(`S{S,L,Y,O}_Z{4,6,8,10}_R{20,30}_T{150,100000}_C{2,5}.tar`). No WebRun
needed. `SL_Z6_R30_T150_C2.tar` (1.36 GB) is in `external/galprop_paperII/`
(gitignored, README). Each archive holds HEALPix NSIDE-128 skymaps at 38
energies (50 MeV–814 GeV) for π⁰ decay, bremsstrahlung (neutral and ionized
gas separately), ICS (total and its three ISRF components), and per-ring
H I / H₂ / H II splits (17 rings). Units verified as **dN/dE in
MeV⁻¹ cm⁻² s⁻¹ sr⁻¹** (E²·dN/dE at the GC peaks at 700 MeV — the pion bump).
`totani.galprop` integrates them onto the analysis bins exactly.

## What was run

Our pipeline — same data, cells, point sources (4FGL-DR4 modelled), Loop I
shells, isotropic — with `diffuse="galprop_totani"`: **one** `gas` template
(π⁰ + bremsstrahlung, ratio fixed) and **one** `ics` template, as Totani
§2.3 (7 templates per bin in place of 19); the 3-ICS split (9); and his
§3.4.4 **"gas rings & 3 ICS"** (13): gas in five Galactocentric groups —
GALPROP annuli 1 | 2–5 | 6–12 | 13 | 14–17 of Ackermann+ 2012 Table 2, i.e.
0–1.5 / 1.5–3.5 / 3.5–8 / 8–10 / 10–50 kpc — summed over H I, H₂, H II and
π⁰ + brems in spectrum space before clipping the small negative entries the
per-ring files carry (≤ 0.6% of a ring's flux; recorded), plus the three ICS
components. Each with Oscar's structured bubble map and with the flat mask. The A2020 model
re-scored alongside on the same run.

⚠ Not a reproduction of Totani's fit: his bubble template is his own (§3.1,
from the 4.3 GeV residual), his point-source handling and fitter are his.
This is *our pipeline with his diffuse model*.

## Results at 20.8 GeV

| diffuse model | bubble | n | halo | 68% | TS | quadrants [l+b+, l−b+, l+b−, l−b−] | Σ/joint |
|---|---|---|---|---|---|---|---|
| A2020 (19) | structured | 19 | 4,612 | [4,026, 5,237] | 23.5 | [783, -313, 1,446, 707] | 0.57 |
| A2020 (19) | flat | 19 | 7,705 | [7,306, 8,290] | 69.1 | [1,567, 730, 2,372, 1,836] | 0.84, all positive |
| **Totani GALPROP** | **structured** | 7 | **1,575** | [915, 2,219] | 5.8 | [-452, 447, 862, 884] | 1.11 |
| **Totani GALPROP** | **flat** | 7 | **5,006** | [4,473, 5,548] | 88.0 | [437, 953, 1,839, 1,455] | 0.94, all positive |
| Totani, 3 ICS | structured | 9 | 1,631 | [684, 2,572] | 3.2 | [-346, 1,192, 698, 830] | 1.45 |
| Totani, 3 ICS | flat | 9 | 5,512 | [4,657, 6,324] | 68.4 | [738, 1,653, 1,979, 1,333] | 1.03, all positive |
| **Totani, 5 gas rings + 3 ICS** | **structured** | 13 | **2,193** | [1,232, 3,104] | 3.5 | [-437, 760, 989, -354] | 0.44 |
| **Totani, 5 gas rings + 3 ICS** | **flat** | 13 | **5,906** | [5,093, 6,749] | 25.4 | [783, 1,341, 1,522, 1,238] | 0.83, all positive |

Injection slopes: a2020_structured 1.00–1.01; a2020_flat 1.00–1.01; totani_structured 0.99–1.18; totani_flat 0.97–1.17; totani_ics3_structured 0.99–1.18; totani_ics3_flat 0.97–1.17; totani_rings3ics_structured 0.68–1.07; totani_rings3ics_flat 0.67–1.21. The 1.17–1.18 are at 1.5 GeV in the GALPROP rows; **the 0.67–0.68 in the five-ring rows** is the inner gas rings going degenerate with the halo template at low energy — the gate fails there, so those bins are not measurements. 1.00 at 20.8 GeV in every row except the ring rows (1.0).

## The halo spectrum, per bin (counts)

| | 1.5 | 2.6 | 4.3 | 7.3 | 12.3 | 20.8 | 35 | 59 | 100 |
|---|---|---|---|---|---|---|---|---|---|
| A2020, structured | 106,612 | 56,683 | 31,305 | 13,639 | 8,697 | 4,612 | 267 | 185 | 71 |
| A2020, flat | 150,226 | 91,952 | 50,570 | 26,187 | 14,660 | 7,705 | 3,843 | 1,494 | 555 |
| **Totani, flat** | -50,883 | -7,558 | 5,418 | 8,456 | 7,206 | 5,006 | 2,188 | 1,054 | 299 |
| Totani, structured | -94,497 | -39,048 | -16,255 | -5,456 | 621 | 1,575 | -214 | -220 | -204 |
| Totani 3-ICS, flat | -21,807 | 11,337 | 9,632 | 8,183 | 7,827 | 5,512 | 3,149 | 892 | 370 |
| Totani rings+3-ICS, flat | -19,384 | 5,038 | 16,255 | 10,911 | 9,069 | 5,906 | 2,722 | 1,030 | 176 |
| Totani rings+3-ICS, structured | -29,076 | -8,817 | -15,051 | -4,092 | 2,112 | 2,193 | -1,121 | -347 | -266 |

🚨 **With his diffuse model and the flat mask, Totani's headline appears in
our pipeline**: the halo is *negative* at 1.5–2.6 GeV (his paper: "consistent
with zero below 2 GeV", and §3.4.4 speaks of "the negative halo fits at the
lowest energies in the baseline analysis"), positive from 4 GeV, and in
E × counts it **peaks at ~20 GeV** (23k, 62k, 89k, 104k, 77k, 62k, 30k from 4.3 to 100 GeV), with TS 88 and the four quadrants all positive and
consistent to 6%. That is his 20 GeV peak, his quadrant consistency, his
low-energy sign. Under A2020 with the flat mask the spectrum has no peak
(E × counts declines monotonically from 225k at 1.5 GeV).

🚨 **Oscar's structured bubble map removes most of it under his model too**:
1,575 [915, 2,219], TS 5.8, and the quadrant that was −313 under A2020 is −452
here. The bubble template, not the diffuse model, is what the 20.8 GeV
number is most sensitive to: flat → structured is −40% under A2020 and −69%
under GALPROP.

## The price of the two-template model

Background-only −2lnL relative to A2020+structured, per bin
(1.5 / 2.6 / 4.3 / 7.3 / 12.3 / 20.8 / 35 / 59 / 100 GeV):

| model | Δ(−2lnL) |
|---|---|
| A2020, flat | -243 / -45 / +54 / +100 / -48 / +10 / +54 / +26 / -7 |
| GALPROP-7, flat | +1,162 / +415 / +164 / +155 / +33 / +54 / +87 / +32 / +1 |
| GALPROP-7, structured | +2,308 / +1,070 / +400 / +193 / +72 / +17 / +26 / +12 / +6 |
| GALPROP-9 (3 ICS), flat | +902 / +366 / +154 / +126 / +5 / +33 / +81 / +29 / -0 |
| GALPROP-13 (5 rings + 3 ICS), flat | +694 / +246 / +99 / +91 / -66 / -19 / +48 / +18 / -5 |
| GALPROP-13 (5 rings + 3 ICS), structured | +1,520 / +722 / +274 / +135 / +24 / +1 / +19 / +8 / +4 |

Mean |residual| per cell over the halo region: 5.93 vs 4.71 at 1.5 GeV, 2.87
vs 2.39 at 4.3, 1.27 vs 1.21 at 20.8 (GALPROP-7 vs A2020-19). The seven-template model is a worse description at every energy,
as twelve fewer free parameters must be; the halo template then absorbs
diffuse residuals of whichever sign the model leaves — negative at low
energy under GALPROP, positive under A2020. **Neither low-energy halo is
physical; both are diffuse-model residual.** What is robust across the two
models is: a positive 20.8 GeV amplitude of 4,600–5,000 with the flat mask
or A2020+structured, and its collapse to ~1,600 with GALPROP+structured.

## What this answers, and what it opens

- **"Can we get Totani's model?" Yes, already have it**: it is the public
  Ackermann+ 2012 output, no GALPROP run needed. His five parameter variants
  are 1.4 GB each, same place (downloading; `totani_galprop_variants.py`).
- **His "gas rings & 3 ICS" variant, built**: with the flat mask the 20.8 GeV
  halo goes 5,006 → 5,906 [5,093, 6,749] and TS 88 → 25; his paper says the
  halo "decreases in bins below 12 GeV but remains almost unchanged at the
  21 GeV bin" — here 4.3–12 GeV go *up* (5,418 → 16,255 at 4.3) and 20.8 GeV
  up 18%. The 13-template model fits far better (Δ(−2lnL) = +694 at 1.5 GeV vs
  +1,162 for the 7-template one, relative to A2020), but the injection
  gate fails at 0.67 in a low-energy bin: the inner gas rings are degenerate
  with the halo there. With the structured map: 2,193 [1,232, 3,104], TS 3.5,
  two quadrants negative. Same conclusion as the two-template version: the
  bubble template decides.
- **The 20 GeV peak is reproducible with his ingredients** and is not
  reproduced with the A2020 ring model. It is a statement about the
  GALPROP two-template diffuse model's residual spectrum — over-prediction at
  1–3 GeV, under-prediction at 7–35 GeV at |b| ≥ 10° — combined with a flat
  bubble mask.
- **The structured bubble map is decisive under both models.** That is the
  thing to settle with Oscar (his energy-range answer, the eROSITA brightness
  map, and the Loop I interaction all bear on it).
- The quadrant test passing under GALPROP+flat (0.94) and failing under
  A2020+structured (0.57) says the quadrant failure is about the *shape
  disagreement between the A2020 ring templates, the structured map and
  Loop I* in the l < 0, b > 0 quadrant — new evidence for open item 7.

Analysis and text: Claude Fable 5.1 (`claude-fable-5-1[1m]`), under Kev's
direction; 158 tests passing.
