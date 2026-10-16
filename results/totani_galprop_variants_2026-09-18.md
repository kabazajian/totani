# Totani's five GALPROP parameter variants — 2026-09-18

His §3.4.4 / Figure 15 (bottom): five models one step from SLZ6R30T150C2 on
the Ackermann+ 2012 grid. Archives from galprop.stanford.edu/PaperIISuppMaterial
(1.2–1.5 GB each, gitignored), the three total maps each needs extracted;
`scripts/totani_galprop_variants.py` → `fits/totani_galprop_variants_2026-09-18.json`;
230 s; peak RSS 6.2 GB. Each fitted in his baseline configuration (gas + ICS,
two free templates per bin) with both bubble templates.

## How different the models are in our ROI (|l| ≤ 60°, 10° ≤ |b| ≤ 60°, 1.5 GeV)

Ratio to the baseline, median over ROI pixels (min–max):

| variant | π⁰ | bremss | ICS |
|---|---|---|---|
| SNR sources (`SS_Z6_R30_T150_C2`) | 0.97 (0.90–1.00) | 0.93 | **0.83 (0.63–0.93)** |
| z_h = 10 kpc (`SL_Z10_…`) | 1.02 | 0.91 | **1.23 (1.05–1.28)** |
| R_h = 20 kpc (`SL_Z6_R20_…`) | 1.002 | 0.998 | 0.999 |
| T_S = 10⁵ K (`…_T100000_C2`) | 0.97 (0.83–1.54) | 0.97 | 1.03 |
| E(B−V) cut 5 mag (`…_T150_C5`) | 0.99 (0.74–2.87) | 0.98 | 0.98 |

R_h = 20 changes the ROI by 0.2% — its effect is in the outer Galaxy — so
its fit is identical to the baseline to the fitter's resolution (that was
checked, not assumed). T_S and the E(B−V) cut change *shapes* (pixel ratios
0.7–2.9) but not levels; SNR sources and z_h = 10 change the ICS level and
morphology.

## Results at 20.8 GeV

| variant | flat mask | | | structured map | | |
|---|---|---|---|---|---|---|
| | halo | TS | quadrants | halo | TS | quadrants |
| baseline SLZ6R30T150C2 | 5,006 [4,473, 5,548] | 88 | 0.94, all + | 1,575 [915, 2,219] | 5.8 | one − |
| **SNR sources** | **6,356 [5,886, 6,857]** | **177** | 0.93, all + | **3,037 [2,472, 3,655]** | **27** | one − |
| z_h = 10 kpc | 4,668 [4,140, 5,251] | 73 | 0.94, all + | 1,237 [565, 1,903] | 3.4 | one − |
| R_h = 20 kpc | 5,006 [4,475, 5,550] | 88 | 0.94, all + | 1,575 [917, 2,221] | 5.8 | one − |
| T_S = 10⁵ K | 5,062 [4,510, 5,597] | 88 | 0.94, all + | 1,631 [958, 2,277] | 6.0 | one − |
| E(B−V) cut 5 mag | 4,949 [4,427, 5,499] | 87 | 0.90, all + | 1,462 [833, 2,135] | 5.2 | one − |

Injection slopes 0.97–1.18 (1.17–1.18 at 12.3 GeV in every GALPROP row —
a low-level degeneracy flag for that bin; 1.00–1.04 at 20.8 GeV).

## The halo spectrum (counts), flat mask

| | 1.5 | 2.6 | 4.3 | 7.3 | 12.3 | 20.8 | 35 | 59 | 100 |
|---|---|---|---|---|---|---|---|---|---|
| baseline | −50,883 | −7,558 | 5,418 | 8,456 | 7,206 | 5,006 | 2,188 | 1,054 | 299 |
| **SNR** | **+26,653** | **+34,010** | 25,285 | 16,639 | 10,560 | 6,356 | 2,856 | 1,146 | 422 |
| z_h = 10 | −72,690 | −17,635 | 602 | 6,274 | 6,460 | 4,668 | 2,055 | 926 | 280 |
| T_S = 10⁵ | −50,883 | −6,298 | 6,020 | 8,456 | 7,330 | 5,062 | 2,188 | 915 | 304 |
| cut 5 mag | −55,729 | −8,817 | 4,816 | 7,911 | 7,082 | 4,949 | 2,162 | 1,042 | 299 |

## Reading, against his §3.4.4

Totani: *"In most cases, the fit results are almost unchanged from the
baseline analysis, but relatively large changes are observed when the
cosmic-ray source distribution is changed to SNR. In this case, the halo
component is positively fitted at all energy bins and is nearly zero at the
lowest energy bin. This suggests that one reason for the negative halo fits
at the lowest energies in the baseline analysis may be the uncertainty of
the gas component of GALPROP."*

**Reproduced, feature for feature, with the flat mask:** four of the five
variants move the 20.8 GeV halo by ≤ 7% and leave the spectrum's shape
alone; SNR sources are the exception — the halo becomes positive in every
bin, the lowest bin goes from −50,883 to +26,653 (small next to the A2020
model's +150,226, i.e. "nearly zero" on his scale), TS at 20.8 GeV doubles to
177, and the 20.8 GeV halo rises 27%. The E×counts spectrum turns from a
peak at 20 GeV (23k, 62k, 89k, **104k**, 77k) into a broad plateau
(109k, 121k, 130k, **132k**, 100k from 4.3 to 35 GeV). His attribution — the
low-energy sign is the GALPROP gas/ICS residual — is what our ICS-level
comparison above says: the SNR model's ICS is 17% lower in the ROI, and the
halo template takes up the difference at every energy.

**And the structured map does to every variant what it does to the
baseline**: 1,237–1,631 for four of them (TS 3–6), 3,037 (TS 27) for SNR; one
negative quadrant in all six. The bubble template is the dominant choice
under every GALPROP parameter set he considered.

## Bottom line

- His GALPROP-parameter systematics are reproduced: small except for the
  source distribution, which is a 27% (flat) / 93% (structured) effect at
  20.8 GeV and flips the sign of the low-energy halo.
- The GALPROP-parameter envelope at 20.8 GeV in his configuration:
  4,668–6,356 (flat), 1,237–3,037 (structured).
- The structured-vs-flat ratio is 0.26–0.48 across all six model sets. No
  GALPROP parameter choice closes that gap.

Analysis and text: Claude Fable 5.1 (`claude-fable-5-1[1m]`), under Kev's
direction; 162 tests passing.
