# Oscar's matched-angular-distance test, on the two patches Zhang+ actually specify

**Data:** `runs/fff2de5bc7e0`, disk-excluded background fit (Totani §3.2).
**Script:** `scripts/patch_test.py`. **Output:** `results/fits/patch_test_2026-09-06.json`.
**Prompted by:** Oscar Macias, 2026-09-05 — "compare regions along these
radio/X-ray features with regions away from them, matched in angular distance
from the Galactic center. A spherical halo should give the same intensity at
the same angular distance. A lobe or shell need not."

## What could be digitized, and what could not

Zhang et al. (arXiv:2408.06312) specify **two** regions numerically, both in
their §3.1, and both inside Totani's ROI:

| patch | definition | what it is |
|---|---|---|
| R1 | 20° < l < 40°, 45° < b < 60° | northern "outer region" |
| R3 | 20° < l < 40°, −45° < b < −30° | southern "outer region" |

"Outer region" = outside the Fermi Bubbles, inside the eROSITA Bubbles.

**Not specified numerically, and therefore not usable:** the comparison patches
R2/R4 (longitude ranges never given), the background R0 ("northeastern high
galactic latitude sector"), the Fermi cap R5, **and both bubble boundary
curves** — those exist only as figure elements drawn from Predehl+2020. The
data-availability statement lists survey URLs only; no region files.

## The trick that makes the test runnable anyway

Reflecting a patch in `l → −l` or `b → −b` leaves θ_GC **exactly** invariant.
So each Zhang box and its reflections form a matched-θ set that needs no
boundary, no digitized curve, and no assumption about where the eROSITA edge
runs. `tests/test_patch_test.py` asserts the invariance rather than assuming it.

Two things had to be handled that the naive version gets wrong:

1. **The raw sky is nowhere near symmetric.** Raw counts between reflected
   patches differ by up to 52σ (R3: 25,364 vs 38,535 in its b-mirror). Oscar's
   null is a null only *after* the asymmetric foreground is modeled; the test
   is on background-model residuals, not on intensities.
2. **A spherical halo does not predict equal counts in reflected patches.**
   Exposure and the extended-source mask break the symmetry θ_GC preserves — by
   27% between the R3 pair (J = 1.56% vs 1.99% of the ROI). At 20.8 GeV that is
   a −33 count expected difference against a ±37 statistical error, i.e. fully
   comparable. It is subtracted, leaving `lobe` = the part no spherical
   component can produce.

⚠ The patches do **not** tile the 10° cell grid — their latitudes cut cells in
half — so this runs on pixel sums. Validation: the pixel-sum path reproduces
the cell-space stack to 7×10⁻¹⁵ on an exact-cell box, with identical counts.

## Result at 20.8 GeV — Totani's energy

| test | lobe (counts) | σ | 95% UL / halo in the same patch |
|---|---|---|---|
| R1 vs l-mirror | +28 ± 31 | +0.90 | < 1.29 |
| R3 vs l-mirror | −1 ± 37 | −0.03 | **< 0.61** |
| R1 vs b-mirror | +5 ± 30 | +0.15 | < 0.91 |
| R3 vs b-mirror | +31 ± 41 | +0.76 | < 0.95 |

**Every reflection test is consistent with zero at 20.8 GeV, and at every
energy above ~12 GeV.** The best constraint is the southern patch east–west:
any lobe-like component in R3 is **< 61% of the halo signal attributed to that
same patch**. Not vacuous, not decisive.

## Where there *is* structure — and why it is probably not the lobes

| bin | E (GeV) | R1 vs l-mir | R1 vs b-mir | R3 vs l-mir | R3 vs b-mir |
|---|---|---|---|---|---|
| 0 | 1.5 | **−9.4σ** | **−8.8σ** | +1.4σ | −2.5σ |
| 1 | 2.6 | **−5.8σ** | **−5.2σ** | +2.1σ | +0.0σ |
| 2 | 4.3 | **−5.0σ** | **−4.9σ** | +2.1σ | +0.8σ |
| 3 | 7.3 | **−4.5σ** | **−4.3σ** | +2.7σ | +1.6σ |
| 4 | 12.3 | −1.5σ | −0.7σ | +1.8σ | +1.4σ |
| 5 | 20.8 | +0.9σ | +0.2σ | −0.0σ | +0.8σ |

Everything significant is (a) below 12 GeV and (b) on **R1**, in both
reflections, with a **negative** sign — R1 is *fainter* than its matched-θ
reflections, not brighter. That is the opposite of extra eROSITA-associated
emission sitting in R1.

### ⚠ Correcting an attribution made before it was checked

An earlier version of this note said the cause was "likely Loop I, since R1
sits on the North Polar Spur." **That does not survive the check and is
withdrawn.** Loop I template content per patch, against the bin-0 residual:

| patch | l | b | shell 1 | shell 2 | σ (bin 0) |
|---|---|---|---|---|---|
| R1 | 20…40 | 45…60 | 75.2 | 99.1 | **−6.36** |
| R1_lmir | −40…−20 | 45…60 | 75.3 | 52.2 | **+5.78** |
| R1_bmir | 20…40 | −60…−45 | 63.6 | 0 | +4.49 |
| R1_lbmir | −40…−20 | −60…−45 | 85.3 | 0 | −3.80 |
| R3 | 20…40 | −45…−30 | 81.0 | 0 | +0.06 |
| R3_lmir | −40…−20 | −45…−30 | 54.1 | 0 | +1.35 |
| R3_bmir | 20…40 | 30…45 | 84.6 | 90.4 | **+4.93** |
| R3_lbmir | −40…−20 | 30…45 | 47.1 | 48.6 | +1.02 |

**R1 and R3_bmir have nearly identical Loop I content (174 vs 175 total, 99 vs
90 in shell 2) and residuals of opposite sign.** Loop I brightness does not
track the residual, per shell or in total. The Spur story was plausible and
wrong.

### What the pattern actually is

The significant low-energy residuals are confined to the **|b| = 45–60° band**,
where they change sign under *both* reflections — (+l, N) −6.4σ, (−l, N) +5.8σ,
(+l, S) +4.5σ, (−l, S) −3.8σ. That is antisymmetric under point reflection
through the Galactic Centre. In the |b| = 30–45° band everything is positive
and mild (+0.1 to +4.9σ).

**No component we have accounts for it, and this note does not claim one.**
The unexciting reading is the one the project already has evidence for: at
1.5 GeV the background model is simply inadequate here — mean |σ| per cell is
4.34 (`halo_search_v3`), and the PS-masked variant drops the bin-0 halo by a
factor of 38 (`systematics`). The effect dying by 12 GeV, as the residuals
approach Poisson, fits that reading. It is worth a look, but it is a
low-energy background problem, not a lobe.

R3 is clean at **every** energy (max 2.7σ).

## What this does and does not settle

- It does **not** test Zhang's actual claim. Their signal is defined relative to
  the eROSITA boundary; a reflection test asks a strictly weaker question
  ("is the sky lopsided at matched θ_GC?"). Both of a pair can sit inside the
  same bubble and cancel. **This is the main limitation** — Zhang's own northern
  cuts put the eROSITA enhancement across −70° ≲ l ≲ +40° at b ≈ +65–70°, which
  is wider than our entire |l| ≤ 60° ROI, so R1's l-mirror is plausibly inside
  the northern bubble too. The southern result is the more meaningful one.
- It does establish that at 20.8 GeV there is **no lopsidedness at all** to
  attribute to anything, at the ±30–40 count level per patch.
- Combined with the N/S bubble-normalization result (2026-09-05, uncommitted):
  the anisotropies our background leaves behind are a **low-energy**
  phenomenon, and the halo amplitude at 20.8 GeV moves by < 1% when they are
  absorbed. The lobe question and the factor-of-5 halo systematic are close to
  orthogonal.

**Next:** the real version needs the eROSITA boundary curve, requested from
Oscar 2026-09-06.
