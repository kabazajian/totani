# Systematics, re-baselined on the structured bubble template

**Baseline changed 2026-09-07 (Kev):** Oscar Macias' structured all-sky
bubble map replaces the Ackermann 2014 flat mask as the baseline template,
on physical grounds — it is the bubbles as the data show them rather than a
flat geometric mask — and because it fits better with the same parameter count
at 4–7 and 35–60 GeV (`bubble_macias_2026-09-07.md`). The flat mask stays in
the table as the variant `bubble_flat`. **Supersedes `systematics_2026-09-04.md`
and the baseline in `halo_search_v3_2026-09-04.md`**, both kept as the record of
the flat-mask baseline. **Data:** `runs/fff2de5bc7e0`, disk-excluded, profile
estimator. **Script:** `scripts/systematics.py` (default `--bubble macias`).
**Output:** `results/fits/systematics_2026-09-07.json`.

## Baseline (structured bubble, disk-excluded)

| bin | E (GeV) | halo counts | 68% | flat-mask value | ratio |
|---|---|---|---|---|---|
| 0 | 1.51 | 106,612 | [101.1k, 112.0k] | 150,226 | 0.71 |
| 1 | 2.55 | 56,683 | [52.2k, 60.2k] | 91,952 | 0.62 |
| 2 | 4.31 | 31,305 | [29.1k, 33.0k] | 50,570 | 0.62 |
| 3 | 7.28 | 13,639 | [12.3k, 14.9k] | 26,187 | 0.52 |
| 4 | 12.29 | 8,697 | [7.8k, 9.6k] | 14,660 | 0.59 |
| **5** | **20.76** | **4,612** | **[4,026, 5,237]** | 7,705 | **0.60** |
| 6 | 35.06 | 267 | [−315, 831] | 3,843 | 0.07 |
| 7 | 59.22 | 185 | [43, 685] | 1,494 | 0.12 |
| 8 | 100.02 | 71 | [−106, 273] | 555 | 0.13 |

TS at 20.8 GeV is 24 (flat: 69). The SED proxy (counts × E) declines from
1.5 to 7 GeV, is flat at ~100 from 7 to 21 GeV, and falls to zero at 35 GeV:
**no 20 GeV peak, and no halo-shaped excess above 35 GeV** (consistent with
zero at 35 and 100 GeV; 59 GeV marginal). Injection slope 0.997–1.013.

## Halo amplitude vs background choice, 20.8 GeV

| variant | bin 0 (1.5) | bin 2 (4.3) | **bin 5 (20.8)** | 68% | bin 8 (100) |
|---|---|---|---|---|---|
| **baseline (structured)** | 106,612 | 31,305 | **4,612** | [4,026, 5,237] | 71 |
| bubble_flat (the old baseline) | 150,226 | 50,570 | 7,705 | [7,306, 8,290] | 555 |
| bubble_NS (structured, N/S split) | 99,343 | 29,499 | 4,443 | [3,844, 5,065] | 38 |
| no Loop I | 130,842 | 35,519 | 4,612 | [4,127, 5,332] | 66 |
| no bubble | 242,300 | 90,906 | **15,242** | [14,745, 15,772] | 1,609 |
| no Loop I + no bubble | 227,762 | 88,497 | 14,961 | [14,465, 15,473] | 1,604 |
| no dust | 142,957 | 43,346 | 5,962 | [5,363, 6,556] | 190 |
| PS masked (95%) | 2,941 | 16,848 | 4,326 | [3,748, 4,883] | −79 |
| PS masked (1° fixed) | 22,869 | 7,186 | **1,346** | [949, 1,742] | −63 |

**Envelope at 20.8 GeV:** 4,326–15,242 = **3.5×** over the physically
motivated variants (was 2.1× on the flat baseline, 7,424–15,242); 1,346–15,242
= **11×** including the 1° point-source mask (was 5.2×). The envelope widened
because the flat mask, now a variant, sits 1.7× above the new baseline while
the no-bubble variant did not move. The 1° mask is a *geometry* change that
removes 44% of the pixels (open item 5, `STATE_2026-09-04`) and is reported
separately for that reason; the 95% PS mask agrees with the baseline (4,326
vs 4,612).

## Quadrant consistency (sum of quadrants / baseline; 1.0 = consistent)

| variant | bin 0 | bin 2 | bin 5 |
|---|---|---|---|
| **baseline (structured)** | **−0.21 (negative)** | −0.10 (negative) | **0.57 (one negative)** |
| bubble_flat | 0.55 (negative) | 0.76 | 0.84 |
| bubble_NS | −0.26 (negative) | −0.14 (negative) | 0.59 (negative) |
| no Loop I | 0.10 (negative) | −0.10 (negative) | 0.54 (negative) |
| no bubble | 0.93 | 1.06 | 0.93 |
| no Loop I + no bubble | 1.10 | 1.10 | 0.95 |
| no dust | −0.29 (negative) | −0.43 (negative) | 0.37 (negative) |
| PS masked (95%) | 1.53 (negative) | −0.75 (negative) | 0.18 (negative) |
| PS masked (1°) | −0.25 (negative) | −0.76 (negative) | −0.60 (negative) |

🚨 **On the structured baseline the quadrant test fails at every energy, and
in every variant that keeps a bubble template.** At 20.8 GeV the four quadrant
amplitudes are [753, −313, 1,446, 707] against a joint 4,612: one is negative
and the sum is 0.57 of the joint value. At 1.5 and 4.3 GeV the quadrants sum
to the *opposite sign* of the joint amplitude. A spherical halo cannot do
this. The only configurations that pass are the ones with no bubble template
at all — which is the configuration in which the halo template is known to be
absorbing the bubbles' photons (`systematics_2026-09-04`, "correcting an
over-read"). Quadrant consistency was necessary-not-sufficient on the flat
baseline; on the structured baseline it is not even met.

`no_dust` still fails the injection gate at bin 4 (slope 0.68; open item 4).
`no_loopI` reaches slope 1.09 in one bin; every other variant is 0.997–1.04.

## Bottom line, re-baselined

- **No 20 GeV peak.** Flat 7–21 GeV in counts × E, declining below, zero above
  35 GeV.
- **The amplitude at 20.8 GeV is 4,612 ± 13% statistically and undetermined to
  a factor 3.5–11 systematically**, dominated by the bubble treatment and the
  point-source treatment.
- **The baseline halo is not quadrant-consistent.** What the joint fit calls a
  halo at 20.8 GeV is a common amplitude over four quadrants that disagree by
  more than their sum, one of them negative.
- **There is no halo-shaped excess above 35 GeV** with a structured bubble
  template.
- The categorical diffuse marginalization (Galp21 at |b| ≤ 60) remains the
  fix for the amplitude systematic; nothing here changes that argument, and
  the Galp21 band-forecast leverage is unaffected by the bubble choice
  (`band_forecast_bubble_2026-09-07.json`).
