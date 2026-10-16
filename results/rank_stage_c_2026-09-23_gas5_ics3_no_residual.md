# Stage C — 2026-09-23 (gas5_ics3_no_residual)

Consensus selection across 129 model rows (gas5_ics3) by the **fraction-of-rows rule**: a candidate is eligible if it clears its own BIC threshold (p·ln 1080 = 63 per template) in ≥ 80% of rows; the eligible candidate that is the largest improvement in the most rows is selected. Median Δ is shown for reference only. Selected: **bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic**

| step | chosen | clears BIC in | best among eligible in | median Δ (ref.) | runner-up (clears / best) |
|---|---|---|---|---|---|
| 1 | bubble_flat | 100% | 64% | 15,305 | bubble_structured (100% / 36%) |
| 2 | erosita_unif_outer_NS | 100% | 51% | 10,105 | erosita_unif_NS (100% / 34%) |
| 3 | loopI_haslam_raw | 100% | 95% | 882 | loopI_shells+nps (100% / 5%) |
| 4 | sun_moon_ic | 100% | 100% | 347 | loopIV (100% / 0%) |
| 5 | (stop) loopIV | 51% | 0% | 64 | loopII (0% / 0%) |

## Each model's own path

| model | variant | path |
|---|---|---|
| SL_Z10_R20_T100000_C2 | gas5_ics3 | bubble_flat → erosita_unif_NS → loopI_haslam_raw → loopIV → sun_moon_ic |
| SL_Z10_R20_T100000_C5 | gas5_ics3 | bubble_flat → erosita_unif_NS → loopI_haslam_raw → loopIV |
| SL_Z10_R20_T150_C2 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SL_Z10_R20_T150_C5 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SL_Z10_R30_T100000_C2 | gas5_ics3 | bubble_flat → erosita_unif_NS → loopI_haslam_raw → loopIV → sun_moon_ic |
| SL_Z10_R30_T100000_C5 | gas5_ics3 | bubble_flat → erosita_unif_NS → loopI_haslam_raw → loopIV |
| SL_Z10_R30_T150_C2 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SL_Z10_R30_T150_C5 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SL_Z4_R20_T100000_C2 | gas5_ics3 | bubble_flat → erosita_unif_NS → loopI_haslam_raw → loopIV → sun_moon_ic |
| SL_Z4_R20_T100000_C5 | gas5_ics3 | bubble_flat → erosita_unif_NS → loopI_haslam_raw → loopIV → sun_moon_ic |
| SL_Z4_R20_T150_C2 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SL_Z4_R20_T150_C5 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SL_Z4_R30_T100000_C2 | gas5_ics3 | bubble_flat → erosita_unif_NS → loopI_haslam_raw → loopIV → sun_moon_ic |
| SL_Z4_R30_T100000_C5 | gas5_ics3 | bubble_flat → erosita_unif_NS → loopI_haslam_raw → loopIV → sun_moon_ic |
| SL_Z4_R30_T150_C2 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SL_Z4_R30_T150_C5 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SL_Z6_R20_T100000_C2 | gas5_ics3 | bubble_flat → erosita_unif_NS → loopI_haslam_raw → loopIV → sun_moon_ic |
| SL_Z6_R20_T100000_C5 | gas5_ics3 | bubble_flat → erosita_outer_NS → loopI_haslam_raw → sun_moon_ic |
| SL_Z6_R20_T150_C2 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SL_Z6_R20_T150_C5 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SL_Z6_R30_T100000_C2 | gas5_ics3 | bubble_flat → erosita_unif_NS → loopI_haslam_raw → loopIV → sun_moon_ic |
| SL_Z6_R30_T100000_C5 | gas5_ics3 | bubble_flat → erosita_outer_NS → loopI_haslam_raw → sun_moon_ic |
| SL_Z6_R30_T150_C2 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SL_Z6_R30_T150_C5 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SL_Z8_R20_T100000_C2 | gas5_ics3 | bubble_flat → erosita_unif_NS → loopI_haslam_raw → loopIV → sun_moon_ic |
| SL_Z8_R20_T100000_C5 | gas5_ics3 | bubble_flat → erosita_unif_NS → loopI_haslam_raw → loopIV |
| SL_Z8_R20_T150_C2 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SL_Z8_R20_T150_C5 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SL_Z8_R30_T100000_C2 | gas5_ics3 | bubble_flat → erosita_unif_NS → loopI_haslam_raw → loopIV → sun_moon_ic |
| SL_Z8_R30_T100000_C5 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic |
| SL_Z8_R30_T150_C2 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SL_Z8_R30_T150_C5 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SO_Z10_R20_T100000_C2 | gas5_ics3 | bubble_flat → erosita_unif_NS → loopI_haslam_raw → loopIV → sun_moon_ic |
| SO_Z10_R20_T100000_C5 | gas5_ics3 | bubble_flat → erosita_unif_NS → loopI_haslam_raw → loopIV |
| SO_Z10_R20_T150_C2 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SO_Z10_R20_T150_C5 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SO_Z10_R30_T100000_C2 | gas5_ics3 | bubble_flat → erosita_unif_NS → loopI_haslam_raw → loopIV → sun_moon_ic |
| SO_Z10_R30_T100000_C5 | gas5_ics3 | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic |
| SO_Z10_R30_T150_C2 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SO_Z10_R30_T150_C5 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SO_Z4_R20_T100000_C2 | gas5_ics3 | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic |
| SO_Z4_R20_T100000_C5 | gas5_ics3 | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic |
| SO_Z4_R20_T150_C2 | gas5_ics3 | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SO_Z4_R20_T150_C5 | gas5_ics3 | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SO_Z4_R30_T100000_C2 | gas5_ics3 | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic |
| SO_Z4_R30_T100000_C5 | gas5_ics3 | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic |
| SO_Z4_R30_T150_C2 | gas5_ics3 | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SO_Z4_R30_T150_C5 | gas5_ics3 | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SO_Z6_R20_T100000_C2 | gas5_ics3 | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic |
| SO_Z6_R20_T100000_C5 | gas5_ics3 | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic |
| SO_Z6_R20_T150_C2 | gas5_ics3 | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SO_Z6_R20_T150_C5 | gas5_ics3 | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SO_Z6_R30_T100000_C2 | gas5_ics3 | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic |
| SO_Z6_R30_T100000_C5 | gas5_ics3 | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic |
| SO_Z6_R30_T150_C2 | gas5_ics3 | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SO_Z6_R30_T150_C5 | gas5_ics3 | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SO_Z8_R20_T100000_C2 | gas5_ics3 | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic |
| SO_Z8_R20_T100000_C5 | gas5_ics3 | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic |
| SO_Z8_R20_T150_C2 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SO_Z8_R20_T150_C5 | gas5_ics3 | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SO_Z8_R30_T100000_C2 | gas5_ics3 | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic |
| SO_Z8_R30_T100000_C5 | gas5_ics3 | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic |
| SO_Z8_R30_T150_C2 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SO_Z8_R30_T150_C5 | gas5_ics3 | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SS_Z10_R20_T100000_C2 | gas5_ics3 | bubble_flat → erosita_unif_NS → loopI_haslam_raw → loopIV → sun_moon_ic |
| SS_Z10_R20_T100000_C5 | gas5_ics3 | bubble_flat → erosita_unif_NS → loopI_haslam_raw → loopIV |
| SS_Z10_R20_T150_C2 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SS_Z10_R20_T150_C5 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SS_Z10_R30_T100000_C2 | gas5_ics3 | bubble_flat → erosita_unif_NS → loopI_haslam_raw → loopIV → sun_moon_ic |
| SS_Z10_R30_T100000_C5 | gas5_ics3 | bubble_flat → erosita_unif_NS → loopI_haslam_raw → loopIV |
| SS_Z10_R30_T150_C2 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SS_Z10_R30_T150_C5 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SS_Z4_R20_T100000_C2 | gas5_ics3 | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic |
| SS_Z4_R20_T100000_C5 | gas5_ics3 | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic |
| SS_Z4_R20_T150_C2 | gas5_ics3 | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SS_Z4_R20_T150_C5 | gas5_ics3 | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SS_Z4_R30_T100000_C2 | gas5_ics3 | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic |
| SS_Z4_R30_T100000_C5 | gas5_ics3 | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic |
| SS_Z4_R30_T150_C2 | gas5_ics3 | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SS_Z4_R30_T150_C5 | gas5_ics3 | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SS_Z6_R20_T100000_C2 | gas5_ics3 | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic |
| SS_Z6_R20_T100000_C5 | gas5_ics3 | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic |
| SS_Z6_R20_T150_C2 | gas5_ics3 | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SS_Z6_R20_T150_C5 | gas5_ics3 | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SS_Z6_R30_T100000_C2 | gas5_ics3 | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic |
| SS_Z6_R30_T100000_C5 | gas5_ics3 | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic |
| SS_Z6_R30_T150_C2 | gas5_ics3 | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SS_Z6_R30_T150_C5 | gas5_ics3 | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SS_Z8_R20_T100000_C2 | gas5_ics3 | bubble_flat → erosita_unif_NS → loopI_haslam_raw → loopIV → sun_moon_ic |
| SS_Z8_R20_T100000_C5 | gas5_ics3 | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic |
| SS_Z8_R20_T150_C2 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SS_Z8_R20_T150_C5 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SS_Z8_R30_T100000_C2 | gas5_ics3 | bubble_flat → erosita_unif_NS → loopI_haslam_raw → loopIV → sun_moon_ic |
| SS_Z8_R30_T100000_C5 | gas5_ics3 | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic |
| SS_Z8_R30_T150_C2 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SS_Z8_R30_T150_C5 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SY_Z10_R20_T100000_C2 | gas5_ics3 | bubble_flat → erosita_unif_NS → loopI_haslam_raw → loopIV → sun_moon_ic |
| SY_Z10_R20_T100000_C5 | gas5_ics3 | bubble_flat → erosita_unif_NS → loopI_haslam_raw → loopIV |
| SY_Z10_R20_T150_C2 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SY_Z10_R20_T150_C5 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SY_Z10_R30_T100000_C2 | gas5_ics3 | bubble_flat → erosita_unif_NS → loopI_haslam_raw → loopIV → sun_moon_ic |
| SY_Z10_R30_T100000_C5 | gas5_ics3 | bubble_flat → erosita_unif_NS → loopI_haslam_raw → loopIV |
| SY_Z10_R30_T150_C2 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SY_Z10_R30_T150_C5 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SY_Z4_R20_T100000_C2 | gas5_ics3 | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic |
| SY_Z4_R20_T100000_C5 | gas5_ics3 | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic |
| SY_Z4_R20_T150_C2 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SY_Z4_R20_T150_C5 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SY_Z4_R30_T100000_C2 | gas5_ics3 | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic |
| SY_Z4_R30_T100000_C5 | gas5_ics3 | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic |
| SY_Z4_R30_T150_C2 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SY_Z4_R30_T150_C5 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SY_Z6_R20_T100000_C2 | gas5_ics3 | bubble_flat → erosita_unif_NS → loopI_haslam_raw → loopIV → sun_moon_ic |
| SY_Z6_R20_T100000_C5 | gas5_ics3 | bubble_flat → erosita_outer_NS → loopI_haslam_raw → sun_moon_ic |
| SY_Z6_R20_T150_C2 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SY_Z6_R20_T150_C5 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SY_Z6_R30_T100000_C2 | gas5_ics3 | bubble_flat → erosita_unif_NS → loopI_haslam_raw → loopIV → sun_moon_ic |
| SY_Z6_R30_T100000_C5 | gas5_ics3 | bubble_flat → erosita_outer_NS → loopI_haslam_raw → sun_moon_ic |
| SY_Z6_R30_T150_C2 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SY_Z6_R30_T150_C5 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SY_Z8_R20_T100000_C2 | gas5_ics3 | bubble_flat → erosita_unif_NS → loopI_haslam_raw → loopIV → sun_moon_ic |
| SY_Z8_R20_T100000_C5 | gas5_ics3 | bubble_flat → erosita_unif_NS → loopI_haslam_raw → loopIV |
| SY_Z8_R20_T150_C2 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SY_Z8_R20_T150_C5 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SY_Z8_R30_T100000_C2 | gas5_ics3 | bubble_flat → erosita_unif_NS → loopI_haslam_raw → loopIV → sun_moon_ic |
| SY_Z8_R30_T100000_C5 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic |
| SY_Z8_R30_T150_C2 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| SY_Z8_R30_T150_C5 | gas5_ics3 | bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV |
| a2020 | - | bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic → loopIV |

## Ranking with the consensus set (bubble_flat, erosita_unif_outer_NS, loopI_haslam_raw, sun_moon_ic)

| rank | model | variant | n | Δ(−2lnL) vs Totani baseline | halo @ 20.8 | TS @ 20.8 | quad Σ/joint | any neg | inj min | N / S |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | gll_iem_v07 | - | 8 | - | 2,025 [1,442, 2,587] | 12.5 | 0.83 | yes | 0.99 | 652 / 1,032 |
| 2 | a2020 | - | 21 | - | 7,424 [6,877, 7,967] | 49.8 | 0.93 | no | 1.00 | 3,147 / 3,740 |
| 3 | SS_Z10_R30_T150_C5 | gas5_ics3 | 15 | - | 5,962 [4,423, 6,726] | 11.3 | 0.61 | no | 1.00 | 2,105 / 1,542 |
| 4 | SY_Z10_R30_T150_C5 | gas5_ics3 | 15 | - | 5,793 [4,119, 6,679] | 8.8 | 0.59 | no | 0.99 | 1,972 / 1,450 |
| 5 | SL_Z10_R30_T150_C5 | gas5_ics3 | 15 | - | 5,906 [4,171, 6,783] | 8.2 | 0.54 | no | 0.99 | 1,853 / 1,331 |
| 6 | SY_Z10_R20_T150_C5 | gas5_ics3 | 15 | - | 5,793 [4,183, 6,677] | 9.1 | 0.60 | no | 1.00 | 1,972 / 1,478 |
| 7 | SS_Z10_R20_T150_C5 | gas5_ics3 | 15 | - | 5,962 [4,465, 6,717] | 11.5 | 0.62 | no | 1.00 | 2,105 / 1,583 |
| 8 | SL_Z10_R20_T150_C5 | gas5_ics3 | 15 | - | 5,906 [4,230, 6,776] | 8.5 | 0.55 | no | 0.99 | 1,853 / 1,372 |
| 9 | SO_Z10_R30_T150_C5 | gas5_ics3 | 15 | - | 6,356 [4,914, 6,997] | 11.4 | 0.56 | no | 0.98 | 2,090 / 1,488 |
| 10 | SS_Z8_R30_T150_C5 | gas5_ics3 | 15 | - | 6,018 [4,797, 6,744] | 14.3 | 0.68 | no | 0.99 | 2,256 / 1,863 |
| 11 | SY_Z8_R30_T150_C5 | gas5_ics3 | 15 | - | 5,849 [4,514, 6,703] | 11.5 | 0.66 | no | 1.00 | 2,123 / 1,759 |
| 12 | SS_Z6_R30_T150_C5 | gas5_ics3 | 15 | - | 6,131 [5,187, 6,810] | 20.4 | 0.80 | no | 0.99 | 2,618 / 2,291 |
| 13 | SY_Z6_R30_T150_C5 | gas5_ics3 | 15 | - | 6,018 [4,977, 6,792] | 17.3 | 0.76 | no | 0.98 | 2,379 / 2,213 |
| 14 | SY_Z6_R20_T150_C5 | gas5_ics3 | 15 | - | 6,018 [4,987, 6,792] | 17.4 | 0.77 | no | 0.98 | 2,379 / 2,227 |
| 15 | SO_Z10_R20_T150_C5 | gas5_ics3 | 15 | - | 6,356 [4,951, 6,989] | 11.6 | 0.57 | no | 1.00 | 2,090 / 1,516 |
| 16 | SS_Z6_R20_T150_C5 | gas5_ics3 | 15 | - | 6,131 [5,197, 6,809] | 20.5 | 0.80 | no | 1.00 | 2,618 / 2,291 |
| 17 | SL_Z8_R30_T150_C5 | gas5_ics3 | 15 | - | 5,906 [4,529, 6,774] | 10.5 | 0.62 | no | 0.98 | 1,989 / 1,653 |
| 18 | SS_Z8_R20_T150_C5 | gas5_ics3 | 15 | - | 6,018 [4,824, 6,742] | 14.4 | 0.68 | no | 0.99 | 2,256 / 1,864 |
| 19 | SL_Z6_R30_T150_C5 | gas5_ics3 | 15 | - | 6,018 [5,000, 6,797] | 16.0 | 0.72 | no | 1.00 | 2,200 / 2,148 |
| 20 | SY_Z8_R20_T150_C5 | gas5_ics3 | 15 | - | 5,849 [4,540, 6,702] | 11.7 | 0.67 | no | 1.00 | 2,138 / 1,772 |
| 21 | SL_Z6_R20_T150_C5 | gas5_ics3 | 15 | - | 6,018 [5,007, 6,793] | 16.1 | 0.72 | no | 0.99 | 2,200 / 2,161 |
| 22 | SL_Z8_R20_T150_C5 | gas5_ics3 | 15 | - | 5,906 [4,550, 6,771] | 10.6 | 0.62 | no | 0.99 | 1,989 / 1,654 |
| 23 | SO_Z8_R30_T150_C5 | gas5_ics3 | 15 | - | 6,412 [5,244, 6,977] | 14.0 | 0.60 | no | 1.00 | 2,075 / 1,769 |
| 24 | SO_Z8_R20_T150_C5 | gas5_ics3 | 15 | - | 6,412 [5,261, 6,972] | 14.1 | 0.60 | no | 0.99 | 2,075 / 1,796 |
| 25 | SO_Z6_R30_T150_C5 | gas5_ics3 | 15 | - | 6,412 [5,627, 7,000] | 20.0 | 0.71 | no | 0.98 | 2,332 / 2,237 |
| 26 | SO_Z6_R20_T150_C5 | gas5_ics3 | 15 | - | 6,412 [5,628, 6,994] | 20.1 | 0.71 | no | 0.98 | 2,317 / 2,250 |
| 27 | SS_Z10_R30_T150_C2 | gas5_ics3 | 15 | - | 5,849 [4,314, 6,889] | 10.8 | 0.64 | no | 0.99 | 2,165 / 1,583 |
| 28 | SY_Z10_R30_T150_C2 | gas5_ics3 | 15 | - | 5,793 [4,012, 6,881] | 8.4 | 0.61 | no | 0.98 | 2,032 / 1,491 |
| 29 | SS_Z8_R30_T150_C2 | gas5_ics3 | 15 | - | 6,131 [4,665, 6,930] | 13.4 | 0.67 | no | 0.98 | 2,301 / 1,837 |
| 30 | SL_Z10_R30_T150_C2 | gas5_ics3 | 15 | - | 5,962 [4,066, 6,956] | 7.8 | 0.55 | no | 0.99 | 1,913 / 1,358 |
| 31 | SY_Z8_R30_T150_C2 | gas5_ics3 | 15 | - | 6,074 [4,411, 6,918] | 10.9 | 0.65 | no | 0.99 | 2,183 / 1,759 |
| 32 | SS_Z10_R20_T150_C2 | gas5_ics3 | 15 | - | 5,906 [4,352, 6,880] | 10.9 | 0.64 | no | 1.00 | 2,150 / 1,624 |
| 33 | SY_Z10_R20_T150_C2 | gas5_ics3 | 15 | - | 5,793 [4,064, 6,877] | 8.6 | 0.62 | no | 0.99 | 2,047 / 1,519 |
| 34 | SO_Z10_R30_T150_C2 | gas5_ics3 | 15 | - | 6,356 [4,819, 7,145] | 10.9 | 0.57 | no | 0.99 | 2,134 / 1,516 |
| 35 | SL_Z10_R20_T150_C2 | gas5_ics3 | 15 | - | 6,018 [4,118, 6,947] | 8.0 | 0.55 | no | 0.99 | 1,928 / 1,399 |
| 36 | SS_Z8_R20_T150_C2 | gas5_ics3 | 15 | - | 6,187 [4,712, 6,917] | 13.7 | 0.68 | no | 1.00 | 2,301 / 1,891 |
| 37 | SL_Z8_R30_T150_C2 | gas5_ics3 | 15 | - | 6,131 [4,477, 6,945] | 10.0 | 0.61 | no | 0.99 | 2,064 / 1,654 |
| 38 | SO_Z10_R20_T150_C2 | gas5_ics3 | 15 | - | 6,412 [4,854, 7,136] | 11.1 | 0.57 | no | 0.98 | 2,119 / 1,557 |
| 39 | SS_Z4_R30_T150_C5 | gas5_ics3 | 15 | - | 6,131 [5,484, 6,796] | 37.5 | 0.84 | no | 1.00 | 2,365 / 2,794 |
| 40 | SS_Z6_R30_T150_C2 | gas5_ics3 | 15 | - | 6,243 [5,157, 6,968] | 19.6 | 0.79 | no | 1.00 | 2,618 / 2,332 |
| 41 | SS_Z4_R20_T150_C5 | gas5_ics3 | 15 | - | 6,131 [5,485, 6,796] | 37.6 | 0.84 | no | 1.00 | 2,350 / 2,794 |
| 42 | SY_Z8_R20_T150_C2 | gas5_ics3 | 15 | - | 6,074 [4,468, 6,902] | 11.1 | 0.66 | no | 0.99 | 2,198 / 1,813 |
| 43 | SY_Z6_R30_T150_C2 | gas5_ics3 | 15 | - | 6,187 [4,985, 6,979] | 16.6 | 0.75 | no | 0.96 | 2,379 / 2,267 |
| 44 | SY_Z4_R30_T150_C5 | gas5_ics3 | 15 | - | 6,074 [5,348, 6,826] | 33.4 | 0.82 | no | 1.00 | 2,260 / 2,703 |
| 45 | SS_Z6_R20_T150_C2 | gas5_ics3 | 15 | - | 6,243 [5,169, 6,967] | 19.7 | 0.79 | no | 1.00 | 2,618 / 2,332 |
| 46 | SY_Z6_R20_T150_C2 | gas5_ics3 | 15 | - | 6,187 [4,996, 6,979] | 16.7 | 0.75 | no | 0.97 | 2,379 / 2,267 |
| 47 | SY_Z4_R20_T150_C5 | gas5_ics3 | 15 | - | 6,074 [5,348, 6,825] | 33.5 | 0.82 | no | 1.00 | 2,260 / 2,703 |
| 48 | SL_Z8_R20_T150_C2 | gas5_ics3 | 15 | - | 6,131 [4,514, 6,932] | 10.2 | 0.61 | no | 0.99 | 2,064 / 1,695 |
| 49 | SL_Z6_R30_T150_C2 | gas5_ics3 | 15 | - | 6,243 [5,037, 6,937] | 15.4 | 0.70 | no | 0.99 | 2,200 / 2,188 |
| 50 | SO_Z8_R30_T150_C2 | gas5_ics3 | 15 | - | 6,524 [5,191, 7,120] | 13.5 | 0.61 | no | 1.00 | 2,150 / 1,823 |
| 51 | SL_Z6_R20_T150_C2 | gas5_ics3 | 15 | - | 6,243 [5,044, 6,934] | 15.5 | 0.71 | no | 0.98 | 2,200 / 2,201 |
| 52 | SO_Z8_R20_T150_C2 | gas5_ics3 | 15 | - | 6,524 [5,214, 7,115] | 13.6 | 0.61 | no | 1.00 | 2,150 / 1,837 |
| 53 | SL_Z4_R30_T150_C5 | gas5_ics3 | 15 | - | 6,131 [5,369, 6,742] | 31.5 | 0.78 | no | 1.00 | 2,111 / 2,677 |
| 54 | SL_Z4_R20_T150_C5 | gas5_ics3 | 15 | - | 6,131 [5,369, 6,741] | 31.6 | 0.78 | no | 1.00 | 2,111 / 2,677 |
| 55 | SO_Z6_R30_T150_C2 | gas5_ics3 | 15 | - | 6,580 [5,601, 7,126] | 19.4 | 0.71 | no | 1.00 | 2,392 / 2,277 |
| 56 | SO_Z6_R20_T150_C2 | gas5_ics3 | 15 | - | 6,580 [5,605, 7,124] | 19.5 | 0.71 | no | 1.00 | 2,392 / 2,277 |
| 57 | SO_Z4_R30_T150_C5 | gas5_ics3 | 15 | - | 6,524 [5,926, 7,067] | 38.4 | 0.80 | no | 1.00 | 2,364 / 2,861 |
| 58 | SO_Z4_R20_T150_C5 | gas5_ics3 | 15 | - | 6,524 [5,925, 7,066] | 38.4 | 0.80 | no | 1.00 | 2,364 / 2,861 |
| 59 | SS_Z4_R30_T150_C2 | gas5_ics3 | 15 | - | 6,299 [5,573, 6,940] | 36.8 | 0.82 | no | 1.00 | 2,364 / 2,793 |
| 60 | SS_Z4_R20_T150_C2 | gas5_ics3 | 15 | - | 6,299 [5,574, 6,940] | 36.8 | 0.82 | no | 1.00 | 2,364 / 2,780 |
| 61 | SY_Z4_R30_T150_C2 | gas5_ics3 | 15 | - | 6,243 [5,479, 6,976] | 32.8 | 0.79 | no | 0.99 | 2,245 / 2,702 |
| 62 | SY_Z4_R20_T150_C2 | gas5_ics3 | 15 | - | 6,243 [5,479, 6,975] | 32.8 | 0.79 | no | 0.99 | 2,245 / 2,702 |
| 63 | SL_Z4_R30_T150_C2 | gas5_ics3 | 15 | - | 6,299 [5,522, 6,865] | 31.1 | 0.76 | no | 1.00 | 2,096 / 2,664 |
| 64 | SL_Z4_R20_T150_C2 | gas5_ics3 | 15 | - | 6,299 [5,522, 6,864] | 31.1 | 0.76 | no | 1.00 | 2,096 / 2,664 |
| 65 | SY_Z8_R30_T100000_C5 | gas5_ics3 | 15 | - | 5,568 [4,274, 6,431] | 12.9 | 0.68 | no | 0.97 | 1,987 / 1,804 |
| 66 | SL_Z8_R30_T100000_C5 | gas5_ics3 | 15 | - | 5,624 [4,135, 6,501] | 10.8 | 0.56 | no | 0.99 | 1,420 / 1,712 |
| 67 | SO_Z4_R30_T150_C2 | gas5_ics3 | 15 | - | 6,637 [5,988, 7,176] | 37.8 | 0.78 | no | 0.99 | 2,349 / 2,848 |
| 68 | SO_Z4_R20_T150_C2 | gas5_ics3 | 15 | - | 6,637 [5,987, 7,176] | 37.8 | 0.78 | no | 1.00 | 2,349 / 2,848 |
| 69 | SS_Z8_R30_T100000_C5 | gas5_ics3 | 15 | - | 5,849 [5,013, 6,534] | 16.7 | 0.74 | no | 1.00 | 2,405 / 1,896 |
| 70 | SY_Z8_R30_T100000_C2 | gas5_ics3 | 15 | - | 5,737 [4,161, 6,628] | 12.1 | 0.68 | no | 0.99 | 2,016 / 1,873 |
| 71 | SS_Z8_R30_T100000_C2 | gas5_ics3 | 15 | - | 5,962 [4,860, 6,690] | 15.7 | 0.73 | no | 0.99 | 2,360 / 1,978 |
| 72 | SL_Z8_R30_T100000_C2 | gas5_ics3 | 15 | - | 5,849 [4,111, 6,722] | 10.2 | 0.57 | no | 0.99 | 1,555 / 1,781 |
| 73 | SY_Z8_R20_T100000_C2 | gas5_ics3 | 15 | - | 5,793 [4,205, 6,632] | 12.3 | 0.68 | no | 0.99 | 2,016 / 1,900 |
| 74 | SY_Z10_R30_T100000_C5 | gas5_ics3 | 15 | - | 5,456 [3,955, 6,359] | 10.4 | 0.61 | no | 0.99 | 1,657 / 1,657 |
| 75 | SL_Z8_R20_T100000_C2 | gas5_ics3 | 15 | - | 5,849 [4,160, 6,724] | 10.4 | 0.57 | no | 0.99 | 1,555 / 1,808 |
| 76 | SS_Z8_R20_T100000_C2 | gas5_ics3 | 15 | - | 5,962 [4,883, 6,691] | 15.9 | 0.73 | no | 0.99 | 2,360 / 1,993 |
| 77 | SL_Z10_R30_T100000_C5 | gas5_ics3 | 15 | - | 5,512 [3,734, 6,442] | 8.5 | 0.47 | no | 0.99 | 1,105 / 1,462 |
| 78 | SS_Z10_R30_T100000_C5 | gas5_ics3 | 15 | - | 5,793 [4,760, 6,474] | 13.8 | 0.69 | no | 0.99 | 2,224 / 1,762 |
| 79 | SS_Z8_R20_T100000_C5 | gas5_ics3 | 15 | - | 5,849 [5,031, 6,529] | 16.9 | 0.75 | no | 0.99 | 2,405 / 1,964 |
| 80 | SY_Z8_R20_T100000_C5 | gas5_ics3 | 15 | - | 5,568 [4,320, 6,422] | 13.2 | 0.68 | no | 0.99 | 1,912 / 1,886 |
| 81 | SO_Z8_R30_T100000_C2 | gas5_ics3 | 15 | - | 6,356 [5,240, 6,966] | 15.6 | 0.65 | no | 0.99 | 2,179 / 1,925 |
| 82 | SS_Z10_R30_T100000_C2 | gas5_ics3 | 15 | - | 5,681 [4,561, 6,628] | 12.9 | 0.70 | no | 1.00 | 2,179 / 1,805 |
| 83 | SY_Z10_R30_T100000_C2 | gas5_ics3 | 15 | - | 5,512 [3,791, 6,559] | 9.6 | 0.62 | no | 0.98 | 1,731 / 1,713 |
| 84 | SL_Z10_R30_T100000_C2 | gas5_ics3 | 15 | - | 5,624 [3,662, 6,669] | 7.9 | 0.50 | no | 0.99 | 1,284 / 1,530 |
| 85 | SY_Z6_R30_T100000_C5 | gas5_ics3 | 15 | - | 5,737 [4,750, 6,535] | 19.3 | 0.84 | no | 0.98 | 2,543 / 2,300 |
| 86 | SL_Z8_R20_T100000_C5 | gas5_ics3 | 15 | - | 5,624 [4,199, 6,485] | 11.1 | 0.56 | no | 0.99 | 1,361 / 1,807 |
| 87 | SY_Z10_R20_T100000_C5 | gas5_ics3 | 15 | - | 5,512 [4,013, 6,364] | 10.7 | 0.61 | no | 0.99 | 1,627 / 1,712 |
| 88 | SS_Z6_R30_T100000_C5 | gas5_ics3 | 15 | - | 5,962 [5,320, 6,609] | 23.5 | 0.87 | no | 1.00 | 2,827 / 2,364 |
| 89 | SO_Z10_R30_T100000_C5 | gas5_ics3 | 15 | - | 6,243 [5,063, 6,817] | 13.8 | 0.60 | no | 1.00 | 2,029 / 1,697 |
| 90 | SL_Z10_R20_T100000_C5 | gas5_ics3 | 15 | - | 5,512 [3,817, 6,445] | 8.8 | 0.48 | no | 0.99 | 1,105 / 1,542 |
| 91 | SY_Z6_R20_T100000_C5 | gas5_ics3 | 15 | - | 5,737 [4,761, 6,537] | 19.4 | 0.84 | no | 0.98 | 2,543 / 2,287 |
| 92 | SL_Z6_R30_T100000_C5 | gas5_ics3 | 15 | - | 5,737 [4,735, 6,589] | 16.9 | 0.74 | no | 0.99 | 2,021 / 2,222 |
| 93 | SO_Z10_R30_T100000_C2 | gas5_ics3 | 15 | - | 6,243 [4,910, 6,944] | 13.0 | 0.60 | no | 1.00 | 1,953 / 1,792 |
| 94 | SS_Z10_R20_T100000_C5 | gas5_ics3 | 15 | - | 5,793 [4,788, 6,474] | 14.1 | 0.70 | no | 0.99 | 2,224 / 1,804 |
| 95 | SS_Z6_R20_T100000_C5 | gas5_ics3 | 15 | - | 5,962 [5,324, 6,610] | 23.5 | 0.87 | no | 0.99 | 2,827 / 2,364 |
| 96 | SL_Z6_R20_T100000_C5 | gas5_ics3 | 15 | - | 5,737 [4,743, 6,590] | 16.9 | 0.74 | no | 0.99 | 2,021 / 2,222 |
| 97 | SS_Z6_R30_T100000_C2 | gas5_ics3 | 15 | - | 6,018 [5,257, 6,753] | 22.5 | 0.86 | no | 0.99 | 2,782 / 2,420 |
| 98 | SO_Z8_R30_T100000_C5 | gas5_ics3 | 15 | - | 6,243 [5,342, 6,833] | 16.5 | 0.65 | no | 1.00 | 2,165 / 1,910 |
| 99 | SO_Z8_R20_T100000_C2 | gas5_ics3 | 15 | - | 6,356 [5,265, 6,956] | 15.8 | 0.65 | no | 0.99 | 2,165 / 1,993 |
| 100 | SY_Z10_R20_T100000_C2 | gas5_ics3 | 15 | - | 5,568 [3,862, 6,560] | 9.9 | 0.63 | no | 0.98 | 1,716 / 1,781 |
| 101 | SY_Z6_R30_T100000_C2 | gas5_ics3 | 15 | - | 5,906 [4,806, 6,718] | 18.5 | 0.83 | no | 0.99 | 2,528 / 2,355 |
| 102 | SS_Z10_R20_T100000_C2 | gas5_ics3 | 15 | - | 5,737 [4,600, 6,626] | 13.2 | 0.71 | no | 0.99 | 2,179 / 1,873 |
| 103 | SS_Z6_R20_T100000_C2 | gas5_ics3 | 15 | - | 6,018 [5,263, 6,754] | 22.6 | 0.87 | no | 0.99 | 2,782 / 2,433 |
| 104 | SL_Z10_R20_T100000_C2 | gas5_ics3 | 15 | - | 5,681 [3,754, 6,668] | 8.2 | 0.50 | no | 0.99 | 1,269 / 1,598 |
| 105 | SY_Z6_R20_T100000_C2 | gas5_ics3 | 15 | - | 5,906 [4,819, 6,720] | 18.6 | 0.83 | no | 0.99 | 2,528 / 2,369 |
| 106 | SO_Z10_R20_T100000_C5 | gas5_ics3 | 15 | - | 6,243 [5,088, 6,813] | 14.0 | 0.60 | no | 1.00 | 2,014 / 1,751 |
| 107 | SO_Z10_R20_T100000_C2 | gas5_ics3 | 15 | - | 6,243 [4,944, 6,939] | 13.2 | 0.61 | no | 1.00 | 1,953 / 1,833 |
| 108 | SL_Z6_R30_T100000_C2 | gas5_ics3 | 15 | - | 5,962 [4,822, 6,796] | 16.2 | 0.72 | no | 0.99 | 2,037 / 2,277 |
| 109 | SL_Z6_R20_T100000_C2 | gas5_ics3 | 15 | - | 5,962 [4,830, 6,796] | 16.3 | 0.73 | no | 0.99 | 2,037 / 2,291 |
| 110 | SO_Z6_R30_T100000_C5 | gas5_ics3 | 15 | - | 6,356 [5,688, 6,891] | 23.0 | 0.77 | no | 1.00 | 2,542 / 2,324 |
| 111 | SO_Z6_R20_T100000_C5 | gas5_ics3 | 15 | - | 6,356 [5,690, 6,889] | 23.1 | 0.77 | no | 1.00 | 2,542 / 2,337 |
| 112 | SO_Z6_R30_T100000_C2 | gas5_ics3 | 15 | - | 6,412 [5,672, 6,981] | 22.3 | 0.77 | no | 1.00 | 2,527 / 2,406 |
| 113 | SO_Z6_R20_T100000_C2 | gas5_ics3 | 15 | - | 6,412 [5,676, 6,979] | 22.3 | 0.77 | no | 1.00 | 2,527 / 2,406 |
| 114 | SO_Z8_R20_T100000_C5 | gas5_ics3 | 15 | - | 6,243 [5,358, 6,827] | 16.7 | 0.66 | no | 1.00 | 2,150 / 1,965 |
| 115 | SS_Z4_R30_T100000_C5 | gas5_ics3 | 15 | - | 6,074 [5,536, 6,676] | 41.9 | 0.92 | no | 0.99 | 2,768 / 2,828 |
| 116 | SS_Z4_R30_T100000_C2 | gas5_ics3 | 15 | - | 6,131 [5,568, 6,810] | 40.9 | 0.90 | no | 1.00 | 2,738 / 2,801 |
| 117 | SY_Z4_R30_T100000_C2 | gas5_ics3 | 15 | - | 6,074 [5,320, 6,842] | 35.9 | 0.87 | no | 0.99 | 2,604 / 2,710 |
| 118 | SS_Z4_R20_T100000_C2 | gas5_ics3 | 15 | - | 6,131 [5,568, 6,809] | 40.9 | 0.90 | no | 1.00 | 2,738 / 2,801 |
| 119 | SY_Z4_R20_T100000_C2 | gas5_ics3 | 15 | - | 6,074 [5,321, 6,842] | 36.0 | 0.87 | no | 0.99 | 2,604 / 2,710 |
| 120 | SL_Z4_R30_T100000_C2 | gas5_ics3 | 15 | - | 6,131 [5,350, 6,833] | 33.2 | 0.83 | no | 0.99 | 2,410 / 2,672 |
| 121 | SL_Z4_R20_T100000_C2 | gas5_ics3 | 15 | - | 6,131 [5,350, 6,832] | 33.2 | 0.83 | no | 0.99 | 2,410 / 2,672 |
| 122 | SS_Z4_R20_T100000_C5 | gas5_ics3 | 15 | - | 6,018 [5,523, 6,665] | 42.2 | 0.93 | no | 1.00 | 2,798 / 2,814 |
| 123 | SO_Z4_R30_T100000_C2 | gas5_ics3 | 15 | - | 6,524 [5,955, 7,073] | 41.7 | 0.85 | no | 1.00 | 2,708 / 2,868 |
| 124 | SO_Z4_R20_T100000_C2 | gas5_ics3 | 15 | - | 6,524 [5,955, 7,072] | 41.7 | 0.85 | no | 1.00 | 2,708 / 2,868 |
| 125 | SY_Z4_R30_T100000_C5 | gas5_ics3 | 15 | - | 5,906 [5,156, 6,662] | 37.1 | 0.91 | no | 0.99 | 2,679 / 2,722 |
| 126 | SY_Z4_R20_T100000_C5 | gas5_ics3 | 15 | - | 5,906 [5,157, 6,662] | 37.1 | 0.91 | no | 0.99 | 2,679 / 2,722 |
| 127 | SL_Z4_R30_T100000_C5 | gas5_ics3 | 15 | - | 5,962 [5,165, 6,687] | 34.1 | 0.87 | no | 0.99 | 2,500 / 2,684 |
| 128 | SL_Z4_R20_T100000_C5 | gas5_ics3 | 15 | - | 5,962 [5,166, 6,687] | 34.2 | 0.87 | no | 0.99 | 2,500 / 2,684 |
| 129 | SO_Z4_R30_T100000_C5 | gas5_ics3 | 15 | - | 6,412 [5,891, 6,955] | 42.8 | 0.88 | no | 0.98 | 2,768 / 2,881 |
| 130 | SO_Z4_R20_T100000_C5 | gas5_ics3 | 15 | - | 6,412 [5,890, 6,955] | 42.9 | 0.88 | no | 0.98 | 2,768 / 2,881 |
