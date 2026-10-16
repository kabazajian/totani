# Stage C — 2026-09-21 (gas1_ics1_no_residual)

Consensus selection across 129 model rows (gas1_ics1) by the **fraction-of-rows rule**: a candidate is eligible if it clears its own BIC threshold (p·ln 1080 = 63 per template) in ≥ 80% of rows; the eligible candidate that is the largest improvement in the most rows is selected. Median Δ is shown for reference only. Selected: **bubble_structured → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic**

| step | chosen | clears BIC in | best among eligible in | median Δ (ref.) | runner-up (clears / best) |
|---|---|---|---|---|---|
| 1 | bubble_structured | 100% | 57% | 23,240 | bubble_flat (100% / 43%) |
| 2 | erosita_unif_NS | 100% | 100% | 9,892 | loopI_shells (100% / 0%) |
| 3 | loopI_haslam_raw | 100% | 100% | 1,205 | loopI_shells (88% / 0%) |
| 4 | sun_moon_ic | 100% | 72% | 208 | loopIV (100% / 28%) |
| 5 | (stop) loopIV | 60% | 0% | 110 | loopII (0% / 0%) |

## Each model's own path

| model | variant | path |
|---|---|---|
| SL_Z10_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z10_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z10_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SL_Z10_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z10_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SL_Z10_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z10_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z10_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z4_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SL_Z4_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SL_Z4_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SL_Z4_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SL_Z4_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SL_Z4_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SL_Z4_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SL_Z4_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SL_Z6_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z6_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z6_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SL_Z6_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z6_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z6_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z6_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SL_Z6_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z8_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SL_Z8_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z8_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z8_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SL_Z8_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SL_Z8_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SL_Z8_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SL_Z8_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SO_Z10_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SO_Z10_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SO_Z10_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z10_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SO_Z10_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SO_Z10_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SO_Z10_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z10_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SO_Z4_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → sun_moon_ic |
| SO_Z4_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z4_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z4_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z4_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → sun_moon_ic |
| SO_Z4_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z4_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z4_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z6_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z6_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SO_Z6_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → sun_moon_ic |
| SO_Z6_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic → loopII |
| SO_Z6_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SO_Z6_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SO_Z6_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → sun_moon_ic → loopII |
| SO_Z6_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → sun_moon_ic |
| SO_Z8_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SO_Z8_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SO_Z8_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → sun_moon_ic |
| SO_Z8_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → sun_moon_ic → loopII |
| SO_Z8_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SO_Z8_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SO_Z8_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → sun_moon_ic |
| SO_Z8_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → sun_moon_ic → loopII |
| SS_Z10_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z10_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SS_Z10_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SS_Z10_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SS_Z10_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z10_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z10_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SS_Z10_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SS_Z4_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SS_Z4_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SS_Z4_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SS_Z4_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SS_Z4_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SS_Z4_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SS_Z4_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → sun_moon_ic |
| SS_Z4_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SS_Z6_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z6_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z6_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SS_Z6_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z6_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z6_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z6_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z6_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → sun_moon_ic |
| SS_Z8_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z8_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z8_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SS_Z8_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z8_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z8_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SS_Z8_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SS_Z8_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SY_Z10_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SY_Z10_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SY_Z10_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SY_Z10_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SY_Z10_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SY_Z10_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SY_Z10_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SY_Z10_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SY_Z4_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SY_Z4_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SY_Z4_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SY_Z4_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SY_Z4_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SY_Z4_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SY_Z4_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SY_Z4_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SY_Z6_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SY_Z6_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SY_Z6_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SY_Z6_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SY_Z6_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SY_Z6_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SY_Z6_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SY_Z6_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SY_Z8_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SY_Z8_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SY_Z8_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SY_Z8_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SY_Z8_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SY_Z8_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SY_Z8_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SY_Z8_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| a2020 | - | bubble_totani → loopI_shells+nps |

## Ranking with the consensus set (bubble_structured, erosita_unif_NS, loopI_haslam_raw, sun_moon_ic)

| rank | model | variant | n | Δ(−2lnL) vs Totani baseline | halo @ 20.8 | TS @ 20.8 | quad Σ/joint | any neg | inj min | N / S |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | gll_iem_v07 | - | 8 | -5,639 | 956 [321, 1,596] | 2.3 | -0.32 | yes | 0.99 | -247 / -61 |
| 2 | a2020 | - | 21 | -3,152 | 3,487 [2,457, 4,608] | 10.1 | -0.12 | yes | 1.00 | -260 / -172 |
| 3 | SO_Z6_R30_T150_C5 | gas1_ics1 | 9 | -1,815 | 2,306 [1,674, 2,973] | 12.9 | -0.01 | yes | 0.97 | 41 / -58 |
| 4 | SO_Z6_R20_T150_C5 | gas1_ics1 | 9 | -1,794 | 2,306 [1,667, 2,966] | 12.8 | -0.01 | yes | 0.99 | 41 / -58 |
| 5 | SO_Z4_R30_T150_C5 | gas1_ics1 | 9 | -1,738 | 2,418 [1,767, 3,067] | 13.9 | 0.08 | yes | 0.99 | 115 / 83 |
| 6 | SO_Z4_R20_T150_C5 | gas1_ics1 | 9 | -1,734 | 2,418 [1,766, 3,066] | 13.9 | 0.08 | yes | 0.99 | 115 / 83 |
| 7 | SO_Z8_R30_T150_C5 | gas1_ics1 | 9 | -1,709 | 2,250 [1,577, 2,878] | 11.8 | -0.08 | yes | 0.99 | -34 / -146 |
| 8 | SO_Z10_R30_T150_C5 | gas1_ics1 | 9 | -1,705 | 2,193 [1,549, 2,848] | 11.5 | -0.11 | yes | 0.99 | -49 / -197 |
| 9 | SO_Z8_R20_T150_C5 | gas1_ics1 | 9 | -1,684 | 2,193 [1,571, 2,872] | 11.8 | -0.10 | yes | 1.00 | -79 / -146 |
| 10 | SO_Z10_R20_T150_C5 | gas1_ics1 | 9 | -1,664 | 2,193 [1,537, 2,837] | 11.4 | -0.13 | yes | 0.99 | -79 / -197 |
| 11 | SS_Z6_R30_T150_C5 | gas1_ics1 | 9 | -1,650 | 1,856 [1,200, 2,511] | 8.1 | -0.19 | yes | 0.99 | -183 / -161 |
| 12 | SS_Z6_R20_T150_C5 | gas1_ics1 | 9 | -1,625 | 1,856 [1,194, 2,505] | 8.0 | -0.19 | yes | 0.99 | -198 / -161 |
| 13 | SS_Z4_R30_T150_C5 | gas1_ics1 | 9 | -1,607 | 1,969 [1,342, 2,651] | 9.4 | -0.06 | yes | 1.00 | -109 / -8 |
| 14 | SS_Z4_R20_T150_C5 | gas1_ics1 | 9 | -1,600 | 1,969 [1,341, 2,650] | 9.4 | -0.07 | yes | 1.00 | -124 / -8 |
| 15 | SS_Z8_R30_T150_C5 | gas1_ics1 | 9 | -1,515 | 1,744 [1,070, 2,385] | 6.9 | -0.30 | yes | 0.99 | -243 / -275 |
| 16 | SO_Z6_R30_T150_C2 | gas1_ics1 | 9 | -1,501 | 2,418 [1,784, 3,081] | 14.2 | 0.10 | yes | 0.97 | 338 / -96 |
| 17 | SS_Z10_R30_T150_C5 | gas1_ics1 | 9 | -1,498 | 1,687 [1,014, 2,332] | 6.5 | -0.38 | yes | 0.99 | -288 / -351 |
| 18 | SO_Z6_R20_T150_C2 | gas1_ics1 | 9 | -1,493 | 2,418 [1,782, 3,078] | 14.2 | 0.11 | yes | 0.98 | 338 / -83 |
| 19 | SS_Z8_R20_T150_C5 | gas1_ics1 | 9 | -1,463 | 1,687 [1,056, 2,372] | 6.8 | -0.32 | yes | 0.99 | -258 / -275 |
| 20 | SO_Z4_R30_T150_C2 | gas1_ics1 | 9 | -1,447 | 2,531 [1,883, 3,181] | 15.4 | 0.20 | yes | 0.99 | 443 / 57 |
| 21 | SO_Z4_R20_T150_C2 | gas1_ics1 | 9 | -1,444 | 2,531 [1,882, 3,180] | 15.4 | 0.20 | yes | 0.99 | 443 / 57 |
| 22 | SS_Z10_R20_T150_C5 | gas1_ics1 | 9 | -1,416 | 1,631 [988, 2,305] | 6.3 | -0.40 | yes | 0.98 | -303 / -351 |
| 23 | SO_Z8_R30_T150_C2 | gas1_ics1 | 9 | -1,406 | 2,306 [1,684, 2,982] | 13.0 | 0.03 | yes | 0.99 | 249 / -172 |
| 24 | SO_Z8_R20_T150_C2 | gas1_ics1 | 9 | -1,384 | 2,306 [1,678, 2,976] | 13.0 | 0.02 | yes | 0.98 | 234 / -185 |
| 25 | SO_Z10_R30_T150_C2 | gas1_ics1 | 9 | -1,383 | 2,306 [1,651, 2,948] | 12.7 | -0.01 | yes | 0.99 | 219 / -249 |
| 26 | SS_Z6_R30_T150_C2 | gas1_ics1 | 9 | -1,348 | 1,969 [1,298, 2,608] | 8.9 | 0.00 | yes | 0.99 | 204 / -199 |
| 27 | SO_Z10_R20_T150_C2 | gas1_ics1 | 9 | -1,347 | 2,306 [1,638, 2,936] | 12.5 | -0.02 | yes | 0.99 | 204 / -249 |
| 28 | SS_Z4_R30_T150_C2 | gas1_ics1 | 9 | -1,334 | 2,081 [1,447, 2,754] | 10.4 | 0.10 | yes | 1.00 | 248 / -46 |
| 29 | SS_Z6_R20_T150_C2 | gas1_ics1 | 9 | -1,328 | 1,969 [1,292, 2,603] | 8.9 | 0.00 | yes | 1.00 | 204 / -199 |
| 30 | SS_Z4_R20_T150_C2 | gas1_ics1 | 9 | -1,327 | 2,081 [1,446, 2,752] | 10.4 | 0.10 | yes | 1.00 | 248 / -46 |
| 31 | SS_Z8_R30_T150_C2 | gas1_ics1 | 9 | -1,266 | 1,856 [1,189, 2,503] | 7.9 | -0.11 | yes | 0.99 | 114 / -313 |
| 32 | SS_Z8_R20_T150_C2 | gas1_ics1 | 9 | -1,203 | 1,800 [1,164, 2,478] | 7.7 | -0.13 | yes | 0.99 | 84 / -313 |
| 33 | SS_Z10_R30_T150_C2 | gas1_ics1 | 9 | -1,198 | 1,800 [1,117, 2,434] | 7.3 | -0.19 | yes | 0.99 | 54 / -389 |
| 34 | SS_Z10_R20_T150_C2 | gas1_ics1 | 9 | -1,125 | 1,744 [1,089, 2,405] | 7.1 | -0.21 | yes | 0.99 | 24 / -389 |
| 35 | SO_Z8_R30_T100000_C2 | gas1_ics1 | 9 | -1,017 | 2,362 [1,726, 3,027] | 13.4 | 0.14 | yes | 0.98 | 415 / -91 |
| 36 | SO_Z6_R30_T100000_C2 | gas1_ics1 | 9 | -945 | 2,418 [1,789, 3,091] | 14.2 | 0.20 | yes | 1.00 | 475 / 10 |
| 37 | SO_Z6_R30_T100000_C5 | gas1_ics1 | 9 | -943 | 2,362 [1,736, 3,041] | 13.5 | 0.09 | yes | 1.00 | 267 / -59 |
| 38 | SO_Z6_R20_T100000_C2 | gas1_ics1 | 9 | -931 | 2,418 [1,785, 3,088] | 14.1 | 0.20 | yes | 1.00 | 475 / 10 |
| 39 | SO_Z6_R20_T100000_C5 | gas1_ics1 | 9 | -926 | 2,362 [1,733, 3,037] | 13.5 | 0.14 | yes | 1.00 | 267 / 62 |
| 40 | SO_Z8_R20_T100000_C2 | gas1_ics1 | 9 | -923 | 2,362 [1,704, 3,006] | 13.2 | 0.14 | yes | 1.00 | 400 / -78 |
| 41 | SO_Z4_R30_T100000_C2 | gas1_ics1 | 9 | -906 | 2,587 [1,939, 3,240] | 16.0 | 0.18 | yes | 0.99 | 460 / 12 |
| 42 | SO_Z10_R30_T100000_C5 | gas1_ics1 | 9 | -903 | 2,250 [1,617, 2,922] | 12.2 | 0.03 | yes | 1.00 | 192 / -116 |
| 43 | SO_Z4_R20_T100000_C2 | gas1_ics1 | 9 | -902 | 2,587 [1,938, 3,239] | 16.0 | 0.24 | yes | 0.99 | 460 / 162 |
| 44 | SO_Z10_R30_T100000_C2 | gas1_ics1 | 9 | -901 | 2,306 [1,661, 2,963] | 12.7 | 0.08 | yes | 0.98 | 355 / -169 |
| 45 | SY_Z6_R30_T150_C5 | gas1_ics1 | 9 | -900 | 1,125 [424, 1,773] | 2.7 | -0.99 | yes | 0.99 | -647 / -466 |
| 46 | SO_Z8_R30_T100000_C5 | gas1_ics1 | 9 | -900 | 2,306 [1,652, 2,956] | 12.6 | 0.07 | yes | 1.00 | 222 / -67 |
| 47 | SY_Z6_R20_T150_C5 | gas1_ics1 | 9 | -875 | 1,069 [423, 1,771] | 2.7 | -1.04 | yes | 0.97 | -647 / -466 |
| 48 | SS_Z8_R30_T100000_C5 | gas1_ics1 | 9 | -873 | 1,856 [1,212, 2,532] | 8.1 | -0.13 | yes | 1.00 | -17 / -222 |
| 49 | SS_Z8_R30_T100000_C2 | gas1_ics1 | 9 | -850 | 1,912 [1,244, 2,564] | 8.4 | -0.05 | yes | 0.99 | 131 / -220 |
| 50 | SY_Z4_R30_T150_C5 | gas1_ics1 | 9 | -840 | 1,294 [609, 1,944] | 3.7 | -0.62 | yes | 0.99 | -558 / -250 |
| 51 | SO_Z10_R20_T100000_C2 | gas1_ics1 | 9 | -836 | 2,306 [1,645, 2,948] | 12.5 | 0.08 | yes | 1.00 | 340 / -155 |
| 52 | SO_Z10_R20_T100000_C5 | gas1_ics1 | 9 | -836 | 2,250 [1,601, 2,905] | 12.0 | 0.04 | yes | 1.00 | 177 / -90 |
| 53 | SY_Z4_R20_T150_C5 | gas1_ics1 | 9 | -833 | 1,294 [608, 1,944] | 3.7 | -0.62 | yes | 0.98 | -558 / -250 |
| 54 | SS_Z4_R30_T100000_C5 | gas1_ics1 | 9 | -803 | 2,137 [1,482, 2,793] | 10.7 | 0.13 | yes | 1.00 | 133 / 136 |
| 55 | SO_Z8_R20_T100000_C5 | gas1_ics1 | 9 | -796 | 2,306 [1,633, 2,938] | 12.4 | 0.05 | yes | 1.00 | 192 / -67 |
| 56 | SS_Z4_R30_T100000_C2 | gas1_ics1 | 9 | -794 | 2,193 [1,528, 2,838] | 11.2 | 0.19 | yes | 1.00 | 311 / 98 |
| 57 | SS_Z8_R20_T100000_C2 | gas1_ics1 | 9 | -792 | 1,912 [1,228, 2,548] | 8.2 | -0.06 | yes | 0.99 | 71 / -179 |
| 58 | SS_Z4_R20_T100000_C2 | gas1_ics1 | 9 | -783 | 2,193 [1,526, 2,836] | 11.2 | 0.18 | yes | 1.00 | 296 / 98 |
| 59 | SS_Z6_R30_T100000_C5 | gas1_ics1 | 9 | -768 | 1,969 [1,292, 2,609] | 8.8 | -0.03 | yes | 1.00 | 28 / -96 |
| 60 | SS_Z6_R30_T100000_C2 | gas1_ics1 | 9 | -764 | 1,969 [1,328, 2,644] | 9.2 | 0.06 | yes | 1.00 | 206 / -79 |
| 61 | SY_Z8_R30_T150_C5 | gas1_ics1 | 9 | -746 | 956 [271, 1,630] | 2.0 | -1.42 | yes | 0.99 | -722 / -631 |
| 62 | SS_Z8_R20_T100000_C5 | gas1_ics1 | 9 | -743 | 1,856 [1,182, 2,502] | 7.8 | -0.11 | yes | 0.99 | -32 / -168 |
| 63 | SY_Z10_R30_T150_C5 | gas1_ics1 | 9 | -739 | 900 [205, 1,571] | 1.7 | -1.63 | yes | 0.99 | -737 / -733 |
| 64 | SS_Z6_R20_T100000_C5 | gas1_ics1 | 9 | -736 | 1,969 [1,284, 2,601] | 8.8 | -0.01 | yes | 1.00 | 28 / -55 |
| 65 | SS_Z6_R20_T100000_C2 | gas1_ics1 | 9 | -733 | 1,969 [1,320, 2,636] | 9.1 | 0.06 | yes | 1.00 | 206 / -93 |
| 66 | SS_Z10_R30_T100000_C2 | gas1_ics1 | 9 | -719 | 1,800 [1,156, 2,479] | 7.6 | -0.07 | yes | 0.99 | 176 / -295 |
| 67 | SS_Z10_R30_T100000_C5 | gas1_ics1 | 9 | -711 | 1,800 [1,124, 2,448] | 7.3 | -0.18 | yes | 0.99 | -47 / -271 |
| 68 | SO_Z4_R30_T100000_C5 | gas1_ics1 | 9 | -688 | 2,475 [1,841, 3,146] | 14.7 | 0.22 | yes | 0.99 | 342 / 214 |
| 69 | SY_Z8_R20_T150_C5 | gas1_ics1 | 9 | -685 | 956 [262, 1,620] | 1.9 | -1.42 | yes | 0.97 | -737 / -618 |
| 70 | SO_Z4_R20_T100000_C5 | gas1_ics1 | 9 | -681 | 2,475 [1,840, 3,145] | 14.7 | 0.14 | yes | 0.99 | 342 / 10 |
| 71 | SY_Z10_R20_T150_C5 | gas1_ics1 | 9 | -651 | 900 [191, 1,557] | 1.6 | -1.64 | yes | 0.98 | -752 / -720 |
| 72 | SY_Z6_R30_T150_C2 | gas1_ics1 | 9 | -610 | 1,181 [519, 1,868] | 3.1 | -0.73 | yes | 0.99 | -350 / -517 |
| 73 | SS_Z10_R20_T100000_C5 | gas1_ics1 | 9 | -601 | 1,744 [1,097, 2,420] | 7.1 | -0.18 | yes | 0.99 | -77 / -244 |
| 74 | SS_Z10_R20_T100000_C2 | gas1_ics1 | 9 | -591 | 1,800 [1,124, 2,447] | 7.3 | -0.08 | yes | 0.99 | 146 / -296 |
| 75 | SS_Z4_R20_T100000_C5 | gas1_ics1 | 9 | -588 | 2,081 [1,445, 2,757] | 10.3 | 0.12 | yes | 1.00 | 117 / 124 |
| 76 | SY_Z6_R20_T150_C2 | gas1_ics1 | 9 | -587 | 1,181 [518, 1,866] | 3.1 | -0.75 | yes | 0.99 | -365 / -517 |
| 77 | SY_Z4_R30_T150_C2 | gas1_ics1 | 9 | -572 | 1,350 [703, 2,038] | 4.2 | -0.37 | yes | 1.00 | -216 / -289 |
| 78 | SY_Z4_R20_T150_C2 | gas1_ics1 | 9 | -566 | 1,350 [703, 2,038] | 4.2 | -0.37 | yes | 0.98 | -216 / -289 |
| 79 | SY_Z8_R30_T150_C2 | gas1_ics1 | 9 | -523 | 1,069 [391, 1,750] | 2.5 | -1.07 | yes | 0.97 | -469 / -669 |
| 80 | SY_Z10_R30_T150_C2 | gas1_ics1 | 9 | -463 | 1,012 [314, 1,680] | 2.1 | -1.27 | yes | 0.98 | -514 / -771 |
| 81 | SY_Z8_R30_T100000_C5 | gas1_ics1 | 9 | -450 | 1,125 [462, 1,830] | 2.8 | -1.15 | yes | 0.99 | -511 / -780 |
| 82 | SY_Z8_R20_T150_C2 | gas1_ics1 | 9 | -428 | 1,069 [366, 1,725] | 2.4 | -1.09 | yes | 0.98 | -499 / -669 |
| 83 | SY_Z8_R30_T100000_C2 | gas1_ics1 | 9 | -378 | 1,181 [484, 1,854] | 2.9 | -0.85 | yes | 0.97 | -318 / -682 |
| 84 | SY_Z10_R20_T150_C2 | gas1_ics1 | 9 | -374 | 956 [292, 1,659] | 2.0 | -1.38 | yes | 0.98 | -544 / -771 |
| 85 | SY_Z8_R20_T100000_C2 | gas1_ics1 | 9 | -317 | 1,181 [476, 1,845] | 2.9 | -0.62 | yes | 0.97 | -333 / -397 |
| 86 | SY_Z6_R30_T100000_C5 | gas1_ics1 | 9 | -300 | 1,237 [545, 1,903] | 3.3 | -0.52 | yes | 1.00 | -451 / -194 |
| 87 | SL_Z6_R30_T150_C5 | gas1_ics1 | 9 | -293 | 675 [7, 1,386] | 1.0 | -2.46 | yes | 0.98 | -977 / -683 |
| 88 | SY_Z6_R20_T100000_C5 | gas1_ics1 | 9 | -269 | 1,237 [541, 1,899] | 3.2 | -0.62 | yes | 0.99 | -466 / -303 |
| 89 | SL_Z6_R20_T150_C5 | gas1_ics1 | 9 | -262 | 675 [4, 1,383] | 1.0 | -2.46 | yes | 0.96 | -977 / -683 |
| 90 | SY_Z10_R30_T100000_C5 | gas1_ics1 | 9 | -261 | 1,069 [362, 1,738] | 2.3 | -1.17 | yes | 0.97 | -556 / -691 |
| 91 | SY_Z4_R30_T100000_C2 | gas1_ics1 | 9 | -253 | 1,462 [806, 2,147] | 4.9 | -0.11 | yes | 0.98 | -49 / -116 |
| 92 | SY_Z8_R20_T100000_C5 | gas1_ics1 | 9 | -249 | 1,125 [424, 1,791] | 2.6 | -1.00 | yes | 1.00 | -511 / -618 |
| 93 | SY_Z4_R20_T100000_C2 | gas1_ics1 | 9 | -244 | 1,462 [805, 2,146] | 4.9 | -0.26 | yes | 0.98 | -228 / -156 |
| 94 | SY_Z6_R30_T100000_C2 | gas1_ics1 | 9 | -243 | 1,237 [570, 1,928] | 3.4 | -0.56 | yes | 1.00 | -244 / -450 |
| 95 | SY_Z10_R30_T100000_C2 | gas1_ics1 | 9 | -222 | 1,069 [386, 1,763] | 2.4 | -1.22 | yes | 0.99 | -408 / -891 |
| 96 | SY_Z6_R20_T100000_C2 | gas1_ics1 | 9 | -212 | 1,237 [566, 1,924] | 3.4 | -0.56 | yes | 1.00 | -244 / -450 |
| 97 | SL_Z10_R30_T150_C5 | gas1_ics1 | 9 | -178 | 506 [-207, 1,196] | 0.5 | -4.23 | yes | 0.97 | -1,156 / -987 |
| 98 | SL_Z4_R30_T150_C5 | gas1_ics1 | 9 | -169 | 900 [196, 1,556] | 1.7 | -1.37 | yes | 0.97 | -827 / -403 |
| 99 | SL_Z4_R20_T150_C5 | gas1_ics1 | 9 | -163 | 900 [195, 1,556] | 1.7 | -1.37 | yes | 0.97 | -827 / -403 |
| 100 | SL_Z8_R30_T150_C5 | gas1_ics1 | 9 | -148 | 562 [-157, 1,238] | 0.6 | -3.55 | yes | 0.96 | -1,126 / -873 |
| 101 | SY_Z10_R20_T100000_C5 | gas1_ics1 | 9 | -138 | 1,012 [340, 1,714] | 2.2 | -1.33 | yes | 0.99 | -482 / -869 |
| 102 | SL_Z8_R30_T100000_C5 | gas1_ics1 | 9 | -129 | 787 [58, 1,464] | 1.2 | -2.22 | yes | 0.96 | -796 / -952 |
| 103 | SL_Z10_R20_T150_C5 | gas1_ics1 | 9 | -94 | 506 [-216, 1,186] | 0.5 | -4.21 | yes | 0.98 | -1,156 / -974 |
| 104 | SL_Z8_R20_T150_C5 | gas1_ics1 | 9 | -92 | 506 [-163, 1,231] | 0.6 | -3.95 | yes | 0.97 | -1,126 / -873 |
| 105 | SY_Z10_R20_T100000_C2 | gas1_ics1 | 9 | -81 | 1,069 [359, 1,735] | 2.3 | -1.18 | yes | 0.97 | -438 / -824 |
| 106 | SL_Z8_R30_T100000_C2 | gas1_ics1 | 9 | -25 | 787 [86, 1,493] | 1.3 | -1.99 | yes | 0.98 | -723 / -840 |
| 107 | SY_Z4_R30_T100000_C5 | gas1_ics1 | 9 | -11 | 1,406 [725, 2,066] | 4.3 | -0.38 | yes | 0.98 | -332 / -199 |
| 108 | SY_Z4_R20_T100000_C5 | gas1_ics1 | 9 | -0 | 1,406 [724, 2,065] | 4.3 | -0.33 | yes | 0.98 | -332 / -131 |
| 109 | SL_Z6_R30_T150_C2 | gas1_ics1 | 9 | +0 | 787 [106, 1,487] | 1.3 | -1.89 | yes | 0.96 | -754 / -733 |
| 110 | SL_Z6_R20_T150_C2 | gas1_ics1 | 9 | +20 | 787 [106, 1,486] | 1.3 | -1.89 | yes | 0.98 | -754 / -733 |
| 111 | SL_Z8_R20_T100000_C2 | gas1_ics1 | 9 | +31 | 787 [82, 1,488] | 1.2 | -2.06 | yes | 0.96 | -738 / -881 |
| 112 | SL_Z10_R30_T100000_C5 | gas1_ics1 | 9 | +59 | 675 [-40, 1,374] | 0.9 | -2.69 | yes | 0.99 | -871 / -944 |
| 113 | SL_Z6_R30_T100000_C5 | gas1_ics1 | 9 | +64 | 844 [134, 1,525] | 1.4 | -1.65 | yes | 0.99 | -751 / -639 |
| 114 | SL_Z8_R30_T150_C2 | gas1_ics1 | 9 | +75 | 675 [-27, 1,368] | 0.9 | -2.66 | yes | 0.98 | -873 / -924 |
| 115 | SL_Z6_R20_T100000_C5 | gas1_ics1 | 9 | +92 | 844 [132, 1,523] | 1.4 | -1.52 | yes | 0.99 | -751 / -531 |
| 116 | SL_Z10_R30_T150_C2 | gas1_ics1 | 9 | +101 | 619 [-96, 1,309] | 0.7 | -3.21 | yes | 0.98 | -948 / -1,038 |
| 117 | SL_Z4_R30_T150_C2 | gas1_ics1 | 9 | +111 | 956 [292, 1,652] | 2.0 | -1.08 | yes | 0.99 | -589 / -441 |
| 118 | SL_Z4_R20_T150_C2 | gas1_ics1 | 9 | +116 | 956 [291, 1,651] | 2.0 | -1.08 | yes | 0.98 | -589 / -441 |
| 119 | SL_Z8_R20_T100000_C5 | gas1_ics1 | 9 | +120 | 731 [13, 1,417] | 1.0 | -2.28 | yes | 0.96 | -826 / -843 |
| 120 | SL_Z10_R30_T100000_C2 | gas1_ics1 | 9 | +139 | 675 [-14, 1,404] | 1.0 | -3.02 | yes | 0.98 | -813 / -1,226 |
| 121 | SL_Z6_R30_T100000_C2 | gas1_ics1 | 9 | +165 | 844 [160, 1,553] | 1.5 | -1.59 | yes | 0.99 | -648 / -691 |
| 122 | SL_Z8_R20_T150_C2 | gas1_ics1 | 9 | +169 | 619 [-51, 1,344] | 0.9 | -2.95 | yes | 0.98 | -903 / -924 |
| 123 | SL_Z10_R20_T100000_C5 | gas1_ics1 | 9 | +178 | 675 [-58, 1,355] | 0.8 | -2.71 | yes | 0.99 | -901 / -931 |
| 124 | SL_Z10_R20_T150_C2 | gas1_ics1 | 9 | +188 | 562 [-113, 1,292] | 0.7 | -3.54 | yes | 0.97 | -963 / -1,025 |
| 125 | SL_Z6_R20_T100000_C2 | gas1_ics1 | 9 | +191 | 844 [159, 1,551] | 1.5 | -1.70 | yes | 0.99 | -648 / -787 |
| 126 | SL_Z4_R30_T100000_C2 | gas1_ics1 | 9 | +236 | 1,069 [399, 1,767] | 2.5 | -0.72 | yes | 1.00 | -438 / -335 |
| 127 | SL_Z4_R20_T100000_C2 | gas1_ics1 | 9 | +239 | 1,069 [398, 1,767] | 2.5 | -0.72 | yes | 1.00 | -438 / -335 |
| 128 | SL_Z10_R20_T100000_C2 | gas1_ics1 | 9 | +270 | 675 [-34, 1,381] | 0.9 | -2.72 | yes | 0.98 | -843 / -995 |
| 129 | SL_Z4_R30_T100000_C5 | gas1_ics1 | 9 | +466 | 1,012 [314, 1,682] | 2.1 | -0.94 | yes | 1.00 | -646 / -309 |
| 130 | SL_Z4_R20_T100000_C5 | gas1_ics1 | 9 | +474 | 1,012 [313, 1,681] | 2.1 | -0.93 | yes | 1.00 | -631 / -309 |
