# Stage C — 2026-09-21 (gas5_ics3_no_residual)

Consensus selection across 129 model rows (gas5_ics3) by the **fraction-of-rows rule**: a candidate is eligible if it clears its own BIC threshold (p·ln 1080 = 63 per template) in ≥ 80% of rows; the eligible candidate that is the largest improvement in the most rows is selected. Median Δ is shown for reference only. Selected: **bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic**

| step | chosen | clears BIC in | best among eligible in | median Δ (ref.) | runner-up (clears / best) |
|---|---|---|---|---|---|
| 1 | bubble_flat | 100% | 65% | 15,305 | bubble_structured (100% / 35%) |
| 2 | erosita_unif_outer_NS | 100% | 53% | 10,157 | erosita_unif_NS (100% / 37%) |
| 3 | loopI_haslam_raw | 100% | 95% | 890 | loopI_shells+nps (100% / 5%) |
| 4 | sun_moon_ic | 100% | 100% | 348 | loopIV (98% / 0%) |
| 5 | (stop) loopIV | 52% | 0% | 64 | loopII (0% / 0%) |

## Each model's own path

| model | variant | path |
|---|---|---|
| SL_Z10_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SL_Z10_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z10_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z10_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z10_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z10_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z10_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z10_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z4_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SL_Z4_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z4_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z4_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z4_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SL_Z4_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z4_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z4_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SL_Z6_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z6_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z6_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z6_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z6_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z6_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z6_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z6_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z8_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z8_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z8_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z8_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z8_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z8_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z8_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z8_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SO_Z10_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z10_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z10_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SO_Z10_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z10_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SO_Z10_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z10_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SO_Z10_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z4_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SO_Z4_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SO_Z4_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SO_Z4_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SO_Z4_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z4_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SO_Z4_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z4_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SO_Z6_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z6_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z6_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SO_Z6_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z6_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z6_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z6_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SO_Z6_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z8_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z8_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z8_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SO_Z8_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SO_Z8_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z8_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z8_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SO_Z8_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SS_Z10_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SS_Z10_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SS_Z10_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z10_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z10_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SS_Z10_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SS_Z10_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z10_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SS_Z4_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z4_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SS_Z4_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z4_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SS_Z4_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SS_Z4_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SS_Z4_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z4_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → erosita_unif_NS → sun_moon_ic |
| SS_Z6_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SS_Z6_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SS_Z6_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z6_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SS_Z6_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SS_Z6_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SS_Z6_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z6_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z8_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SS_Z8_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SS_Z8_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z8_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z8_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SS_Z8_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SS_Z8_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z8_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z10_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SY_Z10_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SY_Z10_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z10_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z10_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z10_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z10_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z10_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z4_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SY_Z4_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SY_Z4_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z4_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SY_Z4_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SY_Z4_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SY_Z4_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z4_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SY_Z6_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z6_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SY_Z6_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z6_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z6_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SY_Z6_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SY_Z6_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z6_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SY_Z8_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SY_Z8_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SY_Z8_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z8_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z8_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z8_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SY_Z8_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z8_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| a2020 | - | bubble_totani → loopI_shells+nps |

## Ranking with the consensus set (bubble_flat, erosita_unif_outer_NS, loopI_haslam_raw, sun_moon_ic)

| rank | model | variant | n | Δ(−2lnL) vs Totani baseline | halo @ 20.8 | TS @ 20.8 | quad Σ/joint | any neg | inj min | N / S |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | gll_iem_v07 | - | 8 | - | 2,025 [1,442, 2,587] | 12.5 | 0.84 | yes | 0.97 | 652 / 1,045 |
| 2 | a2020 | - | 21 | - | 7,424 [6,878, 7,967] | 49.8 | 0.92 | no | 1.00 | 3,147 / 3,711 |
| 3 | SS_Z10_R30_T150_C5 | gas5_ics3 | 15 | - | 5,962 [5,238, 6,726] | 11.3 | 0.61 | no | 0.68 | 2,075 / 1,580 |
| 4 | SY_Z10_R30_T150_C5 | gas5_ics3 | 15 | - | 5,793 [4,955, 6,679] | 8.8 | 0.61 | no | 0.68 | 2,031 / 1,475 |
| 5 | SL_Z10_R30_T150_C5 | gas5_ics3 | 15 | - | 5,456 [4,088, 6,815] | 8.1 | 0.57 | no | 0.68 | 1,823 / 1,306 |
| 6 | SY_Z10_R20_T150_C5 | gas5_ics3 | 15 | - | 5,793 [4,963, 6,677] | 9.1 | 0.59 | no | 0.68 | 1,942 / 1,503 |
| 7 | SS_Z10_R20_T150_C5 | gas5_ics3 | 15 | - | 5,962 [4,514, 6,717] | 11.5 | 0.62 | no | 0.68 | 2,090 / 1,597 |
| 8 | SL_Z10_R20_T150_C5 | gas5_ics3 | 15 | - | 5,849 [5,018, 6,778] | 8.5 | 0.56 | no | 0.68 | 1,942 / 1,346 |
| 9 | SO_Z10_R30_T150_C5 | gas5_ics3 | 15 | - | 6,412 [5,782, 7,000] | 11.4 | 0.56 | no | 0.68 | 2,075 / 1,513 |
| 10 | SS_Z8_R30_T150_C5 | gas5_ics3 | 15 | - | 6,018 [5,298, 6,744] | 14.3 | 0.69 | no | 0.67 | 2,271 / 1,863 |
| 11 | SY_Z8_R30_T150_C5 | gas5_ics3 | 15 | - | 5,849 [5,028, 6,703] | 11.5 | 0.67 | no | 0.67 | 2,167 / 1,772 |
| 12 | SS_Z6_R30_T150_C5 | gas5_ics3 | 15 | - | 6,131 [5,427, 6,810] | 20.4 | 0.80 | no | 0.67 | 2,618 / 2,291 |
| 13 | SY_Z6_R30_T150_C5 | gas5_ics3 | 15 | - | 6,018 [5,194, 6,792] | 17.3 | 0.76 | no | 0.67 | 2,379 / 2,213 |
| 14 | SO_Z10_R20_T150_C5 | gas5_ics3 | 15 | - | 6,412 [5,780, 6,990] | 11.6 | 0.57 | no | 0.68 | 2,075 / 1,554 |
| 15 | SY_Z6_R20_T150_C5 | gas5_ics3 | 15 | - | 6,018 [5,198, 6,792] | 17.4 | 0.76 | no | 0.67 | 2,379 / 2,213 |
| 16 | SS_Z6_R20_T150_C5 | gas5_ics3 | 15 | - | 6,131 [5,429, 6,809] | 20.5 | 0.80 | no | 0.67 | 2,588 / 2,318 |
| 17 | SL_Z8_R30_T150_C5 | gas5_ics3 | 15 | - | 5,906 [5,063, 6,774] | 10.5 | 0.62 | no | 0.67 | 1,989 / 1,678 |
| 18 | SS_Z8_R20_T150_C5 | gas5_ics3 | 15 | - | 6,018 [5,302, 6,742] | 14.4 | 0.69 | no | 0.67 | 2,286 / 1,877 |
| 19 | SL_Z6_R30_T150_C5 | gas5_ics3 | 15 | - | 6,018 [5,221, 6,797] | 16.0 | 0.72 | no | 0.67 | 2,170 / 2,148 |
| 20 | SY_Z8_R20_T150_C5 | gas5_ics3 | 15 | - | 5,849 [5,033, 6,702] | 11.7 | 0.67 | no | 0.67 | 2,182 / 1,747 |
| 21 | SL_Z8_R20_T150_C5 | gas5_ics3 | 15 | - | 5,906 [5,067, 6,771] | 10.6 | 0.62 | no | 0.67 | 1,989 / 1,692 |
| 22 | SL_Z6_R20_T150_C5 | gas5_ics3 | 15 | - | 6,018 [5,222, 6,793] | 16.1 | 0.72 | no | 0.67 | 2,185 / 2,161 |
| 23 | SO_Z8_R30_T150_C5 | gas5_ics3 | 15 | - | 6,412 [5,801, 6,978] | 13.9 | 0.60 | no | 0.67 | 2,060 / 1,769 |
| 24 | SO_Z8_R20_T150_C5 | gas5_ics3 | 15 | - | 6,412 [5,801, 6,973] | 14.1 | 0.60 | no | 0.67 | 2,060 / 1,783 |
| 25 | SO_Z6_R30_T150_C5 | gas5_ics3 | 15 | - | 6,468 [5,647, 7,000] | 20.0 | 0.71 | no | 0.67 | 2,347 / 2,237 |
| 26 | SO_Z6_R20_T150_C5 | gas5_ics3 | 15 | - | 6,243 [5,597, 7,024] | 20.0 | 0.73 | no | 0.67 | 2,332 / 2,250 |
| 27 | SS_Z10_R30_T150_C2 | gas5_ics3 | 15 | - | 6,074 [4,362, 6,895] | 10.8 | 0.63 | no | 0.68 | 2,225 / 1,621 |
| 28 | SY_Z10_R30_T150_C2 | gas5_ics3 | 15 | - | 6,018 [5,156, 6,892] | 8.4 | 0.60 | no | 0.68 | 2,106 / 1,503 |
| 29 | SS_Z8_R30_T150_C2 | gas5_ics3 | 15 | - | 6,074 [4,670, 6,930] | 13.4 | 0.69 | no | 0.67 | 2,331 / 1,850 |
| 30 | SL_Z10_R30_T150_C2 | gas5_ics3 | 15 | - | 6,187 [5,227, 6,961] | 7.8 | 0.52 | no | 0.68 | 1,913 / 1,333 |
| 31 | SY_Z8_R30_T150_C2 | gas5_ics3 | 15 | - | 6,074 [5,227, 6,918] | 10.9 | 0.64 | no | 0.67 | 2,183 / 1,734 |
| 32 | SS_Z10_R20_T150_C2 | gas5_ics3 | 15 | - | 5,793 [4,400, 6,882] | 10.9 | 0.66 | no | 0.68 | 2,150 / 1,662 |
| 33 | SY_Z10_R20_T150_C2 | gas5_ics3 | 15 | - | 5,849 [4,064, 6,877] | 8.6 | 0.61 | no | 0.68 | 2,047 / 1,506 |
| 34 | SO_Z10_R30_T150_C2 | gas5_ics3 | 15 | - | 6,412 [5,146, 7,150] | 10.9 | 0.57 | no | 0.68 | 2,134 / 1,503 |
| 35 | SL_Z10_R20_T150_C2 | gas5_ics3 | 15 | - | 6,187 [5,226, 6,950] | 8.0 | 0.53 | no | 0.68 | 1,928 / 1,374 |
| 36 | SS_Z8_R20_T150_C2 | gas5_ics3 | 15 | - | 6,187 [5,462, 6,917] | 13.7 | 0.68 | no | 0.67 | 2,301 / 1,891 |
| 37 | SL_Z8_R30_T150_C2 | gas5_ics3 | 15 | - | 5,962 [5,260, 6,959] | 10.0 | 0.62 | no | 0.67 | 2,079 / 1,616 |
| 38 | SO_Z10_R20_T150_C2 | gas5_ics3 | 15 | - | 6,412 [5,139, 7,136] | 11.1 | 0.57 | no | 0.68 | 2,134 / 1,544 |
| 39 | SS_Z4_R30_T150_C5 | gas5_ics3 | 15 | - | 6,187 [5,505, 6,797] | 37.6 | 0.83 | no | 0.66 | 2,335 / 2,794 |
| 40 | SS_Z4_R20_T150_C5 | gas5_ics3 | 15 | - | 6,131 [5,506, 6,796] | 37.7 | 0.84 | no | 0.66 | 2,350 / 2,794 |
| 41 | SY_Z8_R20_T150_C2 | gas5_ics3 | 15 | - | 6,074 [5,219, 6,902] | 11.1 | 0.68 | no | 0.67 | 2,302 / 1,839 |
| 42 | SS_Z6_R30_T150_C2 | gas5_ics3 | 15 | - | 6,243 [5,570, 6,968] | 19.6 | 0.80 | no | 0.67 | 2,633 / 2,332 |
| 43 | SY_Z6_R30_T150_C2 | gas5_ics3 | 15 | - | 6,074 [5,357, 6,985] | 16.6 | 0.77 | no | 0.67 | 2,394 / 2,280 |
| 44 | SY_Z4_R30_T150_C5 | gas5_ics3 | 15 | - | 6,074 [5,354, 6,826] | 33.5 | 0.82 | no | 0.66 | 2,290 / 2,703 |
| 45 | SS_Z6_R20_T150_C2 | gas5_ics3 | 15 | - | 6,243 [5,571, 6,967] | 19.7 | 0.80 | no | 0.67 | 2,663 / 2,332 |
| 46 | SY_Z4_R20_T150_C5 | gas5_ics3 | 15 | - | 6,074 [5,354, 6,825] | 33.5 | 0.82 | no | 0.66 | 2,290 / 2,703 |
| 47 | SY_Z6_R20_T150_C2 | gas5_ics3 | 15 | - | 6,243 [5,363, 6,982] | 16.7 | 0.74 | no | 0.67 | 2,379 / 2,253 |
| 48 | SL_Z8_R20_T150_C2 | gas5_ics3 | 15 | - | 6,131 [5,270, 6,932] | 10.2 | 0.62 | no | 0.67 | 2,064 / 1,733 |
| 49 | SO_Z8_R30_T150_C2 | gas5_ics3 | 15 | - | 6,524 [5,932, 7,122] | 13.5 | 0.61 | no | 0.67 | 2,135 / 1,836 |
| 50 | SL_Z6_R30_T150_C2 | gas5_ics3 | 15 | - | 6,243 [5,408, 6,937] | 15.4 | 0.70 | no | 0.67 | 2,245 / 2,133 |
| 51 | SL_Z6_R20_T150_C2 | gas5_ics3 | 15 | - | 6,243 [5,410, 6,934] | 15.5 | 0.71 | no | 0.67 | 2,215 / 2,188 |
| 52 | SO_Z8_R20_T150_C2 | gas5_ics3 | 15 | - | 6,524 [5,931, 7,115] | 13.6 | 0.61 | no | 0.67 | 2,135 / 1,837 |
| 53 | SL_Z4_R30_T150_C5 | gas5_ics3 | 15 | - | 6,131 [5,374, 6,718] | 31.6 | 0.78 | no | 0.66 | 2,096 / 2,677 |
| 54 | SL_Z4_R20_T150_C5 | gas5_ics3 | 15 | - | 6,131 [5,374, 6,717] | 31.6 | 0.79 | no | 0.66 | 2,141 / 2,677 |
| 55 | SO_Z6_R30_T150_C2 | gas5_ics3 | 15 | - | 6,524 [5,986, 7,115] | 19.3 | 0.72 | no | 0.67 | 2,407 / 2,277 |
| 56 | SO_Z6_R20_T150_C2 | gas5_ics3 | 15 | - | 6,468 [5,669, 7,130] | 19.4 | 0.72 | no | 0.67 | 2,392 / 2,277 |
| 57 | SO_Z4_R30_T150_C5 | gas5_ics3 | 15 | - | 6,243 [5,893, 7,127] | 38.3 | 0.83 | no | 0.66 | 2,335 / 2,861 |
| 58 | SO_Z4_R20_T150_C5 | gas5_ics3 | 15 | - | 6,749 [5,902, 7,116] | 38.3 | 0.77 | no | 0.66 | 2,364 / 2,861 |
| 59 | SS_Z4_R30_T150_C2 | gas5_ics3 | 15 | - | 6,243 [5,635, 6,942] | 36.9 | 0.83 | no | 0.66 | 2,364 / 2,793 |
| 60 | SS_Z4_R20_T150_C2 | gas5_ics3 | 15 | - | 6,356 [5,633, 6,943] | 36.9 | 0.81 | no | 0.66 | 2,364 / 2,780 |
| 61 | SY_Z4_R30_T150_C2 | gas5_ics3 | 15 | - | 6,243 [5,510, 6,976] | 32.9 | 0.80 | no | 0.66 | 2,275 / 2,702 |
| 62 | SY_Z4_R20_T150_C2 | gas5_ics3 | 15 | - | 6,243 [5,510, 6,975] | 32.9 | 0.79 | no | 0.66 | 2,245 / 2,702 |
| 63 | SL_Z4_R20_T150_C2 | gas5_ics3 | 15 | - | 6,299 [5,803, 6,829] | 31.1 | 0.75 | no | 0.66 | 2,081 / 2,664 |
| 64 | SL_Z4_R30_T150_C2 | gas5_ics3 | 15 | - | 6,299 [5,546, 6,830] | 31.1 | 0.76 | no | 0.66 | 2,111 / 2,664 |
| 65 | SY_Z8_R30_T100000_C5 | gas5_ics3 | 15 | - | 5,568 [4,742, 6,431] | 12.9 | 0.68 | no | 0.67 | 1,972 / 1,804 |
| 66 | SL_Z8_R30_T100000_C5 | gas5_ics3 | 15 | - | 5,624 [4,728, 6,501] | 10.8 | 0.56 | no | 0.68 | 1,420 / 1,712 |
| 67 | SO_Z4_R30_T150_C2 | gas5_ics3 | 15 | - | 6,356 [6,005, 7,236] | 37.7 | 0.81 | no | 0.66 | 2,305 / 2,848 |
| 68 | SO_Z4_R20_T150_C2 | gas5_ics3 | 15 | - | 6,356 [6,005, 7,235] | 37.7 | 0.82 | no | 0.66 | 2,379 / 2,848 |
| 69 | SS_Z8_R30_T100000_C5 | gas5_ics3 | 15 | - | 5,849 [5,306, 6,443] | 16.7 | 0.74 | no | 0.67 | 2,405 / 1,896 |
| 70 | SY_Z8_R30_T100000_C2 | gas5_ics3 | 15 | - | 5,737 [4,894, 6,628] | 12.1 | 0.66 | no | 0.67 | 1,882 / 1,887 |
| 71 | SS_Z8_R30_T100000_C2 | gas5_ics3 | 15 | - | 5,962 [5,341, 6,690] | 15.7 | 0.73 | no | 0.67 | 2,360 / 1,978 |
| 72 | SL_Z8_R30_T100000_C2 | gas5_ics3 | 15 | - | 5,849 [4,928, 6,722] | 10.2 | 0.56 | no | 0.68 | 1,480 / 1,781 |
| 73 | SY_Z8_R20_T100000_C2 | gas5_ics3 | 15 | - | 5,793 [4,907, 6,632] | 12.3 | 0.64 | no | 0.67 | 1,793 / 1,900 |
| 74 | SY_Z10_R30_T100000_C5 | gas5_ics3 | 15 | - | 5,456 [3,994, 6,359] | 10.4 | 0.63 | no | 0.68 | 1,776 / 1,657 |
| 75 | SL_Z8_R20_T100000_C2 | gas5_ics3 | 15 | - | 5,793 [4,939, 6,724] | 10.4 | 0.56 | no | 0.68 | 1,436 / 1,808 |
| 76 | SS_Z8_R20_T100000_C2 | gas5_ics3 | 15 | - | 5,962 [5,345, 6,691] | 15.9 | 0.73 | no | 0.67 | 2,360 / 1,993 |
| 77 | SL_Z10_R30_T100000_C5 | gas5_ics3 | 15 | - | 5,512 [4,623, 6,442] | 8.5 | 0.46 | no | 0.68 | 1,045 / 1,488 |
| 78 | SS_Z10_R30_T100000_C5 | gas5_ics3 | 15 | - | 5,793 [5,222, 6,311] | 13.8 | 0.64 | no | 0.67 | 1,941 / 1,762 |
| 79 | SS_Z8_R20_T100000_C5 | gas5_ics3 | 15 | - | 5,849 [5,304, 6,378] | 16.9 | 0.74 | no | 0.67 | 2,345 / 1,964 |
| 80 | SY_Z8_R20_T100000_C5 | gas5_ics3 | 15 | - | 5,568 [4,741, 6,422] | 13.2 | 0.63 | no | 0.67 | 1,629 / 1,886 |
| 81 | SS_Z10_R30_T100000_C2 | gas5_ics3 | 15 | - | 5,849 [5,256, 6,638] | 12.9 | 0.67 | no | 0.67 | 2,104 / 1,805 |
| 82 | SO_Z8_R30_T100000_C2 | gas5_ics3 | 15 | - | 6,356 [5,810, 6,967] | 15.6 | 0.65 | no | 0.67 | 2,194 / 1,925 |
| 83 | SY_Z10_R30_T100000_C2 | gas5_ics3 | 15 | - | 5,399 [3,819, 6,560] | 9.6 | 0.59 | no | 0.68 | 1,478 / 1,713 |
| 84 | SL_Z10_R30_T100000_C2 | gas5_ics3 | 15 | - | 5,624 [4,823, 6,680] | 7.9 | 0.50 | no | 0.68 | 1,284 / 1,556 |
| 85 | SY_Z6_R30_T100000_C5 | gas5_ics3 | 15 | - | 5,737 [4,902, 6,535] | 19.3 | 0.76 | no | 0.67 | 2,200 / 2,178 |
| 86 | SL_Z8_R20_T100000_C5 | gas5_ics3 | 15 | - | 5,624 [4,720, 6,485] | 11.1 | 0.56 | no | 0.67 | 1,361 / 1,807 |
| 87 | SY_Z10_R20_T100000_C5 | gas5_ics3 | 15 | - | 5,512 [4,650, 6,364] | 10.7 | 0.67 | no | 0.68 | 1,955 / 1,712 |
| 88 | SS_Z6_R30_T100000_C5 | gas5_ics3 | 15 | - | 5,962 [5,413, 6,473] | 23.5 | 0.85 | no | 0.67 | 2,737 / 2,324 |
| 89 | SL_Z10_R20_T100000_C5 | gas5_ics3 | 15 | - | 5,568 [4,638, 6,445] | 8.8 | 0.48 | no | 0.68 | 1,120 / 1,555 |
| 90 | SY_Z6_R20_T100000_C5 | gas5_ics3 | 15 | - | 5,737 [4,907, 6,537] | 19.4 | 0.82 | no | 0.67 | 2,439 / 2,287 |
| 91 | SO_Z10_R30_T100000_C5 | gas5_ics3 | 15 | - | 6,243 [5,709, 6,743] | 13.8 | 0.61 | no | 0.67 | 2,104 / 1,684 |
| 92 | SO_Z10_R30_T100000_C2 | gas5_ics3 | 15 | - | 6,299 [5,752, 6,947] | 13.0 | 0.57 | no | 0.67 | 1,819 / 1,792 |
| 93 | SL_Z6_R30_T100000_C5 | gas5_ics3 | 15 | - | 5,737 [4,898, 6,589] | 16.9 | 0.73 | no | 0.67 | 1,947 / 2,222 |
| 94 | SS_Z10_R20_T100000_C5 | gas5_ics3 | 15 | - | 5,793 [5,224, 6,310] | 14.1 | 0.67 | no | 0.67 | 2,045 / 1,818 |
| 95 | SL_Z6_R20_T100000_C5 | gas5_ics3 | 15 | - | 5,737 [4,902, 6,590] | 16.9 | 0.73 | no | 0.67 | 1,902 / 2,263 |
| 96 | SS_Z6_R20_T100000_C5 | gas5_ics3 | 15 | - | 5,962 [5,414, 6,474] | 23.5 | 0.84 | no | 0.67 | 2,678 / 2,310 |
| 97 | SS_Z6_R30_T100000_C2 | gas5_ics3 | 15 | - | 6,074 [5,453, 6,754] | 22.5 | 0.87 | no | 0.67 | 2,782 / 2,528 |
| 98 | SO_Z8_R30_T100000_C5 | gas5_ics3 | 15 | - | 6,243 [5,753, 6,784] | 16.5 | 0.64 | no | 0.67 | 2,075 / 1,910 |
| 99 | SO_Z8_R20_T100000_C2 | gas5_ics3 | 15 | - | 6,412 [5,803, 6,907] | 15.8 | 0.65 | no | 0.67 | 2,194 / 1,993 |
| 100 | SY_Z10_R20_T100000_C2 | gas5_ics3 | 15 | - | 5,512 [3,883, 6,561] | 9.9 | 0.62 | no | 0.68 | 1,627 / 1,781 |
| 101 | SY_Z6_R30_T100000_C2 | gas5_ics3 | 15 | - | 5,906 [5,058, 6,718] | 18.5 | 0.77 | no | 0.67 | 2,171 / 2,383 |
| 102 | SS_Z10_R20_T100000_C2 | gas5_ics3 | 15 | - | 5,849 [5,257, 6,631] | 13.2 | 0.69 | no | 0.67 | 2,149 / 1,873 |
| 103 | SS_Z6_R20_T100000_C2 | gas5_ics3 | 15 | - | 5,962 [5,451, 6,514] | 22.6 | 0.85 | no | 0.67 | 2,663 / 2,433 |
| 104 | SL_Z10_R20_T100000_C2 | gas5_ics3 | 15 | - | 5,343 [3,724, 6,682] | 8.2 | 0.54 | no | 0.68 | 1,254 / 1,636 |
| 105 | SY_Z6_R20_T100000_C2 | gas5_ics3 | 15 | - | 5,906 [5,064, 6,720] | 18.6 | 0.84 | no | 0.67 | 2,618 / 2,355 |
| 106 | SO_Z10_R20_T100000_C5 | gas5_ics3 | 15 | - | 6,243 [5,710, 6,743] | 14.0 | 0.61 | no | 0.67 | 2,029 / 1,751 |
| 107 | SO_Z10_R20_T100000_C2 | gas5_ics3 | 15 | - | 6,299 [5,754, 6,940] | 13.2 | 0.59 | no | 0.67 | 1,923 / 1,820 |
| 108 | SL_Z6_R30_T100000_C2 | gas5_ics3 | 15 | - | 5,849 [5,079, 6,800] | 16.2 | 0.74 | no | 0.67 | 2,066 / 2,291 |
| 109 | SO_Z6_R30_T100000_C5 | gas5_ics3 | 15 | - | 6,299 [5,833, 6,863] | 23.0 | 0.78 | no | 0.67 | 2,601 / 2,324 |
| 110 | SL_Z6_R20_T100000_C2 | gas5_ics3 | 15 | - | 5,906 [5,086, 6,796] | 16.3 | 0.74 | no | 0.67 | 2,037 / 2,305 |
| 111 | SO_Z6_R20_T100000_C5 | gas5_ics3 | 15 | - | 6,356 [5,836, 6,862] | 23.1 | 0.73 | no | 0.67 | 2,333 / 2,337 |
| 112 | SO_Z6_R30_T100000_C2 | gas5_ics3 | 15 | - | 6,412 [5,885, 6,911] | 22.3 | 0.76 | no | 0.67 | 2,512 / 2,338 |
| 113 | SO_Z6_R20_T100000_C2 | gas5_ics3 | 15 | - | 6,412 [5,886, 6,912] | 22.3 | 0.73 | no | 0.67 | 2,318 / 2,392 |
| 114 | SO_Z8_R20_T100000_C5 | gas5_ics3 | 15 | - | 6,187 [5,747, 6,789] | 16.6 | 0.66 | no | 0.67 | 2,135 / 1,965 |
| 115 | SS_Z4_R30_T100000_C5 | gas5_ics3 | 15 | - | 6,018 [5,536, 6,581] | 41.9 | 0.99 | no | 0.66 | 2,753 / 3,183 |
| 116 | SS_Z4_R30_T100000_C2 | gas5_ics3 | 15 | - | 6,074 [5,585, 6,713] | 40.9 | 0.93 | no | 0.66 | 2,529 / 3,144 |
| 117 | SY_Z4_R30_T100000_C2 | gas5_ics3 | 15 | - | 6,074 [5,323, 6,842] | 35.9 | 0.87 | no | 0.66 | 2,589 / 2,710 |
| 118 | SS_Z4_R20_T100000_C2 | gas5_ics3 | 15 | - | 6,131 [5,589, 6,809] | 41.0 | 0.96 | no | 0.66 | 2,723 / 3,182 |
| 119 | SY_Z4_R20_T100000_C2 | gas5_ics3 | 15 | - | 6,074 [5,324, 6,842] | 36.0 | 0.88 | no | 0.66 | 2,648 / 2,697 |
| 120 | SL_Z4_R30_T100000_C2 | gas5_ics3 | 15 | - | 6,131 [5,351, 6,833] | 33.2 | 0.83 | no | 0.66 | 2,395 / 2,672 |
| 121 | SL_Z4_R20_T100000_C2 | gas5_ics3 | 15 | - | 6,131 [5,352, 6,832] | 33.2 | 0.84 | no | 0.66 | 2,469 / 2,672 |
| 122 | SS_Z4_R20_T100000_C5 | gas5_ics3 | 15 | - | 5,906 [5,506, 6,586] | 42.1 | 0.93 | no | 0.66 | 2,679 / 2,839 |
| 123 | SO_Z4_R30_T100000_C2 | gas5_ics3 | 15 | - | 6,299 [5,938, 7,023] | 41.6 | 0.88 | no | 0.66 | 2,708 / 2,818 |
| 124 | SO_Z4_R20_T100000_C2 | gas5_ics3 | 15 | - | 6,524 [5,967, 6,994] | 41.7 | 0.86 | no | 0.66 | 2,723 / 2,868 |
| 125 | SY_Z4_R30_T100000_C5 | gas5_ics3 | 15 | - | 5,906 [5,156, 6,662] | 37.1 | 0.96 | no | 0.66 | 2,679 / 2,989 |
| 126 | SY_Z4_R20_T100000_C5 | gas5_ics3 | 15 | - | 5,906 [5,157, 6,662] | 37.1 | 0.91 | no | 0.66 | 2,679 / 2,722 |
| 127 | SL_Z4_R30_T100000_C5 | gas5_ics3 | 15 | - | 5,962 [5,165, 6,687] | 34.1 | 0.88 | no | 0.66 | 2,559 / 2,684 |
| 128 | SL_Z4_R20_T100000_C5 | gas5_ics3 | 15 | - | 5,962 [5,166, 6,687] | 34.2 | 0.88 | no | 0.66 | 2,485 / 2,748 |
| 129 | SO_Z4_R30_T100000_C5 | gas5_ics3 | 15 | - | 6,693 [5,815, 6,993] | 42.6 | 0.89 | no | 0.66 | 2,768 / 3,185 |
| 130 | SO_Z4_R20_T100000_C5 | gas5_ics3 | 15 | - | 6,412 [5,891, 6,916] | 42.9 | 0.88 | no | 0.66 | 2,768 / 2,881 |
