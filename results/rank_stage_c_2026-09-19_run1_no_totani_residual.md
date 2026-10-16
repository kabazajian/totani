# Stage C — 2026-09-19

Consensus selection (median Δ(−2lnL) across 27 model rows; BIC threshold p·ln 1080 = 63 per template): **bubble_structured → erosita_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV**

| step | chosen | median Δ | min | max | threshold | runner-up (median Δ) |
|---|---|---|---|---|---|---|
| 1 | bubble_structured | 19,891 | 11,785 | 27,403 | 63 | bubble_flat (19,610) |
| 2 | erosita_outer_NS | 8,454 | 6,382 | 9,243 | 126 | loopI_haslam_raw (7,821) |
| 3 | loopI_haslam_raw | 1,599 | 1,306 | 1,969 | 63 | loopI_shells+nps (620) |
| 4 | sun_moon_ic | 814 | 359 | 969 | 63 | loopIV (438) |
| 5 | loopIV | 172 | 9 | 368 | 63 | loopII (0) |
| 6 | (stop) loopII | 0 | -86 | 84 | 63 | loopIII (-0) |

## Each model's own path

| model | variant | path |
|---|---|---|
| SL_Z10_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_outer_NS |
| SL_Z10_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z6_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_outer_NS |
| SL_Z6_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z6_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells |
| SL_Z6_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z6_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_outer_NS |
| SL_Z6_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z6_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_outer_NS → loopII |
| SL_Z6_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z4_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells |
| SS_Z4_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z4_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells |
| SS_Z4_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z4_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SS_Z4_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z4_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SS_Z4_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z4_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells |
| SS_Z4_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z4_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells |
| SS_Z4_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z4_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SS_Z4_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z6_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SS_Z6_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| a2020 | - | not a Paper II model code: 'a2020' |

## Ranking with the consensus set (bubble_structured, erosita_outer_NS, loopI_haslam_raw, sun_moon_ic, loopIV)

| rank | model | variant | n | Δ(−2lnL) vs Totani baseline | halo @ 20.8 | TS @ 20.8 | quad Σ/joint | any neg | inj min | N / S |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | gll_iem_v07 | - | 9 | -6,080 | 900 [245, 1,561] | 1.9 | -0.20 | yes | 0.91 | -84 / -99 |
| 2 | a2020 | - | 22 | -3,304 | 2,925 [2,040, 4,221] | 9.2 | -0.17 | yes | 1.00 | -577 / 72 |
| 3 | SL_Z6_R30_T150_C5 | gas5_ics3 | 16 | -1,755 | 1,856 [886, 2,764] | 2.1 | -0.58 | yes | 0.66 | -332 / -741 |
| 4 | SS_Z4_R30_T150_C5 | gas5_ics3 | 16 | -1,696 | 2,193 [1,537, 2,820] | 7.7 | -0.28 | yes | 0.65 | -477 / -128 |
| 5 | SS_Z4_R20_T150_C5 | gas5_ics3 | 16 | -1,690 | 2,250 [1,532, 2,823] | 7.7 | -0.26 | yes | 0.65 | -447 / -128 |
| 6 | SS_Z6_R30_T150_C2 | gas5_ics3 | 16 | -1,590 | 2,362 [1,525, 2,874] | 5.8 | -0.00 | yes | 0.66 | 311 / -313 |
| 7 | SL_Z10_R30_T150_C2 | gas5_ics3 | 16 | -1,581 | 900 [-905, 2,741] | 0.2 | -1.99 | yes | 0.67 | -572 / -1,217 |
| 8 | SL_Z6_R30_T150_C2 | gas5_ics3 | 16 | -1,545 | 2,193 [1,133, 3,039] | 2.0 | -0.38 | yes | 0.66 | -20 / -817 |
| 9 | SL_Z6_R20_T150_C2 | gas5_ics3 | 16 | -1,457 | 1,800 [1,094, 3,068] | 2.2 | -0.50 | yes | 0.66 | -79 / -818 |
| 10 | SS_Z4_R30_T150_C2 | gas5_ics3 | 16 | -1,433 | 2,250 [1,635, 2,920] | 7.4 | -0.28 | yes | 0.65 | -477 / -152 |
| 11 | SS_Z4_R20_T150_C2 | gas5_ics3 | 16 | -1,429 | 2,193 [1,638, 2,936] | 7.5 | -0.32 | yes | 0.65 | -552 / -152 |
| 12 | SS_Z4_R30_T150_C5 | gas1_ics1 | 10 | -1,298 | 1,912 [1,236, 2,559] | 8.3 | -0.01 | yes | 0.94 | 99 / -124 |
| 13 | SS_Z4_R20_T150_C5 | gas1_ics1 | 10 | -1,289 | 1,912 [1,234, 2,558] | 8.3 | -0.01 | yes | 0.94 | 99 / -124 |
| 14 | SS_Z6_R30_T150_C2 | gas1_ics1 | 10 | -1,087 | 1,912 [1,228, 2,556] | 8.2 | -0.10 | yes | 0.94 | 143 / -340 |
| 15 | SS_Z4_R30_T150_C2 | gas1_ics1 | 10 | -1,052 | 2,025 [1,345, 2,667] | 9.3 | 0.06 | yes | 0.94 | 262 / -149 |
| 16 | SS_Z4_R20_T150_C2 | gas1_ics1 | 10 | -1,044 | 2,025 [1,343, 2,665] | 9.3 | 0.06 | yes | 0.94 | 262 / -149 |
| 17 | SS_Z4_R20_T100000_C2 | gas5_ics3 | 16 | -789 | 2,531 [1,904, 3,199] | 10.8 | 0.02 | yes | 0.65 | 120 / -60 |
| 18 | SS_Z4_R30_T100000_C2 | gas5_ics3 | 16 | -758 | 2,475 [1,901, 3,206] | 10.8 | 0.03 | yes | 0.65 | 135 / -60 |
| 19 | SL_Z6_R30_T100000_C2 | gas5_ics3 | 16 | -740 | 1,800 [1,032, 2,813] | 3.4 | -0.34 | yes | 0.66 | 355 / -959 |
| 20 | SS_Z4_R20_T100000_C5 | gas5_ics3 | 16 | -567 | 2,475 [1,825, 3,118] | 11.2 | 0.09 | yes | 0.65 | 76 / 143 |
| 21 | SS_Z4_R30_T100000_C2 | gas1_ics1 | 10 | -289 | 2,081 [1,417, 2,743] | 9.9 | 0.20 | yes | 0.94 | 414 / -2 |
| 22 | SS_Z4_R20_T100000_C2 | gas1_ics1 | 10 | -277 | 2,081 [1,415, 2,740] | 9.9 | 0.20 | yes | 0.94 | 414 / -2 |
| 23 | SL_Z6_R30_T150_C5 | gas1_ics1 | 10 | -262 | 787 [119, 1,499] | 1.3 | -2.19 | yes | 0.96 | -901 / -820 |
| 24 | SS_Z4_R20_T100000_C5 | gas1_ics1 | 10 | -42 | 1,969 [1,328, 2,656] | 9.1 | 0.17 | yes | 0.94 | 310 / 23 |
| 25 | SL_Z10_R30_T150_C2 | gas1_ics1 | 10 | -5 | 731 [22, 1,436] | 1.1 | -2.96 | yes | 0.97 | -798 / -1,365 |
| 26 | SL_Z6_R30_T150_C2 | gas1_ics1 | 10 | +0 | 900 [219, 1,602] | 1.7 | -1.71 | yes | 0.96 | -678 / -859 |
| 27 | SL_Z6_R20_T150_C2 | gas1_ics1 | 10 | +20 | 900 [219, 1,601] | 1.7 | -1.71 | yes | 0.96 | -678 / -859 |
| 28 | SL_Z6_R30_T100000_C2 | gas1_ics1 | 10 | +460 | 956 [251, 1,669] | 1.8 | -1.49 | yes | 0.97 | -527 / -894 |
