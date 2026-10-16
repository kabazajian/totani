# Stage C — 2026-09-19 (no_residual)

Consensus selection (median Δ(−2lnL) across 73 model rows; BIC threshold p·ln 1080 = 63 per template): **bubble_structured → erosita_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV**

| step | chosen | median Δ | min | max | threshold | runner-up (median Δ) |
|---|---|---|---|---|---|---|
| 1 | bubble_structured | 19,926 | 11,785 | 27,403 | 63 | bubble_flat (19,610) |
| 2 | erosita_outer_NS | 8,421 | 6,382 | 9,243 | 126 | loopI_shells+nps (7,784) |
| 3 | loopI_haslam_raw | 1,417 | 1,248 | 1,969 | 63 | loopI_shells+nps (717) |
| 4 | sun_moon_ic | 692 | 359 | 969 | 63 | loopIV (371) |
| 5 | loopIV | 97 | 7 | 368 | 63 | loopII (0) |
| 6 | (stop) loopII | 0 | -99 | 84 | 63 | loopIII (-0) |

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
| SS_Z10_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SS_Z10_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z10_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells |
| SS_Z10_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z10_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_outer_NS |
| SS_Z10_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z10_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_outer_NS |
| SS_Z10_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z10_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SS_Z10_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z10_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells |
| SS_Z10_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z10_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_outer_NS |
| SS_Z10_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
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
| SS_Z4_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells |
| SS_Z4_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z4_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells |
| SS_Z4_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z4_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SS_Z4_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z6_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SS_Z6_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z6_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SS_Z6_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z6_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells |
| SS_Z6_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z6_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SS_Z6_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z6_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SS_Z6_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z6_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SS_Z6_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z6_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SS_Z6_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z6_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells |
| SS_Z6_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z8_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SS_Z8_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z8_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SS_Z8_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z8_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells |
| SS_Z8_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z8_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SS_Z8_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z8_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SS_Z8_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z8_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells |
| SS_Z8_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z8_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells |
| SS_Z8_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z8_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells |
| SS_Z8_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| a2020 | - | bubble_totani → loopI_shells+nps |

## Ranking with the consensus set (bubble_structured, erosita_outer_NS, loopI_haslam_raw, sun_moon_ic, loopIV)

| rank | model | variant | n | Δ(−2lnL) vs Totani baseline | halo @ 20.8 | TS @ 20.8 | quad Σ/joint | any neg | inj min | N / S |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | gll_iem_v07 | - | 9 | -6,080 | 900 [245, 1,561] | 1.9 | -0.20 | yes | 0.91 | -84 / -99 |
| 2 | a2020 | - | 22 | -3,304 | 2,925 [2,040, 4,221] | 9.2 | -0.17 | yes | 1.00 | -577 / 72 |
| 3 | SS_Z10_R30_T150_C5 | gas5_ics3 | 16 | -1,880 | 1,912 [1,250, 2,585] | 3.3 | -0.29 | yes | 0.66 | 104 / -658 |
| 4 | SS_Z8_R30_T150_C5 | gas5_ics3 | 16 | -1,842 | 1,969 [1,316, 2,638] | 4.5 | -0.04 | yes | 0.66 | 402 / -480 |
| 5 | SS_Z10_R20_T150_C5 | gas5_ics3 | 16 | -1,834 | 1,912 [1,232, 2,567] | 3.3 | -0.34 | yes | 0.66 | -17 / -632 |
| 6 | SS_Z8_R20_T150_C5 | gas5_ics3 | 16 | -1,817 | 1,969 [1,309, 2,630] | 4.5 | -0.20 | yes | 0.67 | 73 / -467 |
| 7 | SS_Z6_R30_T150_C5 | gas5_ics3 | 16 | -1,809 | 2,137 [1,440, 2,747] | 6.2 | 0.26 | yes | 0.66 | 700 / -154 |
| 8 | SS_Z6_R20_T150_C5 | gas5_ics3 | 16 | -1,794 | 2,137 [1,437, 2,743] | 6.2 | -0.08 | yes | 0.66 | 73 / -236 |
| 9 | SL_Z6_R30_T150_C5 | gas5_ics3 | 16 | -1,755 | 1,856 [886, 2,764] | 2.1 | -0.58 | yes | 0.66 | -332 / -741 |
| 10 | SS_Z4_R30_T150_C5 | gas5_ics3 | 16 | -1,696 | 2,193 [1,537, 2,820] | 7.7 | -0.28 | yes | 0.65 | -477 / -128 |
| 11 | SS_Z4_R20_T150_C5 | gas5_ics3 | 16 | -1,690 | 2,250 [1,532, 2,823] | 7.7 | -0.26 | yes | 0.65 | -447 / -128 |
| 12 | SS_Z8_R30_T150_C2 | gas5_ics3 | 16 | -1,609 | 2,081 [1,439, 2,768] | 4.0 | -0.18 | yes | 0.66 | 251 / -624 |
| 13 | SS_Z10_R30_T150_C2 | gas5_ics3 | 16 | -1,609 | 1,744 [747, 2,655] | 2.9 | -0.46 | yes | 0.67 | 27 / -831 |
| 14 | SS_Z6_R30_T150_C2 | gas5_ics3 | 16 | -1,590 | 2,362 [1,525, 2,874] | 5.8 | -0.00 | yes | 0.66 | 311 / -313 |
| 15 | SL_Z10_R30_T150_C2 | gas5_ics3 | 16 | -1,581 | 900 [-905, 2,741] | 0.2 | -1.99 | yes | 0.67 | -572 / -1,217 |
| 16 | SS_Z10_R20_T150_C2 | gas5_ics3 | 16 | -1,578 | 1,744 [735, 2,637] | 2.9 | -0.34 | yes | 0.67 | 147 / -736 |
| 17 | SS_Z8_R20_T150_C2 | gas5_ics3 | 16 | -1,571 | 2,137 [1,422, 2,752] | 4.0 | -0.27 | yes | 0.66 | -17 / -571 |
| 18 | SL_Z6_R30_T150_C2 | gas5_ics3 | 16 | -1,545 | 2,193 [1,133, 3,039] | 2.0 | -0.38 | yes | 0.66 | -20 / -817 |
| 19 | SS_Z6_R20_T150_C2 | gas5_ics3 | 16 | -1,500 | 2,193 [1,543, 2,850] | 5.9 | 0.07 | yes | 0.66 | 460 / -300 |
| 20 | SL_Z6_R20_T150_C2 | gas5_ics3 | 16 | -1,457 | 1,800 [1,094, 3,068] | 2.2 | -0.50 | yes | 0.66 | -79 / -818 |
| 21 | SS_Z4_R30_T150_C2 | gas5_ics3 | 16 | -1,433 | 2,250 [1,635, 2,920] | 7.4 | -0.28 | yes | 0.65 | -477 / -152 |
| 22 | SS_Z4_R20_T150_C2 | gas5_ics3 | 16 | -1,429 | 2,193 [1,638, 2,936] | 7.5 | -0.32 | yes | 0.65 | -552 / -152 |
| 23 | SS_Z6_R30_T150_C5 | gas1_ics1 | 10 | -1,361 | 1,800 [1,123, 2,452] | 7.3 | -0.17 | yes | 0.94 | -5 / -302 |
| 24 | SS_Z6_R20_T150_C5 | gas1_ics1 | 10 | -1,335 | 1,800 [1,117, 2,446] | 7.2 | -0.18 | yes | 0.94 | -20 / -302 |
| 25 | SS_Z4_R30_T150_C5 | gas1_ics1 | 10 | -1,298 | 1,912 [1,236, 2,559] | 8.3 | -0.01 | yes | 0.94 | 99 / -124 |
| 26 | SS_Z4_R20_T150_C5 | gas1_ics1 | 10 | -1,289 | 1,912 [1,234, 2,558] | 8.3 | -0.01 | yes | 0.94 | 99 / -124 |
| 27 | SS_Z10_R30_T150_C5 | gas1_ics1 | 10 | -1,263 | 1,631 [967, 2,304] | 6.0 | -0.40 | yes | 0.95 | -154 / -505 |
| 28 | SS_Z8_R30_T150_C5 | gas1_ics1 | 10 | -1,245 | 1,687 [1,011, 2,345] | 6.4 | -0.31 | yes | 0.95 | -110 / -416 |
| 29 | SS_Z8_R20_T150_C5 | gas1_ics1 | 10 | -1,189 | 1,687 [996, 2,331] | 6.3 | -0.33 | yes | 0.95 | -140 / -416 |
| 30 | SS_Z10_R20_T150_C5 | gas1_ics1 | 10 | -1,179 | 1,631 [939, 2,277] | 5.8 | -0.42 | yes | 0.95 | -184 / -505 |
| 31 | SS_Z6_R30_T150_C2 | gas1_ics1 | 10 | -1,087 | 1,912 [1,228, 2,556] | 8.2 | -0.10 | yes | 0.94 | 143 / -340 |
| 32 | SS_Z6_R20_T150_C2 | gas1_ics1 | 10 | -1,065 | 1,912 [1,221, 2,550] | 8.1 | -0.10 | yes | 0.94 | 143 / -328 |
| 33 | SS_Z4_R30_T150_C2 | gas1_ics1 | 10 | -1,052 | 2,025 [1,345, 2,667] | 9.3 | 0.06 | yes | 0.94 | 262 / -149 |
| 34 | SS_Z4_R20_T150_C2 | gas1_ics1 | 10 | -1,044 | 2,025 [1,343, 2,665] | 9.3 | 0.06 | yes | 0.94 | 262 / -149 |
| 35 | SS_Z8_R30_T150_C2 | gas1_ics1 | 10 | -1,036 | 1,800 [1,138, 2,471] | 7.4 | -0.21 | yes | 0.94 | 69 / -454 |
| 36 | SS_Z10_R30_T150_C2 | gas1_ics1 | 10 | -995 | 1,744 [1,078, 2,414] | 6.9 | -0.29 | yes | 0.95 | 24 / -531 |
| 37 | SS_Z8_R30_T100000_C5 | gas5_ics3 | 16 | -972 | 2,362 [1,686, 3,030] | 7.7 | 0.02 | yes | 0.66 | 582 / -539 |
| 38 | SS_Z8_R20_T150_C2 | gas1_ics1 | 10 | -966 | 1,800 [1,112, 2,446] | 7.2 | -0.23 | yes | 0.94 | 39 / -454 |
| 39 | SS_Z8_R30_T100000_C2 | gas5_ics3 | 16 | -933 | 2,418 [1,725, 3,076] | 7.1 | -0.08 | yes | 0.52 | 447 / -632 |
| 40 | SS_Z8_R20_T100000_C2 | gas5_ics3 | 16 | -925 | 2,418 [1,717, 3,067] | 7.1 | -0.02 | yes | 0.66 | 521 / -579 |
| 41 | SS_Z10_R20_T150_C2 | gas1_ics1 | 10 | -920 | 1,744 [1,048, 2,385] | 6.6 | -0.31 | yes | 0.95 | -6 / -531 |
| 42 | SS_Z10_R30_T100000_C2 | gas5_ics3 | 16 | -892 | 2,306 [1,643, 3,010] | 5.8 | -0.13 | yes | 0.66 | 417 / -714 |
| 43 | SS_Z8_R20_T100000_C5 | gas5_ics3 | 16 | -859 | 2,362 [1,667, 3,009] | 7.6 | 0.01 | yes | 0.66 | 492 / -473 |
| 44 | SS_Z6_R30_T100000_C5 | gas5_ics3 | 16 | -851 | 2,418 [1,756, 3,078] | 9.5 | 0.26 | yes | 0.66 | 850 / -231 |
| 45 | SS_Z6_R20_T100000_C5 | gas5_ics3 | 16 | -822 | 2,418 [1,752, 3,073] | 9.5 | 0.30 | yes | 0.66 | 865 / -136 |
| 46 | SS_Z6_R30_T100000_C2 | gas5_ics3 | 16 | -809 | 2,475 [1,798, 3,125] | 9.0 | 0.28 | yes | 0.66 | 983 / -284 |
| 47 | SS_Z4_R20_T100000_C2 | gas5_ics3 | 16 | -789 | 2,531 [1,904, 3,199] | 10.8 | 0.02 | yes | 0.65 | 120 / -60 |
| 48 | SS_Z10_R20_T100000_C2 | gas5_ics3 | 16 | -789 | 2,306 [1,624, 2,988] | 5.8 | -0.38 | yes | 0.65 | -183 / -702 |
| 49 | SS_Z6_R20_T100000_C2 | gas5_ics3 | 16 | -783 | 2,475 [1,793, 3,120] | 9.0 | 0.32 | yes | 0.66 | 1,028 / -230 |
| 50 | SS_Z4_R30_T100000_C5 | gas5_ics3 | 16 | -765 | 2,475 [1,856, 3,150] | 11.3 | -0.03 | yes | 0.65 | -28 / -48 |
| 51 | SS_Z10_R20_T100000_C5 | gas5_ics3 | 16 | -763 | 2,250 [1,589, 2,944] | 6.5 | -0.27 | yes | 0.66 | -3 / -609 |
| 52 | SS_Z4_R30_T100000_C2 | gas5_ics3 | 16 | -758 | 2,475 [1,901, 3,206] | 10.8 | 0.03 | yes | 0.65 | 135 / -60 |
| 53 | SL_Z6_R30_T100000_C2 | gas5_ics3 | 16 | -740 | 1,800 [1,032, 2,813] | 3.4 | -0.34 | yes | 0.66 | 355 / -959 |
| 54 | SS_Z4_R20_T100000_C5 | gas5_ics3 | 16 | -567 | 2,475 [1,825, 3,118] | 11.2 | 0.09 | yes | 0.65 | 76 / 143 |
| 55 | SS_Z8_R30_T100000_C2 | gas1_ics1 | 10 | -368 | 1,856 [1,181, 2,520] | 7.7 | -0.08 | yes | 0.94 | 220 / -371 |
| 56 | SS_Z8_R30_T100000_C5 | gas1_ics1 | 10 | -358 | 1,800 [1,140, 2,480] | 7.4 | -0.09 | yes | 0.94 | 176 / -333 |
| 57 | SS_Z8_R20_T100000_C2 | gas1_ics1 | 10 | -311 | 1,856 [1,163, 2,504] | 7.5 | -0.08 | yes | 0.94 | 205 / -358 |
| 58 | SS_Z4_R30_T100000_C2 | gas1_ics1 | 10 | -289 | 2,081 [1,417, 2,743] | 9.9 | 0.20 | yes | 0.94 | 414 / -2 |
| 59 | SS_Z4_R20_T100000_C2 | gas1_ics1 | 10 | -277 | 2,081 [1,415, 2,740] | 9.9 | 0.20 | yes | 0.94 | 414 / -2 |
| 60 | SS_Z4_R30_T100000_C5 | gas1_ics1 | 10 | -271 | 2,025 [1,368, 2,694] | 9.5 | 0.16 | yes | 0.94 | 295 / 23 |
| 61 | SL_Z6_R30_T150_C5 | gas1_ics1 | 10 | -262 | 787 [119, 1,499] | 1.3 | -2.19 | yes | 0.96 | -901 / -820 |
| 62 | SS_Z6_R30_T100000_C2 | gas1_ics1 | 10 | -253 | 1,912 [1,246, 2,581] | 8.3 | 0.02 | yes | 0.94 | 265 / -219 |
| 63 | SS_Z10_R30_T100000_C2 | gas1_ics1 | 10 | -249 | 1,800 [1,100, 2,445] | 7.0 | -0.17 | yes | 0.95 | 160 / -460 |
| 64 | SS_Z6_R30_T100000_C5 | gas1_ics1 | 10 | -229 | 1,856 [1,205, 2,539] | 7.9 | 0.02 | yes | 0.94 | 236 / -193 |
| 65 | SS_Z8_R20_T100000_C5 | gas1_ics1 | 10 | -224 | 1,800 [1,109, 2,449] | 7.1 | -0.10 | yes | 0.95 | 146 / -333 |
| 66 | SS_Z6_R20_T100000_C2 | gas1_ics1 | 10 | -223 | 1,912 [1,238, 2,573] | 8.2 | 0.02 | yes | 0.94 | 265 / -219 |
| 67 | SS_Z6_R20_T100000_C5 | gas1_ics1 | 10 | -196 | 1,856 [1,196, 2,531] | 7.9 | 0.01 | yes | 0.94 | 206 / -193 |
| 68 | SS_Z10_R20_T100000_C2 | gas1_ics1 | 10 | -118 | 1,744 [1,067, 2,411] | 6.7 | -0.19 | yes | 0.95 | 115 / -447 |
| 69 | SS_Z10_R20_T100000_C5 | gas1_ics1 | 10 | -94 | 1,687 [1,032, 2,375] | 6.5 | -0.19 | yes | 0.95 | 87 / -409 |
| 70 | SS_Z4_R20_T100000_C5 | gas1_ics1 | 10 | -42 | 1,969 [1,328, 2,656] | 9.1 | 0.17 | yes | 0.94 | 310 / 23 |
| 71 | SL_Z10_R30_T150_C2 | gas1_ics1 | 10 | -5 | 731 [22, 1,436] | 1.1 | -2.96 | yes | 0.97 | -798 / -1,365 |
| 72 | SL_Z6_R30_T150_C2 | gas1_ics1 | 10 | +0 | 900 [219, 1,602] | 1.7 | -1.71 | yes | 0.96 | -678 / -859 |
| 73 | SL_Z6_R20_T150_C2 | gas1_ics1 | 10 | +20 | 900 [219, 1,601] | 1.7 | -1.71 | yes | 0.96 | -678 / -859 |
| 74 | SL_Z6_R30_T100000_C2 | gas1_ics1 | 10 | +460 | 956 [251, 1,669] | 1.8 | -1.49 | yes | 0.97 | -527 / -894 |
