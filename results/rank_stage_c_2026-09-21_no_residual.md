# Stage C — 2026-09-21 (no_residual)

Consensus selection (median Δ(−2lnL) across 257 model rows; BIC threshold p·ln 1080 = 63 per template): **bubble_flat → erosita_unif_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV**

| step | chosen | median Δ | min | max | threshold | runner-up (median Δ) |
|---|---|---|---|---|---|---|
| 1 | bubble_flat | 19,194 | 12,490 | 26,944 | 63 | bubble_structured (18,996) |
| 2 | erosita_unif_outer_NS | 10,381 | 8,256 | 12,214 | 126 | erosita_outer_NS (10,366) |
| 3 | loopI_haslam_raw | 990 | 761 | 1,922 | 63 | loopI_shells+nps (666) |
| 4 | sun_moon_ic | 348 | 140 | 744 | 63 | loopIV (210) |
| 5 | loopIV | 88 | 18 | 201 | 63 | loopII (1) |
| 6 | (stop) loopII | 1 | -24 | 20 | 63 | loopIII (-0) |

## Each model's own path

| model | variant | path |
|---|---|---|
| SL_Z10_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z10_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SL_Z10_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z10_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z10_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SL_Z10_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z10_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z10_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z10_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SL_Z10_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z10_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z10_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z10_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z10_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z10_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z10_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z4_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SL_Z4_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SL_Z4_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SL_Z4_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z4_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SL_Z4_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z4_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SL_Z4_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z4_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SL_Z4_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SL_Z4_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SL_Z4_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z4_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SL_Z4_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z4_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SL_Z4_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SL_Z6_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z6_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z6_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z6_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z6_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SL_Z6_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z6_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z6_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z6_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z6_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z6_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z6_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z6_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SL_Z6_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z6_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z6_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z8_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SL_Z8_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z8_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z8_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z8_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z8_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z8_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SL_Z8_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z8_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SL_Z8_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z8_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SL_Z8_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z8_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SL_Z8_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z8_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SL_Z8_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SO_Z10_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SO_Z10_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z10_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SO_Z10_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z10_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z10_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SO_Z10_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SO_Z10_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z10_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SO_Z10_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SO_Z10_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SO_Z10_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z10_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z10_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SO_Z10_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SO_Z10_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z4_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → sun_moon_ic |
| SO_Z4_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SO_Z4_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z4_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SO_Z4_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z4_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SO_Z4_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z4_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SO_Z4_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → sun_moon_ic |
| SO_Z4_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z4_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z4_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SO_Z4_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z4_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z4_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z4_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SO_Z6_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z6_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z6_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SO_Z6_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z6_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → sun_moon_ic |
| SO_Z6_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SO_Z6_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic → loopII |
| SO_Z6_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z6_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SO_Z6_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z6_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SO_Z6_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z6_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → sun_moon_ic → loopII |
| SO_Z6_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SO_Z6_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → sun_moon_ic |
| SO_Z6_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z8_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SO_Z8_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z8_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SO_Z8_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z8_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → sun_moon_ic |
| SO_Z8_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SO_Z8_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → sun_moon_ic → loopII |
| SO_Z8_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SO_Z8_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SO_Z8_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z8_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SO_Z8_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z8_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → sun_moon_ic |
| SO_Z8_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SO_Z8_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → sun_moon_ic → loopII |
| SO_Z8_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SS_Z10_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z10_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SS_Z10_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SS_Z10_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SS_Z10_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SS_Z10_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z10_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SS_Z10_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z10_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z10_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SS_Z10_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z10_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SS_Z10_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SS_Z10_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z10_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SS_Z10_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SS_Z4_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SS_Z4_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z4_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SS_Z4_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SS_Z4_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SS_Z4_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z4_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SS_Z4_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SS_Z4_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SS_Z4_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SS_Z4_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SS_Z4_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SS_Z4_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → sun_moon_ic |
| SS_Z4_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z4_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SS_Z4_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → erosita_unif_NS → sun_moon_ic |
| SS_Z6_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z6_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SS_Z6_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z6_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SS_Z6_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SS_Z6_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z6_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z6_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SS_Z6_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z6_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SS_Z6_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z6_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SS_Z6_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z6_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z6_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → sun_moon_ic |
| SS_Z6_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z8_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z8_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SS_Z8_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z8_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SS_Z8_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SS_Z8_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z8_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z8_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z8_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z8_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SS_Z8_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SS_Z8_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SS_Z8_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SS_Z8_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z8_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SS_Z8_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z10_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SY_Z10_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SY_Z10_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SY_Z10_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SY_Z10_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SY_Z10_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z10_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SY_Z10_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z10_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SY_Z10_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z10_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SY_Z10_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z10_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SY_Z10_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z10_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SY_Z10_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z4_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SY_Z4_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SY_Z4_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SY_Z4_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SY_Z4_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SY_Z4_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z4_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SY_Z4_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SY_Z4_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SY_Z4_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SY_Z4_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SY_Z4_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SY_Z4_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SY_Z4_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z4_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SY_Z4_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SY_Z6_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SY_Z6_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z6_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SY_Z6_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SY_Z6_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SY_Z6_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z6_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SY_Z6_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z6_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SY_Z6_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SY_Z6_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SY_Z6_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SY_Z6_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SY_Z6_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z6_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SY_Z6_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SY_Z8_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SY_Z8_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SY_Z8_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SY_Z8_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SY_Z8_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SY_Z8_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z8_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SY_Z8_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z8_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SY_Z8_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z8_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SY_Z8_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SY_Z8_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SY_Z8_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z8_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SY_Z8_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| a2020 | - | bubble_totani → loopI_shells+nps |

## Ranking with the consensus set (bubble_flat, erosita_unif_outer_NS, loopI_haslam_raw, sun_moon_ic, loopIV)

| rank | model | variant | n | Δ(−2lnL) vs Totani baseline | halo @ 20.8 | TS @ 20.8 | quad Σ/joint | any neg | inj min | N / S |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | gll_iem_v07 | - | 9 | -4,186 | 2,025 [1,429, 2,578] | 12.3 | 0.82 | yes | 0.90 | 637 / 1,020 |
| 2 | a2020 | - | 22 | -1,817 | 7,255 [6,675, 7,797] | 48.4 | 1.01 | no | 0.98 | 3,385 / 3,947 |
| 3 | SS_Z10_R30_T150_C5 | gas5_ics3 | 16 | -1,354 | 5,962 [5,227, 6,716] | 11.2 | 0.63 | no | 0.66 | 2,180 / 1,593 |
| 4 | SY_Z10_R30_T150_C5 | gas5_ics3 | 16 | -1,347 | 5,849 [4,984, 6,708] | 9.2 | 0.60 | no | 0.66 | 2,031 / 1,488 |
| 5 | SL_Z10_R30_T150_C5 | gas5_ics3 | 16 | -1,330 | 5,906 [5,036, 6,726] | 8.6 | 0.54 | no | 0.66 | 1,883 / 1,293 |
| 6 | SY_Z10_R20_T150_C5 | gas5_ics3 | 16 | -1,306 | 5,849 [4,991, 6,705] | 9.4 | 0.61 | no | 0.66 | 2,061 / 1,529 |
| 7 | SS_Z10_R20_T150_C5 | gas5_ics3 | 16 | -1,305 | 5,962 [5,223, 6,706] | 11.4 | 0.63 | no | 0.66 | 2,165 / 1,596 |
| 8 | SL_Z10_R20_T150_C5 | gas5_ics3 | 16 | -1,295 | 5,906 [5,042, 6,716] | 8.9 | 0.54 | no | 0.66 | 1,883 / 1,308 |
| 9 | SO_Z10_R30_T150_C5 | gas5_ics3 | 16 | -1,288 | 6,299 [5,694, 6,917] | 11.3 | 0.58 | no | 0.66 | 2,164 / 1,513 |
| 10 | SS_Z8_R30_T150_C5 | gas5_ics3 | 16 | -1,283 | 6,018 [5,278, 6,726] | 14.0 | 0.68 | no | 0.65 | 2,256 / 1,863 |
| 11 | SY_Z8_R30_T150_C5 | gas5_ics3 | 16 | -1,273 | 5,906 [5,049, 6,725] | 11.7 | 0.66 | no | 0.66 | 2,108 / 1,784 |
| 12 | SS_Z6_R30_T150_C5 | gas5_ics3 | 16 | -1,268 | 6,074 [5,396, 6,780] | 19.7 | 0.79 | no | 0.65 | 2,529 / 2,291 |
| 13 | SY_Z6_R30_T150_C5 | gas5_ics3 | 16 | -1,268 | 6,018 [5,203, 6,788] | 17.2 | 0.78 | no | 0.65 | 2,469 / 2,213 |
| 14 | SY_Z6_R20_T150_C5 | gas5_ics3 | 16 | -1,257 | 6,018 [5,206, 6,787] | 17.3 | 0.76 | no | 0.65 | 2,364 / 2,227 |
| 15 | SO_Z10_R20_T150_C5 | gas5_ics3 | 16 | -1,253 | 6,299 [5,690, 6,905] | 11.5 | 0.58 | no | 0.66 | 2,149 / 1,516 |
| 16 | SS_Z6_R20_T150_C5 | gas5_ics3 | 16 | -1,252 | 6,074 [5,397, 6,778] | 19.8 | 0.79 | no | 0.65 | 2,529 / 2,291 |
| 17 | SL_Z8_R30_T150_C5 | gas5_ics3 | 16 | -1,250 | 5,962 [5,081, 6,703] | 10.9 | 0.62 | no | 0.65 | 1,989 / 1,678 |
| 18 | SS_Z8_R20_T150_C5 | gas5_ics3 | 16 | -1,249 | 6,018 [5,281, 6,722] | 14.1 | 0.69 | no | 0.65 | 2,256 / 1,877 |
| 19 | SY_Z8_R20_T150_C5 | gas5_ics3 | 16 | -1,247 | 5,906 [5,054, 6,723] | 11.9 | 0.67 | no | 0.65 | 2,167 / 1,798 |
| 20 | SL_Z6_R30_T150_C5 | gas5_ics3 | 16 | -1,245 | 6,074 [5,227, 6,694] | 16.1 | 0.72 | no | 0.65 | 2,215 / 2,148 |
| 21 | SL_Z8_R20_T150_C5 | gas5_ics3 | 16 | -1,229 | 5,962 [5,084, 6,698] | 11.0 | 0.62 | no | 0.65 | 2,018 / 1,692 |
| 22 | SL_Z6_R20_T150_C5 | gas5_ics3 | 16 | -1,228 | 6,018 [5,227, 6,690] | 16.2 | 0.73 | no | 0.65 | 2,215 / 2,161 |
| 23 | SO_Z8_R30_T150_C5 | gas5_ics3 | 16 | -1,208 | 6,299 [5,703, 6,886] | 13.6 | 0.64 | no | 0.65 | 2,239 / 1,782 |
| 24 | SO_Z8_R20_T150_C5 | gas5_ics3 | 16 | -1,185 | 6,299 [5,702, 6,880] | 13.8 | 0.64 | no | 0.65 | 2,254 / 1,809 |
| 25 | SO_Z6_R30_T150_C5 | gas5_ics3 | 16 | -1,170 | 6,299 [5,750, 6,874] | 19.2 | 0.73 | no | 0.65 | 2,362 / 2,237 |
| 26 | SO_Z6_R20_T150_C5 | gas5_ics3 | 16 | -1,146 | 6,299 [5,746, 6,868] | 19.3 | 0.73 | no | 0.65 | 2,377 / 2,250 |
| 27 | SS_Z10_R30_T150_C2 | gas5_ics3 | 16 | -1,071 | 6,018 [4,515, 6,894] | 10.8 | 0.64 | no | 0.66 | 2,254 / 1,597 |
| 28 | SY_Z10_R30_T150_C2 | gas5_ics3 | 16 | -1,065 | 6,018 [4,250, 6,920] | 8.8 | 0.61 | no | 0.66 | 2,121 / 1,529 |
| 29 | SS_Z8_R30_T150_C2 | gas5_ics3 | 16 | -1,056 | 6,187 [5,452, 6,915] | 13.3 | 0.66 | no | 0.65 | 2,256 / 1,850 |
| 30 | SL_Z10_R30_T150_C2 | gas5_ics3 | 16 | -1,050 | 6,131 [4,321, 6,865] | 8.3 | 0.53 | no | 0.66 | 1,913 / 1,307 |
| 31 | SY_Z8_R30_T150_C2 | gas5_ics3 | 16 | -1,040 | 6,074 [5,250, 6,929] | 11.1 | 0.65 | no | 0.66 | 2,198 / 1,759 |
| 32 | SS_Z10_R20_T150_C2 | gas5_ics3 | 16 | -1,033 | 6,018 [4,548, 6,882] | 10.9 | 0.64 | no | 0.66 | 2,210 / 1,612 |
| 33 | SY_Z10_R20_T150_C2 | gas5_ics3 | 16 | -1,031 | 6,018 [4,300, 6,910] | 9.0 | 0.61 | no | 0.66 | 2,121 / 1,531 |
| 34 | SO_Z10_R30_T150_C2 | gas5_ics3 | 16 | -1,021 | 6,468 [5,831, 7,072] | 10.8 | 0.58 | no | 0.66 | 2,209 / 1,516 |
| 35 | SL_Z10_R20_T150_C2 | gas5_ics3 | 16 | -1,020 | 6,131 [5,253, 6,853] | 8.5 | 0.54 | no | 0.66 | 1,928 / 1,386 |
| 36 | SS_Z8_R20_T150_C2 | gas5_ics3 | 16 | -1,010 | 6,187 [5,444, 6,900] | 13.5 | 0.68 | no | 0.65 | 2,331 / 1,904 |
| 37 | SL_Z8_R30_T150_C2 | gas5_ics3 | 16 | -1,006 | 6,187 [5,300, 6,842] | 10.5 | 0.61 | no | 0.65 | 2,093 / 1,692 |
| 38 | SS_Z4_R30_T150_C5 | gas5_ics3 | 16 | -998 | 6,131 [5,450, 6,746] | 35.9 | 0.85 | no | 0.64 | 2,409 / 2,794 |
| 39 | SO_Z10_R20_T150_C2 | gas5_ics3 | 16 | -993 | 6,412 [5,279, 7,057] | 11.0 | 0.58 | no | 0.66 | 2,194 / 1,557 |
| 40 | SS_Z4_R20_T150_C5 | gas5_ics3 | 16 | -989 | 6,018 [5,446, 6,749] | 35.9 | 0.86 | no | 0.64 | 2,364 / 2,794 |
| 41 | SY_Z8_R20_T150_C2 | gas5_ics3 | 16 | -984 | 6,074 [5,242, 6,912] | 11.4 | 0.66 | no | 0.66 | 2,213 / 1,826 |
| 42 | SY_Z4_R30_T150_C5 | gas5_ics3 | 16 | -984 | 6,074 [5,340, 6,741] | 32.6 | 0.80 | no | 0.64 | 2,200 / 2,664 |
| 43 | SS_Z6_R30_T150_C2 | gas5_ics3 | 16 | -982 | 6,243 [5,541, 6,940] | 19.1 | 0.78 | no | 0.65 | 2,543 / 2,332 |
| 44 | SY_Z4_R20_T150_C5 | gas5_ics3 | 16 | -979 | 6,074 [5,340, 6,740] | 32.6 | 0.82 | no | 0.64 | 2,320 / 2,677 |
| 45 | SY_Z6_R30_T150_C2 | gas5_ics3 | 16 | -978 | 6,187 [5,374, 6,942] | 16.6 | 0.76 | no | 0.65 | 2,439 / 2,280 |
| 46 | SS_Z6_R20_T150_C2 | gas5_ics3 | 16 | -971 | 6,243 [5,542, 6,938] | 19.1 | 0.79 | no | 0.65 | 2,558 / 2,345 |
| 47 | SY_Z6_R20_T150_C2 | gas5_ics3 | 16 | -969 | 6,187 [5,377, 6,939] | 16.7 | 0.75 | no | 0.65 | 2,364 / 2,267 |
| 48 | SL_Z8_R20_T150_C2 | gas5_ics3 | 16 | -957 | 6,187 [5,289, 6,828] | 10.6 | 0.61 | no | 0.65 | 2,093 / 1,669 |
| 49 | SL_Z6_R30_T150_C2 | gas5_ics3 | 16 | -952 | 6,299 [5,770, 6,810] | 15.6 | 0.70 | no | 0.65 | 2,230 / 2,175 |
| 50 | SO_Z8_R30_T150_C2 | gas5_ics3 | 16 | -948 | 6,412 [5,834, 7,031] | 13.2 | 0.62 | no | 0.66 | 2,180 / 1,823 |
| 51 | SL_Z6_R20_T150_C2 | gas5_ics3 | 16 | -945 | 6,299 [5,767, 6,807] | 15.6 | 0.70 | no | 0.65 | 2,215 / 2,201 |
| 52 | SL_Z4_R30_T150_C5 | gas5_ics3 | 16 | -939 | 6,018 [5,515, 6,633] | 31.0 | 0.79 | no | 0.65 | 2,126 / 2,652 |
| 53 | SL_Z4_R20_T150_C5 | gas5_ics3 | 16 | -935 | 6,018 [5,516, 6,631] | 31.0 | 0.79 | no | 0.65 | 2,096 / 2,652 |
| 54 | SO_Z8_R20_T150_C2 | gas5_ics3 | 16 | -931 | 6,412 [5,832, 7,023] | 13.4 | 0.64 | no | 0.66 | 2,269 / 1,837 |
| 55 | SO_Z6_R30_T150_C2 | gas5_ics3 | 16 | -906 | 6,468 [5,873, 7,011] | 18.7 | 0.73 | no | 0.65 | 2,452 / 2,277 |
| 56 | SO_Z6_R20_T150_C2 | gas5_ics3 | 16 | -898 | 6,412 [5,872, 7,008] | 18.7 | 0.73 | no | 0.65 | 2,422 / 2,277 |
| 57 | SO_Z4_R30_T150_C5 | gas5_ics3 | 16 | -851 | 6,299 [5,808, 6,945] | 36.6 | 0.83 | no | 0.64 | 2,394 / 2,861 |
| 58 | SO_Z4_R20_T150_C5 | gas5_ics3 | 16 | -847 | 6,580 [5,776, 6,975] | 36.5 | 0.80 | no | 0.64 | 2,394 / 2,861 |
| 59 | SS_Z4_R30_T150_C2 | gas5_ics3 | 16 | -731 | 6,243 [5,585, 6,891] | 35.2 | 0.83 | no | 0.64 | 2,379 / 2,793 |
| 60 | SS_Z4_R20_T150_C2 | gas5_ics3 | 16 | -725 | 6,243 [5,586, 6,890] | 35.2 | 0.83 | no | 0.64 | 2,379 / 2,780 |
| 61 | SY_Z4_R30_T150_C2 | gas5_ics3 | 16 | -705 | 6,243 [5,498, 6,863] | 32.0 | 0.80 | no | 0.64 | 2,305 / 2,689 |
| 62 | SS_Z10_R30_T150_C5 | gas1_ics1 | 10 | -704 | 5,343 [4,828, 5,882] | 104.5 | 0.90 | no | 0.90 | 1,558 / 3,240 |
| 63 | SY_Z4_R20_T150_C2 | gas5_ics3 | 16 | -701 | 6,243 [5,498, 6,862] | 32.1 | 0.80 | no | 0.64 | 2,319 / 2,689 |
| 64 | SS_Z8_R30_T150_C5 | gas1_ics1 | 10 | -667 | 5,399 [4,889, 5,941] | 107.3 | 0.91 | no | 0.90 | 1,618 / 3,279 |
| 65 | SS_Z6_R30_T150_C5 | gas1_ics1 | 10 | -662 | 5,512 [5,007, 6,057] | 112.7 | 0.91 | no | 0.89 | 1,737 / 3,292 |
| 66 | SL_Z4_R30_T150_C2 | gas5_ics3 | 16 | -652 | 6,187 [5,631, 6,751] | 30.6 | 0.76 | no | 0.65 | 2,066 / 2,651 |
| 67 | SS_Z10_R20_T150_C5 | gas1_ics1 | 10 | -651 | 5,343 [4,809, 5,863] | 103.7 | 0.90 | no | 0.90 | 1,543 / 3,240 |
| 68 | SL_Z4_R20_T150_C2 | gas5_ics3 | 16 | -649 | 6,187 [5,630, 6,750] | 30.6 | 0.77 | no | 0.65 | 2,125 / 2,651 |
| 69 | SS_Z6_R20_T150_C5 | gas1_ics1 | 10 | -647 | 5,512 [5,003, 6,053] | 112.5 | 0.91 | no | 0.89 | 1,737 / 3,292 |
| 70 | SY_Z6_R30_T150_C5 | gas1_ics1 | 10 | -639 | 4,837 [4,301, 5,382] | 80.9 | 0.91 | no | 0.91 | 1,244 / 3,161 |
| 71 | SS_Z8_R20_T150_C5 | gas1_ics1 | 10 | -631 | 5,399 [4,881, 5,933] | 107.0 | 0.90 | no | 0.90 | 1,618 / 3,253 |
| 72 | SY_Z6_R20_T150_C5 | gas1_ics1 | 10 | -618 | 4,837 [4,301, 5,382] | 80.9 | 0.91 | no | 0.91 | 1,244 / 3,161 |
| 73 | SO_Z4_R30_T150_C2 | gas5_ics3 | 16 | -607 | 6,468 [5,923, 7,047] | 36.0 | 0.81 | no | 0.64 | 2,409 / 2,848 |
| 74 | SO_Z4_R20_T150_C2 | gas5_ics3 | 16 | -604 | 6,468 [5,922, 7,046] | 36.0 | 0.81 | no | 0.64 | 2,394 / 2,848 |
| 75 | SY_Z10_R30_T150_C5 | gas1_ics1 | 10 | -592 | 4,612 [4,076, 5,176] | 71.3 | 0.91 | no | 0.91 | 1,080 / 3,109 |
| 76 | SY_Z8_R30_T150_C5 | gas1_ics1 | 10 | -573 | 4,724 [4,153, 5,246] | 74.7 | 0.90 | no | 0.91 | 1,125 / 3,122 |
| 77 | SY_Z8_R30_T100000_C5 | gas5_ics3 | 16 | -539 | 5,568 [4,759, 6,455] | 14.5 | 0.70 | no | 0.65 | 2,106 / 1,804 |
| 78 | SL_Z8_R30_T100000_C5 | gas5_ics3 | 16 | -526 | 5,624 [4,747, 6,520] | 11.1 | 0.64 | no | 0.66 | 1,912 / 1,712 |
| 79 | SY_Z8_R20_T150_C5 | gas1_ics1 | 10 | -524 | 4,668 [4,148, 5,240] | 74.6 | 0.91 | no | 0.91 | 1,110 / 3,122 |
| 80 | SY_Z10_R20_T150_C5 | gas1_ics1 | 10 | -518 | 4,612 [4,070, 5,168] | 71.3 | 0.90 | no | 0.91 | 1,050 / 3,096 |
| 81 | SS_Z8_R30_T100000_C5 | gas5_ics3 | 16 | -502 | 5,793 [5,259, 6,337] | 16.3 | 0.73 | no | 0.65 | 2,345 / 1,896 |
| 82 | SY_Z8_R30_T100000_C2 | gas5_ics3 | 16 | -475 | 5,793 [4,920, 6,654] | 12.3 | 0.68 | no | 0.66 | 2,061 / 1,873 |
| 83 | SS_Z8_R30_T100000_C2 | gas5_ics3 | 16 | -470 | 5,906 [5,301, 6,685] | 15.4 | 0.74 | no | 0.65 | 2,375 / 1,978 |
| 84 | SL_Z8_R30_T100000_C2 | gas5_ics3 | 16 | -455 | 5,849 [4,949, 6,742] | 10.6 | 0.64 | no | 0.66 | 1,942 / 1,781 |
| 85 | SY_Z8_R20_T100000_C2 | gas5_ics3 | 16 | -437 | 5,793 [4,932, 6,658] | 12.5 | 0.69 | no | 0.66 | 2,121 / 1,900 |
| 86 | SS_Z4_R30_T150_C5 | gas1_ics1 | 10 | -433 | 5,681 [5,133, 6,182] | 118.3 | 0.92 | no | 0.89 | 1,886 / 3,331 |
| 87 | SS_Z4_R20_T150_C5 | gas1_ics1 | 10 | -427 | 5,681 [5,132, 6,182] | 118.3 | 0.92 | no | 0.89 | 1,886 / 3,331 |
| 88 | SY_Z4_R30_T150_C5 | gas1_ics1 | 10 | -424 | 5,006 [4,484, 5,553] | 88.9 | 0.92 | no | 0.90 | 1,409 / 3,200 |
| 89 | SL_Z8_R20_T100000_C2 | gas5_ics3 | 16 | -423 | 5,849 [4,960, 6,742] | 10.8 | 0.64 | no | 0.66 | 1,957 / 1,808 |
| 90 | SS_Z8_R20_T100000_C2 | gas5_ics3 | 16 | -422 | 5,962 [5,303, 6,684] | 15.6 | 0.74 | no | 0.65 | 2,404 / 1,993 |
| 91 | SY_Z10_R30_T100000_C5 | gas5_ics3 | 16 | -417 | 5,512 [4,661, 6,388] | 12.0 | 0.65 | no | 0.66 | 1,910 / 1,657 |
| 92 | SY_Z4_R20_T150_C5 | gas1_ics1 | 10 | -417 | 5,006 [4,484, 5,553] | 88.9 | 0.92 | no | 0.90 | 1,409 / 3,200 |
| 93 | SL_Z10_R30_T100000_C5 | gas5_ics3 | 16 | -403 | 5,568 [4,646, 6,466] | 8.9 | 0.57 | no | 0.66 | 1,686 / 1,488 |
| 94 | SO_Z8_R30_T100000_C2 | gas5_ics3 | 16 | -392 | 6,243 [5,714, 6,755] | 15.2 | 0.67 | no | 0.65 | 2,283 / 1,925 |
| 95 | SS_Z10_R30_T100000_C5 | gas5_ics3 | 16 | -391 | 5,737 [5,184, 6,273] | 13.7 | 0.71 | no | 0.65 | 2,254 / 1,803 |
| 96 | SS_Z10_R30_T100000_C2 | gas5_ics3 | 16 | -391 | 5,849 [5,224, 6,637] | 12.8 | 0.68 | no | 0.66 | 2,179 / 1,819 |
| 97 | SS_Z8_R20_T100000_C5 | gas5_ics3 | 16 | -390 | 5,793 [5,255, 6,332] | 16.5 | 0.76 | no | 0.65 | 2,420 / 1,964 |
| 98 | SY_Z10_R30_T100000_C2 | gas5_ics3 | 16 | -387 | 5,681 [3,970, 6,596] | 9.9 | 0.65 | no | 0.66 | 1,954 / 1,713 |
| 99 | SY_Z8_R20_T100000_C5 | gas5_ics3 | 16 | -387 | 5,568 [4,758, 6,446] | 13.3 | 0.73 | no | 0.65 | 2,165 / 1,886 |
| 100 | SL_Z10_R30_T100000_C2 | gas5_ics3 | 16 | -356 | 5,793 [4,857, 6,698] | 8.4 | 0.58 | no | 0.66 | 1,776 / 1,568 |
| 101 | SS_Z8_R30_T150_C2 | gas1_ics1 | 10 | -353 | 5,512 [4,973, 6,024] | 110.7 | 0.92 | no | 0.90 | 1,811 / 3,253 |
| 102 | SY_Z6_R30_T100000_C5 | gas5_ics3 | 16 | -350 | 5,737 [4,912, 6,550] | 19.2 | 0.85 | no | 0.65 | 2,603 / 2,300 |
| 103 | SL_Z8_R20_T100000_C5 | gas5_ics3 | 16 | -344 | 5,624 [4,740, 6,505] | 11.4 | 0.66 | no | 0.66 | 1,927 / 1,807 |
| 104 | SS_Z10_R30_T150_C2 | gas1_ics1 | 10 | -344 | 5,456 [4,904, 5,958] | 107.4 | 0.91 | no | 0.90 | 1,736 / 3,240 |
| 105 | SO_Z10_R30_T150_C5 | gas1_ics1 | 10 | -344 | 5,849 [5,314, 6,357] | 127.9 | 0.87 | no | 0.89 | 1,693 / 3,371 |
| 106 | SL_Z6_R30_T150_C5 | gas1_ics1 | 10 | -341 | 4,443 [3,915, 5,026] | 65.4 | 0.92 | no | 0.92 | 961 / 3,109 |
| 107 | SY_Z10_R20_T100000_C5 | gas5_ics3 | 16 | -329 | 5,512 [4,672, 6,392] | 12.3 | 0.67 | no | 0.66 | 1,955 / 1,712 |
| 108 | SL_Z10_R20_T100000_C5 | gas5_ics3 | 16 | -324 | 5,568 [4,661, 6,468] | 9.2 | 0.59 | no | 0.66 | 1,731 / 1,555 |
| 109 | SO_Z10_R30_T100000_C5 | gas5_ics3 | 16 | -321 | 6,131 [5,625, 6,667] | 13.6 | 0.65 | no | 0.66 | 2,163 / 1,806 |
| 110 | SO_Z10_R20_T150_C5 | gas1_ics1 | 10 | -317 | 5,849 [5,306, 6,349] | 127.5 | 0.86 | no | 0.89 | 1,663 / 3,358 |
| 111 | SO_Z10_R30_T100000_C2 | gas5_ics3 | 16 | -316 | 6,018 [5,220, 6,705] | 14.6 | 0.65 | no | 0.66 | 2,102 / 1,792 |
| 112 | SS_Z8_R20_T150_C2 | gas1_ics1 | 10 | -314 | 5,456 [4,957, 6,010] | 110.0 | 0.93 | no | 0.90 | 1,796 / 3,253 |
| 113 | SL_Z6_R20_T150_C5 | gas1_ics1 | 10 | -314 | 4,443 [3,914, 5,025] | 65.4 | 0.91 | no | 0.92 | 916 / 3,109 |
| 114 | SS_Z6_R30_T100000_C5 | gas5_ics3 | 16 | -313 | 5,849 [5,349, 6,415] | 22.8 | 0.87 | no | 0.65 | 2,767 / 2,324 |
| 115 | SY_Z6_R20_T100000_C5 | gas5_ics3 | 16 | -311 | 5,737 [4,917, 6,552] | 19.3 | 0.84 | no | 0.65 | 2,543 / 2,274 |
| 116 | SL_Z6_R30_T100000_C5 | gas5_ics3 | 16 | -306 | 5,737 [4,909, 6,602] | 17.0 | 0.78 | no | 0.65 | 2,364 / 2,101 |
| 117 | SS_Z6_R30_T150_C2 | gas1_ics1 | 10 | -306 | 5,624 [5,079, 6,128] | 115.7 | 0.93 | no | 0.89 | 1,916 / 3,292 |
| 118 | SS_Z10_R20_T150_C2 | gas1_ics1 | 10 | -303 | 5,399 [4,884, 5,937] | 106.6 | 0.92 | no | 0.90 | 1,722 / 3,226 |
| 119 | SY_Z6_R30_T150_C2 | gas1_ics1 | 10 | -300 | 4,893 [4,371, 5,453] | 83.0 | 0.93 | no | 0.91 | 1,423 / 3,147 |
| 120 | SS_Z10_R20_T100000_C5 | gas5_ics3 | 16 | -298 | 5,737 [5,183, 6,271] | 13.9 | 0.70 | no | 0.65 | 2,239 / 1,791 |
| 121 | SS_Z6_R20_T150_C2 | gas1_ics1 | 10 | -296 | 5,624 [5,075, 6,124] | 115.5 | 0.93 | no | 0.89 | 1,916 / 3,292 |
| 122 | SO_Z8_R30_T150_C5 | gas1_ics1 | 10 | -295 | 5,849 [5,347, 6,391] | 129.2 | 0.87 | no | 0.89 | 1,722 / 3,384 |
| 123 | SY_Z8_R30_T150_C2 | gas1_ics1 | 10 | -290 | 4,781 [4,240, 5,332] | 77.2 | 0.93 | no | 0.91 | 1,318 / 3,109 |
| 124 | SO_Z8_R20_T100000_C2 | gas5_ics3 | 16 | -290 | 6,243 [5,710, 6,750] | 15.5 | 0.69 | no | 0.65 | 2,329 / 1,993 |
| 125 | SS_Z6_R20_T100000_C5 | gas5_ics3 | 16 | -289 | 5,849 [5,349, 6,415] | 22.8 | 0.88 | no | 0.65 | 2,782 / 2,392 |
| 126 | SY_Z10_R20_T100000_C2 | gas5_ics3 | 16 | -289 | 5,737 [4,831, 6,594] | 11.7 | 0.64 | no | 0.66 | 1,910 / 1,781 |
| 127 | SL_Z6_R20_T100000_C5 | gas5_ics3 | 16 | -288 | 5,737 [4,913, 6,602] | 17.1 | 0.78 | no | 0.65 | 2,364 / 2,114 |
| 128 | SL_Z10_R30_T150_C5 | gas1_ics1 | 10 | -287 | 4,275 [3,688, 4,825] | 56.5 | 0.85 | yes | 0.92 | 588 / 3,056 |
| 129 | SS_Z6_R30_T100000_C2 | gas5_ics3 | 16 | -287 | 6,018 [5,398, 6,737] | 21.9 | 0.85 | no | 0.65 | 2,737 / 2,406 |
| 130 | SS_Z10_R20_T100000_C2 | gas5_ics3 | 16 | -287 | 5,849 [5,222, 6,628] | 13.0 | 0.70 | no | 0.66 | 2,209 / 1,873 |
| 131 | SS_Z6_R20_T100000_C2 | gas5_ics3 | 16 | -282 | 6,018 [5,398, 6,737] | 22.0 | 0.87 | no | 0.65 | 2,797 / 2,447 |
| 132 | SY_Z6_R20_T150_C2 | gas1_ics1 | 10 | -281 | 4,893 [4,371, 5,453] | 83.0 | 0.93 | no | 0.91 | 1,423 / 3,147 |
| 133 | SO_Z8_R30_T100000_C5 | gas5_ics3 | 16 | -281 | 6,187 [5,662, 6,702] | 16.1 | 0.68 | no | 0.65 | 2,269 / 1,910 |
| 134 | SO_Z8_R20_T150_C5 | gas1_ics1 | 10 | -279 | 5,849 [5,344, 6,388] | 129.0 | 0.87 | no | 0.89 | 1,722 / 3,358 |
| 135 | SL_Z10_R20_T100000_C2 | gas5_ics3 | 16 | -275 | 5,793 [4,866, 6,694] | 8.7 | 0.59 | no | 0.66 | 1,776 / 1,623 |
| 136 | SY_Z6_R30_T100000_C2 | gas5_ics3 | 16 | -273 | 5,906 [5,076, 6,736] | 18.4 | 0.83 | no | 0.65 | 2,603 / 2,328 |
| 137 | SY_Z10_R30_T150_C2 | gas1_ics1 | 10 | -258 | 4,724 [4,159, 5,260] | 73.7 | 0.92 | no | 0.91 | 1,259 / 3,096 |
| 138 | SL_Z8_R30_T150_C5 | gas1_ics1 | 10 | -257 | 4,331 [3,754, 4,881] | 59.2 | 0.83 | yes | 0.92 | 529 / 3,056 |
| 139 | SY_Z6_R20_T100000_C2 | gas5_ics3 | 16 | -253 | 5,906 [5,081, 6,737] | 18.5 | 0.84 | no | 0.65 | 2,707 / 2,260 |
| 140 | SO_Z10_R20_T100000_C2 | gas5_ics3 | 16 | -251 | 6,187 [5,660, 6,704] | 13.0 | 0.64 | no | 0.66 | 2,117 / 1,833 |
| 141 | SO_Z10_R20_T100000_C5 | gas5_ics3 | 16 | -250 | 6,131 [5,623, 6,751] | 13.8 | 0.64 | no | 0.65 | 2,102 / 1,819 |
| 142 | SO_Z6_R30_T150_C5 | gas1_ics1 | 10 | -245 | 5,962 [5,433, 6,476] | 133.0 | 0.87 | no | 0.88 | 1,812 / 3,397 |
| 143 | SL_Z6_R20_T100000_C2 | gas5_ics3 | 16 | -235 | 5,962 [5,100, 6,771] | 16.5 | 0.79 | no | 0.65 | 2,409 / 2,291 |
| 144 | SO_Z6_R20_T150_C5 | gas1_ics1 | 10 | -231 | 5,962 [5,428, 6,472] | 132.7 | 0.87 | no | 0.88 | 1,812 / 3,397 |
| 145 | SL_Z6_R30_T100000_C2 | gas5_ics3 | 16 | -229 | 5,962 [5,096, 6,772] | 16.4 | 0.78 | no | 0.65 | 2,394 / 2,264 |
| 146 | SO_Z6_R30_T100000_C5 | gas5_ics3 | 16 | -218 | 6,243 [5,731, 6,769] | 22.3 | 0.78 | no | 0.65 | 2,571 / 2,324 |
| 147 | SY_Z8_R20_T150_C2 | gas1_ics1 | 10 | -214 | 4,781 [4,227, 5,318] | 76.8 | 0.92 | no | 0.91 | 1,303 / 3,109 |
| 148 | SL_Z10_R20_T150_C5 | gas1_ics1 | 10 | -214 | 4,275 [3,686, 4,822] | 56.6 | 0.89 | no | 0.92 | 767 / 3,043 |
| 149 | SO_Z6_R30_T100000_C2 | gas5_ics3 | 16 | -213 | 6,299 [5,778, 6,882] | 21.6 | 0.79 | no | 0.65 | 2,601 / 2,406 |
| 150 | SL_Z8_R20_T150_C5 | gas1_ics1 | 10 | -209 | 4,331 [3,751, 4,877] | 59.2 | 0.86 | no | 0.92 | 648 / 3,056 |
| 151 | SO_Z6_R20_T100000_C5 | gas5_ics3 | 16 | -205 | 6,243 [5,730, 6,769] | 22.4 | 0.79 | no | 0.65 | 2,571 / 2,337 |
| 152 | SO_Z6_R20_T100000_C2 | gas5_ics3 | 16 | -198 | 6,299 [5,777, 6,815] | 21.6 | 0.80 | no | 0.65 | 2,616 / 2,420 |
| 153 | SY_Z10_R20_T150_C2 | gas1_ics1 | 10 | -189 | 4,668 [4,146, 5,246] | 73.3 | 0.93 | no | 0.91 | 1,244 / 3,083 |
| 154 | SO_Z8_R20_T100000_C5 | gas5_ics3 | 16 | -178 | 6,187 [5,658, 6,699] | 16.3 | 0.69 | no | 0.65 | 2,313 / 1,965 |
| 155 | SY_Z4_R30_T150_C2 | gas1_ics1 | 10 | -107 | 5,062 [4,555, 5,625] | 91.3 | 0.95 | no | 0.90 | 1,602 / 3,186 |
| 156 | SY_Z4_R20_T150_C2 | gas1_ics1 | 10 | -101 | 5,062 [4,555, 5,625] | 91.3 | 0.95 | no | 0.90 | 1,602 / 3,186 |
| 157 | SS_Z4_R30_T150_C2 | gas1_ics1 | 10 | -98 | 5,737 [5,212, 6,259] | 121.7 | 0.94 | no | 0.89 | 2,065 / 3,317 |
| 158 | SL_Z4_R30_T150_C5 | gas1_ics1 | 10 | -98 | 4,668 [4,114, 5,205] | 73.5 | 0.91 | no | 0.91 | 1,095 / 3,160 |
| 159 | SS_Z4_R30_T100000_C5 | gas5_ics3 | 16 | -94 | 5,962 [5,451, 6,500] | 40.3 | 0.93 | no | 0.64 | 2,768 / 2,777 |
| 160 | SS_Z4_R20_T150_C2 | gas1_ics1 | 10 | -94 | 5,737 [5,211, 6,258] | 121.7 | 0.94 | no | 0.89 | 2,050 / 3,317 |
| 161 | SL_Z4_R20_T150_C5 | gas1_ics1 | 10 | -92 | 4,668 [4,114, 5,205] | 73.5 | 0.91 | no | 0.91 | 1,095 / 3,160 |
| 162 | SS_Z4_R30_T100000_C2 | gas5_ics3 | 16 | -91 | 6,131 [5,510, 6,776] | 39.3 | 0.91 | no | 0.64 | 2,768 / 2,789 |
| 163 | SY_Z4_R30_T100000_C2 | gas5_ics3 | 16 | -87 | 6,074 [5,321, 6,828] | 35.2 | 0.90 | no | 0.64 | 2,559 / 2,925 |
| 164 | SY_Z4_R20_T100000_C2 | gas5_ics3 | 16 | -80 | 6,074 [5,321, 6,827] | 35.2 | 0.89 | no | 0.64 | 2,723 / 2,710 |
| 165 | SS_Z4_R20_T100000_C2 | gas5_ics3 | 16 | -60 | 6,131 [5,510, 6,776] | 39.4 | 0.90 | no | 0.64 | 2,738 / 2,801 |
| 166 | SL_Z4_R30_T100000_C2 | gas5_ics3 | 16 | -1 | 6,131 [5,347, 6,733] | 32.8 | 0.81 | no | 0.65 | 2,380 / 2,608 |
| 167 | SL_Z6_R30_T150_C2 | gas1_ics1 | 10 | +0 | 4,556 [3,990, 5,103] | 67.2 | 0.90 | no | 0.91 | 991 / 3,095 |
| 168 | SL_Z4_R20_T100000_C2 | gas5_ics3 | 16 | +2 | 6,074 [5,345, 6,734] | 32.8 | 0.82 | no | 0.65 | 2,380 / 2,583 |
| 169 | SO_Z4_R30_T150_C5 | gas1_ics1 | 10 | +13 | 6,018 [5,505, 6,550] | 135.8 | 0.89 | no | 0.88 | 1,931 / 3,422 |
| 170 | SO_Z4_R20_T150_C5 | gas1_ics1 | 10 | +15 | 6,018 [5,504, 6,550] | 135.8 | 0.89 | no | 0.88 | 1,931 / 3,422 |
| 171 | SL_Z6_R20_T150_C2 | gas1_ics1 | 10 | +17 | 4,556 [3,992, 5,103] | 67.3 | 0.92 | no | 0.91 | 1,080 / 3,095 |
| 172 | SY_Z8_R30_T100000_C5 | gas1_ics1 | 10 | +19 | 4,837 [4,262, 5,361] | 77.4 | 0.94 | no | 0.91 | 1,425 / 3,110 |
| 173 | SL_Z8_R30_T150_C2 | gas1_ics1 | 10 | +27 | 4,387 [3,851, 4,980] | 61.7 | 0.91 | no | 0.92 | 946 / 3,056 |
| 174 | SO_Z10_R30_T150_C2 | gas1_ics1 | 10 | +34 | 5,906 [5,388, 6,428] | 131.4 | 0.88 | no | 0.89 | 1,856 / 3,345 |
| 175 | SL_Z10_R30_T150_C2 | gas1_ics1 | 10 | +45 | 4,331 [3,773, 4,912] | 58.5 | 0.95 | no | 0.92 | 1,065 / 3,030 |
| 176 | SO_Z10_R20_T150_C2 | gas1_ics1 | 10 | +52 | 5,906 [5,379, 6,420] | 130.9 | 0.88 | no | 0.89 | 1,841 / 3,345 |
| 177 | SO_Z8_R30_T150_C2 | gas1_ics1 | 10 | +69 | 5,962 [5,426, 6,468] | 132.9 | 0.88 | no | 0.89 | 1,901 / 3,358 |
| 178 | SO_Z4_R30_T100000_C2 | gas5_ics3 | 16 | +78 | 6,412 [5,854, 6,972] | 40.0 | 0.88 | no | 0.64 | 2,723 / 2,894 |
| 179 | SO_Z8_R20_T150_C2 | gas1_ics1 | 10 | +81 | 5,962 [5,422, 6,464] | 132.7 | 0.88 | no | 0.89 | 1,901 / 3,358 |
| 180 | SO_Z4_R20_T100000_C2 | gas5_ics3 | 16 | +82 | 6,356 [5,851, 6,974] | 40.0 | 0.88 | no | 0.64 | 2,693 / 2,868 |
| 181 | SY_Z8_R30_T100000_C2 | gas1_ics1 | 10 | +94 | 4,837 [4,284, 5,383] | 77.6 | 0.93 | no | 0.91 | 1,395 / 3,084 |
| 182 | SS_Z4_R20_T100000_C5 | gas5_ics3 | 16 | +96 | 5,962 [5,435, 6,485] | 40.5 | 0.93 | no | 0.64 | 2,723 / 2,814 |
| 183 | SL_Z8_R30_T100000_C5 | gas1_ics1 | 10 | +98 | 4,443 [3,882, 5,017] | 62.0 | 0.98 | no | 0.92 | 1,290 / 3,045 |
| 184 | SY_Z4_R30_T100000_C5 | gas5_ics3 | 16 | +99 | 5,906 [5,152, 6,658] | 36.2 | 0.94 | no | 0.64 | 2,664 / 2,875 |
| 185 | SL_Z8_R20_T150_C2 | gas1_ics1 | 10 | +105 | 4,387 [3,837, 4,965] | 61.3 | 0.92 | no | 0.92 | 975 / 3,056 |
| 186 | SY_Z4_R20_T100000_C5 | gas5_ics3 | 16 | +108 | 5,906 [5,153, 6,658] | 36.3 | 0.91 | no | 0.64 | 2,768 / 2,583 |
| 187 | SL_Z10_R20_T150_C2 | gas1_ics1 | 10 | +117 | 4,331 [3,764, 4,902] | 58.4 | 0.94 | no | 0.92 | 1,035 / 3,030 |
| 188 | SO_Z6_R30_T150_C2 | gas1_ics1 | 10 | +131 | 6,018 [5,514, 6,556] | 136.9 | 0.90 | no | 0.88 | 1,991 / 3,397 |
| 189 | SO_Z6_R20_T150_C2 | gas1_ics1 | 10 | +134 | 6,018 [5,513, 6,555] | 136.8 | 0.90 | no | 0.88 | 1,991 / 3,397 |
| 190 | SY_Z8_R20_T100000_C2 | gas1_ics1 | 10 | +143 | 4,837 [4,281, 5,379] | 77.7 | 0.93 | no | 0.91 | 1,394 / 3,084 |
| 191 | SY_Z10_R30_T100000_C5 | gas1_ics1 | 10 | +180 | 4,724 [4,174, 5,279] | 73.6 | 0.94 | no | 0.91 | 1,350 / 3,098 |
| 192 | SL_Z4_R30_T100000_C5 | gas5_ics3 | 16 | +188 | 5,906 [5,160, 6,617] | 33.7 | 0.86 | no | 0.65 | 2,455 / 2,646 |
| 193 | SL_Z4_R20_T100000_C5 | gas5_ics3 | 16 | +194 | 5,906 [5,160, 6,617] | 33.7 | 0.87 | no | 0.65 | 2,455 / 2,672 |
| 194 | SY_Z8_R20_T100000_C5 | gas1_ics1 | 10 | +197 | 4,781 [4,242, 5,339] | 76.8 | 0.94 | no | 0.91 | 1,395 / 3,110 |
| 195 | SL_Z8_R30_T100000_C2 | gas1_ics1 | 10 | +205 | 4,499 [3,909, 5,047] | 62.3 | 0.99 | no | 0.92 | 1,439 / 3,032 |
| 196 | SY_Z6_R30_T100000_C5 | gas1_ics1 | 10 | +212 | 4,893 [4,365, 5,452] | 82.2 | 0.95 | no | 0.91 | 1,499 / 3,149 |
| 197 | SS_Z8_R30_T100000_C5 | gas1_ics1 | 10 | +218 | 5,512 [4,962, 6,017] | 109.9 | 0.93 | no | 0.89 | 1,888 / 3,240 |
| 198 | SY_Z10_R30_T100000_C2 | gas1_ics1 | 10 | +222 | 4,724 [4,194, 5,302] | 73.8 | 0.94 | no | 0.91 | 1,364 / 3,071 |
| 199 | SL_Z4_R30_T150_C2 | gas1_ics1 | 10 | +228 | 4,724 [4,188, 5,280] | 75.6 | 0.91 | no | 0.91 | 1,125 / 3,160 |
| 200 | SL_Z4_R20_T150_C2 | gas1_ics1 | 10 | +232 | 4,724 [4,188, 5,279] | 75.6 | 0.94 | no | 0.91 | 1,274 / 3,160 |
| 201 | SY_Z6_R20_T100000_C5 | gas1_ics1 | 10 | +237 | 4,893 [4,364, 5,450] | 82.2 | 0.95 | no | 0.91 | 1,499 / 3,136 |
| 202 | SS_Z8_R30_T100000_C2 | gas1_ics1 | 10 | +242 | 5,512 [4,990, 6,044] | 110.7 | 0.94 | no | 0.89 | 1,947 / 3,240 |
| 203 | SL_Z8_R20_T100000_C2 | gas1_ics1 | 10 | +251 | 4,499 [3,909, 5,046] | 62.5 | 0.96 | no | 0.92 | 1,305 / 3,032 |
| 204 | SY_Z6_R30_T100000_C2 | gas1_ics1 | 10 | +271 | 4,949 [4,388, 5,477] | 82.6 | 0.97 | no | 0.91 | 1,648 / 3,136 |
| 205 | SL_Z10_R30_T100000_C5 | gas1_ics1 | 10 | +274 | 4,387 [3,793, 4,939] | 58.6 | 0.93 | no | 0.92 | 1,081 / 3,019 |
| 206 | SS_Z8_R20_T100000_C2 | gas1_ics1 | 10 | +277 | 5,512 [4,980, 6,034] | 110.3 | 0.94 | no | 0.89 | 1,947 / 3,226 |
| 207 | SY_Z10_R20_T100000_C5 | gas1_ics1 | 10 | +285 | 4,724 [4,162, 5,266] | 73.4 | 0.94 | no | 0.91 | 1,335 / 3,084 |
| 208 | SO_Z4_R30_T100000_C5 | gas5_ics3 | 16 | +291 | 6,299 [5,764, 6,805] | 41.1 | 0.89 | no | 0.64 | 2,738 / 2,843 |
| 209 | SY_Z6_R20_T100000_C2 | gas1_ics1 | 10 | +295 | 4,949 [4,387, 5,475] | 82.6 | 0.95 | no | 0.91 | 1,544 / 3,136 |
| 210 | SO_Z4_R20_T100000_C5 | gas5_ics3 | 16 | +298 | 6,243 [5,762, 6,806] | 41.1 | 0.91 | no | 0.64 | 2,798 / 2,856 |
| 211 | SL_Z6_R30_T100000_C5 | gas1_ics1 | 10 | +304 | 4,556 [3,985, 5,104] | 66.5 | 0.99 | no | 0.91 | 1,424 / 3,084 |
| 212 | SS_Z8_R20_T100000_C5 | gas1_ics1 | 10 | +316 | 5,456 [4,944, 5,999] | 109.0 | 0.94 | no | 0.89 | 1,873 / 3,240 |
| 213 | SS_Z10_R30_T100000_C2 | gas1_ics1 | 10 | +323 | 5,456 [4,911, 5,968] | 106.9 | 0.93 | no | 0.90 | 1,887 / 3,214 |
| 214 | SL_Z8_R20_T100000_C5 | gas1_ics1 | 10 | +324 | 4,443 [3,858, 4,992] | 61.4 | 0.97 | no | 0.92 | 1,275 / 3,045 |
| 215 | SL_Z6_R20_T100000_C5 | gas1_ics1 | 10 | +328 | 4,556 [3,984, 5,103] | 66.5 | 0.87 | no | 0.91 | 888 / 3,084 |
| 216 | SS_Z10_R30_T100000_C5 | gas1_ics1 | 10 | +330 | 5,399 [4,886, 5,943] | 106.2 | 0.94 | no | 0.90 | 1,828 / 3,227 |
| 217 | SY_Z10_R20_T100000_C2 | gas1_ics1 | 10 | +342 | 4,724 [4,179, 5,286] | 73.5 | 0.94 | no | 0.91 | 1,379 / 3,058 |
| 218 | SL_Z10_R30_T100000_C2 | gas1_ics1 | 10 | +355 | 4,387 [3,819, 4,967] | 58.8 | 0.99 | no | 0.92 | 1,319 / 3,005 |
| 219 | SO_Z4_R30_T150_C2 | gas1_ics1 | 10 | +374 | 6,131 [5,593, 6,637] | 140.0 | 0.90 | no | 0.88 | 2,110 / 3,409 |
| 220 | SO_Z4_R20_T150_C2 | gas1_ics1 | 10 | +375 | 6,131 [5,592, 6,636] | 140.0 | 0.90 | no | 0.88 | 2,110 / 3,409 |
| 221 | SL_Z10_R20_T100000_C5 | gas1_ics1 | 10 | +379 | 4,331 [3,786, 4,930] | 58.6 | 0.95 | no | 0.92 | 1,096 / 3,019 |
| 222 | SL_Z6_R30_T100000_C2 | gas1_ics1 | 10 | +405 | 4,556 [4,011, 5,132] | 66.8 | 0.97 | no | 0.91 | 1,334 / 3,070 |
| 223 | SS_Z10_R20_T100000_C5 | gas1_ics1 | 10 | +407 | 5,399 [4,869, 5,926] | 105.5 | 0.93 | no | 0.90 | 1,798 / 3,227 |
| 224 | SS_Z6_R30_T100000_C5 | gas1_ics1 | 10 | +410 | 5,568 [5,050, 6,102] | 114.0 | 0.94 | no | 0.89 | 1,977 / 3,265 |
| 225 | SS_Z10_R20_T100000_C2 | gas1_ics1 | 10 | +416 | 5,399 [4,890, 5,948] | 106.0 | 0.94 | no | 0.90 | 1,857 / 3,201 |
| 226 | SS_Z6_R30_T100000_C2 | gas1_ics1 | 10 | +419 | 5,624 [5,080, 6,132] | 115.0 | 0.94 | no | 0.89 | 2,036 / 3,265 |
| 227 | SY_Z4_R30_T100000_C2 | gas1_ics1 | 10 | +421 | 5,118 [4,609, 5,682] | 92.7 | 0.97 | no | 0.90 | 1,768 / 3,187 |
| 228 | SL_Z6_R20_T100000_C2 | gas1_ics1 | 10 | +427 | 4,556 [4,011, 5,132] | 66.9 | 0.97 | no | 0.91 | 1,334 / 3,070 |
| 229 | SY_Z4_R20_T100000_C2 | gas1_ics1 | 10 | +429 | 5,118 [4,609, 5,682] | 92.7 | 0.97 | no | 0.90 | 1,768 / 3,187 |
| 230 | SS_Z6_R20_T100000_C5 | gas1_ics1 | 10 | +429 | 5,568 [5,045, 6,097] | 113.8 | 0.94 | no | 0.89 | 1,962 / 3,265 |
| 231 | SS_Z6_R20_T100000_C2 | gas1_ics1 | 10 | +436 | 5,624 [5,075, 6,127] | 114.8 | 0.94 | no | 0.89 | 2,036 / 3,265 |
| 232 | SL_Z10_R20_T100000_C2 | gas1_ics1 | 10 | +470 | 4,387 [3,809, 4,956] | 58.8 | 1.01 | no | 0.92 | 1,453 / 2,993 |
| 233 | SS_Z4_R30_T100000_C5 | gas1_ics1 | 10 | +571 | 5,737 [5,212, 6,260] | 121.6 | 0.95 | no | 0.89 | 2,141 / 3,317 |
| 234 | SS_Z4_R30_T100000_C2 | gas1_ics1 | 10 | +594 | 5,793 [5,251, 6,299] | 123.1 | 0.95 | no | 0.89 | 2,216 / 3,303 |
| 235 | SL_Z4_R30_T100000_C2 | gas1_ics1 | 10 | +596 | 4,781 [4,246, 5,342] | 76.8 | 0.90 | no | 0.91 | 1,171 / 3,135 |
| 236 | SL_Z4_R20_T100000_C2 | gas1_ics1 | 10 | +599 | 4,781 [4,246, 5,342] | 76.8 | 0.99 | no | 0.91 | 1,618 / 3,135 |
| 237 | SS_Z4_R20_T100000_C2 | gas1_ics1 | 10 | +601 | 5,793 [5,249, 6,298] | 123.0 | 0.95 | no | 0.89 | 2,216 / 3,303 |
| 238 | SO_Z8_R30_T100000_C2 | gas1_ics1 | 10 | +602 | 5,962 [5,434, 6,476] | 132.9 | 0.91 | no | 0.88 | 2,067 / 3,331 |
| 239 | SY_Z4_R30_T100000_C5 | gas1_ics1 | 10 | +634 | 5,062 [4,550, 5,623] | 90.7 | 0.96 | no | 0.90 | 1,664 / 3,187 |
| 240 | SY_Z4_R20_T100000_C5 | gas1_ics1 | 10 | +644 | 5,062 [4,550, 5,623] | 90.7 | 0.96 | no | 0.90 | 1,664 / 3,187 |
| 241 | SO_Z10_R30_T100000_C5 | gas1_ics1 | 10 | +658 | 5,849 [5,343, 6,388] | 128.8 | 0.90 | no | 0.89 | 1,948 / 3,318 |
| 242 | SO_Z10_R30_T100000_C2 | gas1_ics1 | 10 | +665 | 5,906 [5,377, 6,420] | 130.2 | 0.90 | no | 0.89 | 2,007 / 3,318 |
| 243 | SO_Z8_R20_T100000_C2 | gas1_ics1 | 10 | +674 | 5,962 [5,420, 6,464] | 132.1 | 0.90 | no | 0.88 | 2,067 / 3,318 |
| 244 | SO_Z8_R30_T100000_C5 | gas1_ics1 | 10 | +694 | 5,906 [5,380, 6,424] | 130.3 | 0.90 | no | 0.88 | 1,992 / 3,331 |
| 245 | SO_Z10_R20_T100000_C5 | gas1_ics1 | 10 | +702 | 5,849 [5,333, 6,377] | 128.2 | 0.90 | no | 0.89 | 1,948 / 3,318 |
| 246 | SO_Z10_R20_T100000_C2 | gas1_ics1 | 10 | +707 | 5,906 [5,366, 6,410] | 129.6 | 0.90 | no | 0.89 | 2,007 / 3,292 |
| 247 | SS_Z4_R20_T100000_C5 | gas1_ics1 | 10 | +761 | 5,737 [5,189, 6,239] | 120.3 | 0.95 | no | 0.89 | 2,126 / 3,317 |
| 248 | SO_Z6_R30_T100000_C5 | gas1_ics1 | 10 | +771 | 5,962 [5,457, 6,503] | 133.7 | 0.91 | no | 0.88 | 2,067 / 3,370 |
| 249 | SO_Z8_R20_T100000_C5 | gas1_ics1 | 10 | +779 | 5,906 [5,368, 6,414] | 129.6 | 0.90 | no | 0.88 | 1,978 / 3,331 |
| 250 | SO_Z6_R20_T100000_C5 | gas1_ics1 | 10 | +779 | 5,962 [5,455, 6,501] | 133.6 | 0.91 | no | 0.88 | 2,067 / 3,370 |
| 251 | SO_Z6_R30_T100000_C2 | gas1_ics1 | 10 | +784 | 6,018 [5,500, 6,544] | 135.6 | 0.91 | no | 0.88 | 2,141 / 3,343 |
| 252 | SO_Z6_R20_T100000_C2 | gas1_ics1 | 10 | +790 | 6,018 [5,498, 6,542] | 135.5 | 0.91 | no | 0.88 | 2,141 / 3,343 |
| 253 | SL_Z4_R30_T100000_C5 | gas1_ics1 | 10 | +802 | 4,724 [4,183, 5,279] | 75.0 | 0.92 | no | 0.91 | 1,187 / 3,148 |
| 254 | SL_Z4_R20_T100000_C5 | gas1_ics1 | 10 | +809 | 4,724 [4,183, 5,279] | 75.0 | 0.92 | no | 0.91 | 1,201 / 3,148 |
| 255 | SO_Z4_R30_T100000_C2 | gas1_ics1 | 10 | +1,026 | 6,131 [5,617, 6,662] | 140.8 | 0.93 | no | 0.88 | 2,290 / 3,381 |
| 256 | SO_Z4_R20_T100000_C2 | gas1_ics1 | 10 | +1,029 | 6,131 [5,617, 6,662] | 140.8 | 0.93 | no | 0.88 | 2,290 / 3,381 |
| 257 | SO_Z4_R30_T100000_C5 | gas1_ics1 | 10 | +1,196 | 6,074 [5,544, 6,591] | 137.1 | 0.92 | no | 0.88 | 2,201 / 3,381 |
| 258 | SO_Z4_R20_T100000_C5 | gas1_ics1 | 10 | +1,200 | 6,074 [5,543, 6,591] | 137.1 | 0.92 | no | 0.88 | 2,201 / 3,381 |
