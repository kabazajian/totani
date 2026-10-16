# Stage C — 2026-09-20 (no_residual)

Consensus selection (median Δ(−2lnL) across 217 model rows; BIC threshold p·ln 1080 = 63 per template): **bubble_flat → erosita_outer_NS → loopI_haslam_raw → sun_moon_ic → loopIV**

| step | chosen | median Δ | min | max | threshold | runner-up (median Δ) |
|---|---|---|---|---|---|---|
| 1 | bubble_flat | 19,194 | 12,490 | 26,944 | 63 | bubble_structured (18,996) |
| 2 | erosita_outer_NS | 10,308 | 8,300 | 12,210 | 126 | loopI_shells+nps (8,591) |
| 3 | loopI_haslam_raw | 984 | 757 | 1,907 | 63 | loopI_shells+nps (631) |
| 4 | sun_moon_ic | 337 | 140 | 743 | 63 | loopIV (219) |
| 5 | loopIV | 94 | 20 | 207 | 63 | loopII (1) |
| 6 | (stop) loopII | 1 | -28 | 28 | 63 | loopIII (-0) |

## Each model's own path

| model | variant | path |
|---|---|---|
| SL_Z10_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells |
| SL_Z10_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SL_Z10_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells |
| SL_Z10_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z10_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_outer_NS |
| SL_Z10_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z10_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_outer_NS → loopII |
| SL_Z10_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z10_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells |
| SL_Z10_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z10_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells |
| SL_Z10_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z10_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_outer_NS |
| SL_Z10_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z10_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_outer_NS |
| SL_Z10_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z4_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SL_Z4_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SL_Z4_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SL_Z4_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z4_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SL_Z4_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z4_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_outer_NS |
| SL_Z4_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z4_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SL_Z4_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SL_Z4_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SL_Z4_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z4_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SL_Z4_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z4_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_outer_NS |
| SL_Z4_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SL_Z6_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells |
| SL_Z6_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z6_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells |
| SL_Z6_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z6_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_outer_NS |
| SL_Z6_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z6_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_outer_NS → loopII |
| SL_Z6_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z6_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells |
| SL_Z6_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z6_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells |
| SL_Z6_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z6_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_outer_NS |
| SL_Z6_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z6_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_outer_NS → loopII |
| SL_Z6_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z8_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells |
| SL_Z8_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z8_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells |
| SL_Z8_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z8_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells |
| SL_Z8_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z8_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_outer_NS → loopII |
| SL_Z8_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z8_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells |
| SL_Z8_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z8_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells |
| SL_Z8_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z8_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_outer_NS |
| SL_Z8_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z8_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_outer_NS → loopII |
| SL_Z8_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SO_Z4_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → sun_moon_ic |
| SO_Z4_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SO_Z4_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SO_Z4_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SO_Z4_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SO_Z4_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SO_Z4_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SO_Z4_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SO_Z4_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → sun_moon_ic |
| SO_Z4_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SO_Z4_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SO_Z4_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SO_Z4_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SO_Z4_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SO_Z4_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SO_Z4_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SO_Z6_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SO_Z6_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → loopII |
| SO_Z6_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SO_Z6_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SO_Z6_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells |
| SO_Z6_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SO_Z6_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SO_Z6_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
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
| SS_Z10_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SS_Z10_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
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
| SY_Z10_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells |
| SY_Z10_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z10_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells |
| SY_Z10_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z10_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SY_Z10_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z10_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_outer_NS |
| SY_Z10_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z10_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells |
| SY_Z10_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z10_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells |
| SY_Z10_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z10_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_outer_NS |
| SY_Z10_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z10_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_outer_NS |
| SY_Z10_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z4_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells |
| SY_Z4_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z4_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SY_Z4_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z4_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells |
| SY_Z4_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z4_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_outer_NS |
| SY_Z4_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z4_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells |
| SY_Z4_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z4_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SY_Z4_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z4_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells |
| SY_Z4_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z4_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_outer_NS |
| SY_Z4_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z6_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells |
| SY_Z6_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z6_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells |
| SY_Z6_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z6_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SY_Z6_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z6_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_outer_NS → loopII |
| SY_Z6_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z6_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells |
| SY_Z6_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z6_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells |
| SY_Z6_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z6_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SY_Z6_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z6_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_outer_NS |
| SY_Z6_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z8_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells |
| SY_Z8_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z8_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells |
| SY_Z8_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z8_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SY_Z8_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z8_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_outer_NS |
| SY_Z8_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z8_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells |
| SY_Z8_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z8_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells |
| SY_Z8_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z8_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_outer_NS |
| SY_Z8_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z8_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_outer_NS |
| SY_Z8_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| a2020 | - | bubble_totani → loopI_shells+nps |

## Ranking with the consensus set (bubble_flat, erosita_outer_NS, loopI_haslam_raw, sun_moon_ic, loopIV)

| rank | model | variant | n | Δ(−2lnL) vs Totani baseline | halo @ 20.8 | TS @ 20.8 | quad Σ/joint | any neg | inj min | N / S |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | gll_iem_v07 | - | 9 | -4,200 | 2,025 [1,429, 2,578] | 12.3 | 0.84 | yes | 0.90 | 682 / 1,020 |
| 2 | a2020 | - | 22 | -1,851 | 7,255 [6,669, 7,791] | 48.1 | 1.01 | no | 0.98 | 3,400 / 3,920 |
| 3 | SS_Z10_R30_T150_C5 | gas5_ics3 | 16 | -1,352 | 5,962 [5,219, 6,707] | 11.3 | 0.63 | no | 0.66 | 2,194 / 1,580 |
| 4 | SY_Z10_R30_T150_C5 | gas5_ics3 | 16 | -1,346 | 5,849 [4,975, 6,700] | 9.2 | 0.60 | no | 0.66 | 2,076 / 1,412 |
| 5 | SL_Z10_R30_T150_C5 | gas5_ics3 | 16 | -1,328 | 5,906 [5,028, 6,717] | 8.7 | 0.53 | no | 0.66 | 1,868 / 1,280 |
| 6 | SY_Z10_R20_T150_C5 | gas5_ics3 | 16 | -1,304 | 5,849 [4,982, 6,696] | 9.5 | 0.60 | no | 0.66 | 2,076 / 1,452 |
| 7 | SS_Z10_R20_T150_C5 | gas5_ics3 | 16 | -1,303 | 5,962 [5,215, 6,697] | 11.5 | 0.63 | no | 0.66 | 2,194 / 1,570 |
| 8 | SL_Z10_R20_T150_C5 | gas5_ics3 | 16 | -1,294 | 5,906 [5,034, 6,708] | 8.9 | 0.54 | no | 0.66 | 1,883 / 1,283 |
| 9 | SS_Z8_R30_T150_C5 | gas5_ics3 | 16 | -1,282 | 6,018 [5,270, 6,718] | 14.0 | 0.67 | no | 0.66 | 2,181 / 1,863 |
| 10 | SY_Z8_R30_T150_C5 | gas5_ics3 | 16 | -1,273 | 5,849 [5,040, 6,716] | 11.8 | 0.68 | no | 0.65 | 2,182 / 1,784 |
| 11 | SY_Z6_R30_T150_C5 | gas5_ics3 | 16 | -1,270 | 6,018 [5,193, 6,779] | 17.2 | 0.77 | no | 0.65 | 2,409 / 2,227 |
| 12 | SS_Z6_R30_T150_C5 | gas5_ics3 | 16 | -1,269 | 6,074 [5,387, 6,770] | 19.7 | 0.77 | no | 0.65 | 2,409 / 2,291 |
| 13 | SY_Z6_R20_T150_C5 | gas5_ics3 | 16 | -1,259 | 6,018 [5,196, 6,777] | 17.2 | 0.77 | no | 0.65 | 2,409 / 2,227 |
| 14 | SS_Z6_R20_T150_C5 | gas5_ics3 | 16 | -1,253 | 6,074 [5,388, 6,769] | 19.8 | 0.78 | no | 0.65 | 2,424 / 2,305 |
| 15 | SL_Z8_R30_T150_C5 | gas5_ics3 | 16 | -1,250 | 5,962 [5,072, 6,695] | 10.9 | 0.62 | no | 0.66 | 2,033 / 1,678 |
| 16 | SS_Z8_R20_T150_C5 | gas5_ics3 | 16 | -1,248 | 6,018 [5,272, 6,714] | 14.2 | 0.66 | no | 0.65 | 2,107 / 1,878 |
| 17 | SY_Z8_R20_T150_C5 | gas5_ics3 | 16 | -1,247 | 5,906 [5,044, 6,714] | 11.9 | 0.67 | no | 0.65 | 2,212 / 1,734 |
| 18 | SL_Z6_R30_T150_C5 | gas5_ics3 | 16 | -1,246 | 6,018 [5,217, 6,685] | 16.0 | 0.73 | no | 0.65 | 2,260 / 2,148 |
| 19 | SL_Z6_R20_T150_C5 | gas5_ics3 | 16 | -1,230 | 6,074 [5,217, 6,682] | 16.1 | 0.73 | no | 0.65 | 2,275 / 2,161 |
| 20 | SL_Z8_R20_T150_C5 | gas5_ics3 | 16 | -1,227 | 5,962 [5,075, 6,690] | 11.0 | 0.62 | no | 0.65 | 2,033 / 1,641 |
| 21 | SO_Z6_R20_T150_C5 | gas5_ics3 | 16 | -1,146 | 6,356 [5,739, 6,876] | 19.3 | 0.72 | no | 0.65 | 2,317 / 2,250 |
| 22 | SS_Z10_R30_T150_C2 | gas5_ics3 | 16 | -1,067 | 6,018 [4,535, 6,887] | 10.8 | 0.65 | no | 0.66 | 2,299 / 1,609 |
| 23 | SY_Z10_R30_T150_C2 | gas5_ics3 | 16 | -1,063 | 6,018 [5,178, 6,911] | 8.8 | 0.61 | no | 0.66 | 2,136 / 1,529 |
| 24 | SS_Z8_R30_T150_C2 | gas5_ics3 | 16 | -1,054 | 6,187 [5,443, 6,906] | 13.3 | 0.67 | no | 0.66 | 2,256 / 1,863 |
| 25 | SL_Z10_R30_T150_C2 | gas5_ics3 | 16 | -1,047 | 6,131 [5,246, 6,857] | 8.3 | 0.55 | no | 0.66 | 1,913 / 1,434 |
| 26 | SY_Z8_R30_T150_C2 | gas5_ics3 | 16 | -1,038 | 6,074 [5,241, 6,920] | 11.2 | 0.66 | no | 0.66 | 2,257 / 1,772 |
| 27 | SS_Z10_R20_T150_C2 | gas5_ics3 | 16 | -1,029 | 6,131 [5,382, 6,875] | 11.0 | 0.65 | no | 0.66 | 2,299 / 1,688 |
| 28 | SY_Z10_R20_T150_C2 | gas5_ics3 | 16 | -1,029 | 6,018 [4,320, 6,901] | 9.1 | 0.61 | no | 0.66 | 2,136 / 1,557 |
| 29 | SL_Z10_R20_T150_C2 | gas5_ics3 | 16 | -1,017 | 6,131 [5,245, 6,845] | 8.5 | 0.56 | no | 0.66 | 1,957 / 1,450 |
| 30 | SS_Z8_R20_T150_C2 | gas5_ics3 | 16 | -1,008 | 6,187 [5,435, 6,891] | 13.5 | 0.67 | no | 0.65 | 2,271 / 1,904 |
| 31 | SL_Z8_R30_T150_C2 | gas5_ics3 | 16 | -1,005 | 6,187 [5,290, 6,834] | 10.5 | 0.61 | no | 0.66 | 2,078 / 1,704 |
| 32 | SS_Z4_R30_T150_C5 | gas5_ics3 | 16 | -1,003 | 6,074 [5,442, 6,736] | 35.6 | 0.85 | no | 0.64 | 2,364 / 2,794 |
| 33 | SS_Z4_R20_T150_C5 | gas5_ics3 | 16 | -997 | 6,074 [5,443, 6,736] | 35.6 | 0.86 | no | 0.64 | 2,409 / 2,794 |
| 34 | SY_Z4_R30_T150_C5 | gas5_ics3 | 16 | -989 | 6,074 [5,330, 6,732] | 32.3 | 0.83 | no | 0.65 | 2,349 / 2,715 |
| 35 | SY_Z4_R20_T150_C5 | gas5_ics3 | 16 | -984 | 6,074 [5,331, 6,731] | 32.4 | 0.82 | no | 0.65 | 2,290 / 2,677 |
| 36 | SY_Z8_R20_T150_C2 | gas5_ics3 | 16 | -983 | 6,074 [5,233, 6,903] | 11.4 | 0.68 | no | 0.66 | 2,287 / 1,839 |
| 37 | SS_Z6_R30_T150_C2 | gas5_ics3 | 16 | -982 | 6,243 [5,532, 6,930] | 19.0 | 0.76 | no | 0.65 | 2,409 / 2,344 |
| 38 | SY_Z6_R30_T150_C2 | gas5_ics3 | 16 | -978 | 6,187 [5,364, 6,932] | 16.5 | 0.76 | no | 0.65 | 2,424 / 2,267 |
| 39 | SS_Z6_R20_T150_C2 | gas5_ics3 | 16 | -971 | 6,243 [5,533, 6,928] | 19.1 | 0.75 | no | 0.65 | 2,350 / 2,345 |
| 40 | SY_Z6_R20_T150_C2 | gas5_ics3 | 16 | -970 | 6,187 [5,367, 6,930] | 16.6 | 0.76 | no | 0.65 | 2,439 / 2,253 |
| 41 | SL_Z8_R20_T150_C2 | gas5_ics3 | 16 | -956 | 6,131 [5,280, 6,820] | 10.7 | 0.62 | no | 0.66 | 2,078 / 1,745 |
| 42 | SL_Z6_R30_T150_C2 | gas5_ics3 | 16 | -952 | 6,187 [5,405, 6,802] | 15.5 | 0.72 | no | 0.65 | 2,275 / 2,188 |
| 43 | SL_Z6_R20_T150_C2 | gas5_ics3 | 16 | -946 | 6,187 [5,406, 6,799] | 15.6 | 0.72 | no | 0.65 | 2,275 / 2,188 |
| 44 | SL_Z4_R30_T150_C5 | gas5_ics3 | 16 | -944 | 6,018 [5,344, 6,625] | 30.8 | 0.81 | no | 0.65 | 2,141 / 2,715 |
| 45 | SL_Z4_R20_T150_C5 | gas5_ics3 | 16 | -940 | 6,018 [5,504, 6,624] | 30.8 | 0.81 | no | 0.65 | 2,141 / 2,715 |
| 46 | SO_Z6_R20_T150_C2 | gas5_ics3 | 16 | -896 | 6,412 [5,864, 6,999] | 18.7 | 0.73 | no | 0.65 | 2,392 / 2,277 |
| 47 | SO_Z4_R30_T150_C5 | gas5_ics3 | 16 | -855 | 6,299 [5,801, 6,936] | 36.3 | 0.83 | no | 0.64 | 2,364 / 2,861 |
| 48 | SO_Z4_R20_T150_C5 | gas5_ics3 | 16 | -852 | 6,299 [5,801, 6,935] | 36.3 | 0.83 | no | 0.64 | 2,364 / 2,861 |
| 49 | SS_Z4_R30_T150_C2 | gas5_ics3 | 16 | -736 | 6,243 [5,577, 6,882] | 34.9 | 0.83 | no | 0.64 | 2,379 / 2,793 |
| 50 | SS_Z4_R20_T150_C2 | gas5_ics3 | 16 | -731 | 6,243 [5,577, 6,881] | 35.0 | 0.83 | no | 0.64 | 2,409 / 2,793 |
| 51 | SY_Z4_R30_T150_C2 | gas5_ics3 | 16 | -711 | 6,243 [5,488, 6,854] | 31.8 | 0.79 | no | 0.65 | 2,245 / 2,702 |
| 52 | SY_Z4_R20_T150_C2 | gas5_ics3 | 16 | -707 | 6,187 [5,487, 6,854] | 31.8 | 0.80 | no | 0.65 | 2,275 / 2,689 |
| 53 | SS_Z10_R30_T150_C5 | gas1_ics1 | 10 | -705 | 5,343 [4,821, 5,875] | 104.2 | 0.90 | no | 0.90 | 1,558 / 3,253 |
| 54 | SS_Z8_R30_T150_C5 | gas1_ics1 | 10 | -669 | 5,399 [4,882, 5,934] | 107.0 | 0.91 | no | 0.90 | 1,618 / 3,279 |
| 55 | SS_Z6_R30_T150_C5 | gas1_ics1 | 10 | -666 | 5,512 [4,999, 6,049] | 112.4 | 0.91 | no | 0.89 | 1,737 / 3,292 |
| 56 | SL_Z4_R30_T150_C2 | gas5_ics3 | 16 | -657 | 6,131 [5,622, 6,745] | 30.3 | 0.78 | no | 0.65 | 2,140 / 2,664 |
| 57 | SL_Z4_R20_T150_C2 | gas5_ics3 | 16 | -655 | 6,187 [5,622, 6,743] | 30.3 | 0.76 | no | 0.65 | 2,140 / 2,587 |
| 58 | SS_Z6_R20_T150_C5 | gas1_ics1 | 10 | -652 | 5,512 [4,996, 6,046] | 112.2 | 0.91 | no | 0.89 | 1,737 / 3,292 |
| 59 | SS_Z10_R20_T150_C5 | gas1_ics1 | 10 | -651 | 5,343 [4,802, 5,857] | 103.4 | 0.89 | no | 0.90 | 1,528 / 3,240 |
| 60 | SY_Z6_R30_T150_C5 | gas1_ics1 | 10 | -640 | 4,837 [4,294, 5,376] | 80.6 | 0.91 | no | 0.91 | 1,244 / 3,161 |
| 61 | SS_Z8_R20_T150_C5 | gas1_ics1 | 10 | -634 | 5,399 [4,873, 5,926] | 106.6 | 0.90 | no | 0.90 | 1,603 / 3,253 |
| 62 | SY_Z6_R20_T150_C5 | gas1_ics1 | 10 | -619 | 4,837 [4,294, 5,376] | 80.7 | 0.91 | no | 0.91 | 1,229 / 3,161 |
| 63 | SO_Z4_R30_T150_C2 | gas5_ics3 | 16 | -611 | 6,524 [5,914, 7,042] | 35.7 | 0.80 | no | 0.64 | 2,379 / 2,861 |
| 64 | SO_Z4_R20_T150_C2 | gas5_ics3 | 16 | -608 | 6,524 [5,913, 7,041] | 35.7 | 0.80 | no | 0.64 | 2,379 / 2,861 |
| 65 | SY_Z10_R30_T150_C5 | gas1_ics1 | 10 | -590 | 4,612 [4,070, 5,170] | 71.0 | 0.90 | no | 0.91 | 1,065 / 3,096 |
| 66 | SY_Z8_R30_T150_C5 | gas1_ics1 | 10 | -572 | 4,668 [4,147, 5,240] | 74.4 | 0.91 | no | 0.91 | 1,125 / 3,122 |
| 67 | SY_Z8_R30_T100000_C5 | gas5_ics3 | 16 | -552 | 5,568 [4,751, 6,447] | 13.0 | 0.71 | no | 0.65 | 2,150 / 1,804 |
| 68 | SL_Z8_R30_T100000_C5 | gas5_ics3 | 16 | -537 | 5,624 [4,739, 6,512] | 11.1 | 0.64 | no | 0.66 | 1,881 / 1,712 |
| 69 | SY_Z8_R20_T150_C5 | gas1_ics1 | 10 | -523 | 4,668 [4,142, 5,234] | 74.3 | 0.91 | no | 0.91 | 1,110 / 3,122 |
| 70 | SY_Z10_R20_T150_C5 | gas1_ics1 | 10 | -516 | 4,612 [4,064, 5,163] | 71.0 | 0.90 | no | 0.91 | 1,050 / 3,096 |
| 71 | SS_Z8_R30_T100000_C5 | gas5_ics3 | 16 | -515 | 5,793 [5,252, 6,330] | 16.3 | 0.75 | no | 0.65 | 2,434 / 1,923 |
| 72 | SY_Z8_R30_T100000_C2 | gas5_ics3 | 16 | -481 | 5,793 [4,912, 6,646] | 14.0 | 0.69 | no | 0.66 | 2,150 / 1,873 |
| 73 | SS_Z8_R30_T100000_C2 | gas5_ics3 | 16 | -477 | 5,906 [5,293, 6,677] | 15.5 | 0.76 | no | 0.65 | 2,509 / 1,978 |
| 74 | SL_Z8_R30_T100000_C2 | gas5_ics3 | 16 | -462 | 5,849 [4,942, 6,735] | 10.6 | 0.63 | no | 0.66 | 1,927 / 1,781 |
| 75 | SY_Z8_R20_T100000_C2 | gas5_ics3 | 16 | -444 | 5,793 [4,925, 6,650] | 12.5 | 0.70 | no | 0.66 | 2,135 / 1,900 |
| 76 | SS_Z4_R30_T150_C5 | gas1_ics1 | 10 | -439 | 5,624 [5,126, 6,175] | 118.0 | 0.93 | no | 0.89 | 1,886 / 3,331 |
| 77 | SS_Z4_R20_T150_C5 | gas1_ics1 | 10 | -433 | 5,624 [5,125, 6,174] | 117.9 | 0.93 | no | 0.89 | 1,886 / 3,331 |
| 78 | SL_Z8_R20_T100000_C2 | gas5_ics3 | 16 | -431 | 5,849 [4,952, 6,735] | 10.8 | 0.64 | no | 0.66 | 1,957 / 1,808 |
| 79 | SY_Z10_R30_T100000_C5 | gas5_ics3 | 16 | -429 | 5,512 [4,653, 6,380] | 10.6 | 0.66 | no | 0.66 | 1,954 / 1,657 |
| 80 | SY_Z4_R30_T150_C5 | gas1_ics1 | 10 | -428 | 5,006 [4,477, 5,547] | 88.6 | 0.92 | no | 0.91 | 1,394 / 3,213 |
| 81 | SS_Z8_R20_T100000_C2 | gas5_ics3 | 16 | -427 | 5,906 [5,296, 6,676] | 15.6 | 0.76 | no | 0.65 | 2,494 / 2,005 |
| 82 | SY_Z4_R20_T150_C5 | gas1_ics1 | 10 | -421 | 5,006 [4,477, 5,547] | 88.6 | 0.92 | no | 0.91 | 1,394 / 3,200 |
| 83 | SL_Z10_R30_T100000_C5 | gas5_ics3 | 16 | -414 | 5,568 [4,639, 6,459] | 8.9 | 0.56 | no | 0.66 | 1,700 / 1,424 |
| 84 | SS_Z10_R30_T100000_C5 | gas5_ics3 | 16 | -405 | 5,737 [5,176, 6,267] | 13.7 | 0.72 | no | 0.65 | 2,298 / 1,844 |
| 85 | SS_Z8_R20_T100000_C5 | gas5_ics3 | 16 | -404 | 5,793 [5,247, 6,325] | 16.5 | 0.76 | no | 0.65 | 2,419 / 1,964 |
| 86 | SS_Z10_R30_T100000_C2 | gas5_ics3 | 16 | -397 | 5,849 [5,217, 6,630] | 12.8 | 0.71 | no | 0.66 | 2,283 / 1,845 |
| 87 | SY_Z8_R20_T100000_C5 | gas5_ics3 | 16 | -397 | 5,568 [4,750, 6,438] | 14.9 | 0.72 | no | 0.65 | 2,106 / 1,886 |
| 88 | SY_Z10_R30_T100000_C2 | gas5_ics3 | 16 | -396 | 5,681 [4,814, 6,589] | 9.9 | 0.64 | no | 0.66 | 1,924 / 1,713 |
| 89 | SL_Z10_R30_T100000_C2 | gas5_ics3 | 16 | -363 | 5,793 [4,850, 6,691] | 8.4 | 0.58 | no | 0.66 | 1,775 / 1,568 |
| 90 | SL_Z8_R20_T100000_C5 | gas5_ics3 | 16 | -356 | 5,624 [4,732, 6,497] | 11.4 | 0.66 | no | 0.66 | 1,912 / 1,807 |
| 91 | SS_Z8_R30_T150_C2 | gas1_ics1 | 10 | -355 | 5,512 [4,966, 6,018] | 110.4 | 0.92 | no | 0.90 | 1,811 / 3,253 |
| 92 | SY_Z6_R30_T100000_C5 | gas5_ics3 | 16 | -347 | 5,737 [4,903, 6,541] | 19.1 | 0.84 | no | 0.65 | 2,588 / 2,259 |
| 93 | SS_Z10_R30_T150_C2 | gas1_ics1 | 10 | -345 | 5,399 [4,897, 5,951] | 107.1 | 0.92 | no | 0.90 | 1,736 / 3,240 |
| 94 | SY_Z10_R20_T100000_C5 | gas5_ics3 | 16 | -341 | 5,512 [4,665, 6,384] | 10.9 | 0.67 | no | 0.66 | 1,984 / 1,712 |
| 95 | SL_Z6_R30_T150_C5 | gas1_ics1 | 10 | -341 | 4,443 [3,909, 5,021] | 65.1 | 0.86 | no | 0.92 | 693 / 3,121 |
| 96 | SY_Z6_R20_T100000_C5 | gas5_ics3 | 16 | -337 | 5,737 [4,908, 6,543] | 19.2 | 0.86 | no | 0.65 | 2,603 / 2,341 |
| 97 | SL_Z10_R20_T100000_C5 | gas5_ics3 | 16 | -336 | 5,568 [4,654, 6,461] | 9.2 | 0.59 | no | 0.66 | 1,730 / 1,555 |
| 98 | SS_Z6_R30_T100000_C5 | gas5_ics3 | 16 | -326 | 5,906 [5,342, 6,408] | 22.7 | 0.86 | no | 0.65 | 2,782 / 2,324 |
| 99 | SL_Z6_R30_T100000_C5 | gas5_ics3 | 16 | -318 | 5,737 [4,901, 6,593] | 16.9 | 0.80 | no | 0.65 | 2,409 / 2,156 |
| 100 | SS_Z8_R20_T150_C2 | gas1_ics1 | 10 | -317 | 5,456 [4,950, 6,002] | 109.7 | 0.92 | no | 0.90 | 1,781 / 3,253 |
| 101 | SL_Z6_R20_T150_C5 | gas1_ics1 | 10 | -314 | 4,443 [3,908, 5,020] | 65.1 | 0.91 | no | 0.92 | 916 / 3,109 |
| 102 | SS_Z10_R20_T100000_C5 | gas5_ics3 | 16 | -312 | 5,737 [5,176, 6,264] | 13.9 | 0.72 | no | 0.65 | 2,298 / 1,818 |
| 103 | SS_Z6_R30_T150_C2 | gas1_ics1 | 10 | -310 | 5,568 [5,071, 6,121] | 115.3 | 0.94 | no | 0.89 | 1,916 / 3,292 |
| 104 | SS_Z10_R20_T150_C2 | gas1_ics1 | 10 | -303 | 5,399 [4,877, 5,930] | 106.3 | 0.91 | no | 0.90 | 1,707 / 3,226 |
| 105 | SY_Z6_R30_T150_C2 | gas1_ics1 | 10 | -301 | 4,893 [4,365, 5,447] | 82.7 | 0.93 | no | 0.91 | 1,408 / 3,147 |
| 106 | SS_Z6_R20_T150_C2 | gas1_ics1 | 10 | -300 | 5,568 [5,068, 6,117] | 115.2 | 0.94 | no | 0.89 | 1,916 / 3,292 |
| 107 | SL_Z6_R20_T100000_C5 | gas5_ics3 | 16 | -300 | 5,737 [4,905, 6,593] | 17.0 | 0.83 | no | 0.65 | 2,409 / 2,324 |
| 108 | SS_Z6_R20_T100000_C5 | gas5_ics3 | 16 | -299 | 5,849 [5,342, 6,407] | 22.8 | 0.90 | no | 0.65 | 2,901 / 2,351 |
| 109 | SS_Z6_R30_T100000_C2 | gas5_ics3 | 16 | -297 | 6,018 [5,390, 6,729] | 21.8 | 0.89 | no | 0.65 | 2,901 / 2,447 |
| 110 | SY_Z10_R20_T100000_C2 | gas5_ics3 | 16 | -296 | 5,681 [4,824, 6,587] | 11.8 | 0.66 | no | 0.66 | 1,969 / 1,781 |
| 111 | SS_Z10_R20_T100000_C2 | gas5_ics3 | 16 | -292 | 5,849 [5,215, 6,621] | 13.1 | 0.72 | no | 0.66 | 2,313 / 1,873 |
| 112 | SY_Z8_R30_T150_C2 | gas1_ics1 | 10 | -289 | 4,781 [4,234, 5,327] | 76.9 | 0.92 | no | 0.91 | 1,303 / 3,109 |
| 113 | SL_Z10_R20_T100000_C2 | gas5_ics3 | 16 | -284 | 5,793 [4,859, 6,687] | 8.7 | 0.58 | no | 0.66 | 1,776 / 1,611 |
| 114 | SL_Z10_R30_T150_C5 | gas1_ics1 | 10 | -284 | 4,275 [3,682, 4,821] | 56.3 | 0.83 | yes | 0.92 | 499 / 3,056 |
| 115 | SY_Z6_R20_T150_C2 | gas1_ics1 | 10 | -283 | 4,893 [4,365, 5,447] | 82.7 | 0.93 | no | 0.91 | 1,408 / 3,147 |
| 116 | SY_Z6_R30_T100000_C2 | gas5_ics3 | 16 | -276 | 5,906 [5,067, 6,727] | 18.4 | 0.84 | no | 0.65 | 2,618 / 2,328 |
| 117 | SS_Z6_R20_T100000_C2 | gas5_ics3 | 16 | -265 | 6,018 [5,391, 6,729] | 21.9 | 0.88 | no | 0.65 | 2,856 / 2,447 |
| 118 | SY_Z10_R30_T150_C2 | gas1_ics1 | 10 | -256 | 4,724 [4,153, 5,254] | 73.4 | 0.92 | no | 0.91 | 1,259 / 3,096 |
| 119 | SY_Z6_R20_T100000_C2 | gas5_ics3 | 16 | -256 | 5,906 [5,072, 6,728] | 18.4 | 0.84 | no | 0.65 | 2,588 / 2,369 |
| 120 | SL_Z8_R30_T150_C5 | gas1_ics1 | 10 | -254 | 4,331 [3,748, 4,876] | 58.9 | 0.83 | yes | 0.92 | 558 / 3,056 |
| 121 | SL_Z6_R30_T100000_C2 | gas5_ics3 | 16 | -238 | 5,962 [5,088, 6,764] | 16.4 | 0.78 | no | 0.65 | 2,409 / 2,224 |
| 122 | SO_Z6_R20_T150_C5 | gas1_ics1 | 10 | -236 | 5,962 [5,420, 6,465] | 132.3 | 0.87 | no | 0.88 | 1,812 / 3,397 |
| 123 | SL_Z6_R20_T100000_C2 | gas5_ics3 | 16 | -223 | 5,962 [5,091, 6,763] | 16.4 | 0.79 | no | 0.65 | 2,424 / 2,278 |
| 124 | SO_Z6_R20_T100000_C5 | gas5_ics3 | 16 | -218 | 6,243 [5,723, 6,762] | 22.3 | 0.79 | no | 0.65 | 2,616 / 2,337 |
| 125 | SY_Z8_R20_T150_C2 | gas1_ics1 | 10 | -213 | 4,781 [4,220, 5,313] | 76.5 | 0.92 | no | 0.91 | 1,288 / 3,109 |
| 126 | SL_Z10_R20_T150_C5 | gas1_ics1 | 10 | -210 | 4,275 [3,680, 4,817] | 56.4 | 0.81 | yes | 0.92 | 439 / 3,043 |
| 127 | SL_Z8_R20_T150_C5 | gas1_ics1 | 10 | -207 | 4,331 [3,745, 4,873] | 58.9 | 0.81 | yes | 0.92 | 454 / 3,056 |
| 128 | SO_Z6_R20_T100000_C2 | gas5_ics3 | 16 | -200 | 6,299 [5,769, 6,808] | 21.6 | 0.81 | no | 0.65 | 2,631 / 2,501 |
| 129 | SY_Z10_R20_T150_C2 | gas1_ics1 | 10 | -187 | 4,668 [4,139, 5,240] | 73.1 | 0.92 | no | 0.91 | 1,229 / 3,083 |
| 130 | SS_Z4_R30_T100000_C5 | gas5_ics3 | 16 | -115 | 5,962 [5,443, 6,632] | 40.0 | 0.93 | no | 0.64 | 2,738 / 2,815 |
| 131 | SY_Z4_R30_T150_C2 | gas1_ics1 | 10 | -111 | 5,062 [4,548, 5,619] | 91.0 | 0.94 | no | 0.91 | 1,587 / 3,186 |
| 132 | SY_Z4_R20_T150_C2 | gas1_ics1 | 10 | -106 | 5,062 [4,548, 5,619] | 91.0 | 0.94 | no | 0.91 | 1,587 / 3,186 |
| 133 | SS_Z4_R30_T150_C2 | gas1_ics1 | 10 | -105 | 5,737 [5,205, 6,252] | 121.4 | 0.94 | no | 0.89 | 2,050 / 3,317 |
| 134 | SL_Z4_R30_T150_C5 | gas1_ics1 | 10 | -101 | 4,668 [4,108, 5,200] | 73.3 | 0.88 | no | 0.91 | 931 / 3,160 |
| 135 | SS_Z4_R20_T150_C2 | gas1_ics1 | 10 | -100 | 5,737 [5,204, 6,251] | 121.3 | 0.94 | no | 0.89 | 2,050 / 3,317 |
| 136 | SY_Z4_R30_T100000_C2 | gas5_ics3 | 16 | -100 | 6,074 [5,312, 6,819] | 34.9 | 0.86 | no | 0.64 | 2,544 / 2,684 |
| 137 | SS_Z4_R20_T100000_C2 | gas5_ics3 | 16 | -98 | 6,131 [5,502, 6,768] | 39.1 | 0.90 | no | 0.64 | 2,693 / 2,814 |
| 138 | SL_Z4_R20_T150_C5 | gas1_ics1 | 10 | -96 | 4,668 [4,108, 5,200] | 73.3 | 0.88 | no | 0.91 | 931 / 3,160 |
| 139 | SY_Z4_R20_T100000_C2 | gas5_ics3 | 16 | -93 | 6,074 [5,312, 6,819] | 34.9 | 0.86 | no | 0.64 | 2,514 / 2,684 |
| 140 | SS_Z4_R30_T100000_C2 | gas5_ics3 | 16 | -83 | 6,131 [5,503, 6,769] | 39.0 | 0.91 | no | 0.64 | 2,738 / 2,814 |
| 141 | SL_Z4_R30_T100000_C2 | gas5_ics3 | 16 | -14 | 6,131 [5,338, 6,726] | 32.5 | 0.82 | no | 0.65 | 2,380 / 2,672 |
| 142 | SL_Z4_R20_T100000_C2 | gas5_ics3 | 16 | -12 | 6,131 [5,338, 6,726] | 32.5 | 0.82 | no | 0.65 | 2,380 / 2,672 |
| 143 | SL_Z6_R30_T150_C2 | gas1_ics1 | 10 | +0 | 4,556 [3,985, 5,098] | 67.0 | 0.89 | no | 0.92 | 976 / 3,095 |
| 144 | SO_Z4_R30_T150_C5 | gas1_ics1 | 10 | +8 | 6,018 [5,497, 6,543] | 135.4 | 0.89 | no | 0.88 | 1,931 / 3,422 |
| 145 | SO_Z4_R20_T150_C5 | gas1_ics1 | 10 | +10 | 6,018 [5,497, 6,542] | 135.4 | 0.89 | no | 0.88 | 1,931 / 3,422 |
| 146 | SY_Z8_R30_T100000_C5 | gas1_ics1 | 10 | +11 | 4,781 [4,254, 5,353] | 77.1 | 0.95 | no | 0.91 | 1,425 / 3,110 |
| 147 | SL_Z6_R20_T150_C2 | gas1_ics1 | 10 | +17 | 4,556 [3,986, 5,099] | 67.0 | 0.91 | no | 0.92 | 1,035 / 3,095 |
| 148 | SL_Z8_R30_T150_C2 | gas1_ics1 | 10 | +29 | 4,387 [3,846, 4,976] | 61.4 | 0.87 | no | 0.92 | 752 / 3,056 |
| 149 | SL_Z10_R30_T150_C2 | gas1_ics1 | 10 | +48 | 4,331 [3,767, 4,908] | 58.2 | 0.94 | no | 0.92 | 1,050 / 3,030 |
| 150 | SO_Z4_R30_T100000_C2 | gas5_ics3 | 16 | +62 | 6,356 [5,833, 6,872] | 39.7 | 0.89 | no | 0.64 | 2,767 / 2,868 |
| 151 | SO_Z4_R20_T100000_C2 | gas5_ics3 | 16 | +67 | 6,412 [5,846, 6,965] | 39.8 | 0.87 | no | 0.64 | 2,738 / 2,868 |
| 152 | SS_Z4_R20_T100000_C5 | gas5_ics3 | 16 | +79 | 5,962 [5,427, 6,478] | 40.2 | 0.94 | no | 0.64 | 2,783 / 2,814 |
| 153 | SY_Z4_R30_T100000_C5 | gas5_ics3 | 16 | +83 | 5,906 [5,144, 6,649] | 35.9 | 0.92 | no | 0.64 | 2,693 / 2,736 |
| 154 | SY_Z8_R30_T100000_C2 | gas1_ics1 | 10 | +88 | 4,837 [4,276, 5,376] | 77.3 | 0.92 | no | 0.91 | 1,379 / 3,084 |
| 155 | SY_Z4_R20_T100000_C5 | gas5_ics3 | 16 | +92 | 5,906 [5,144, 6,649] | 36.0 | 0.88 | no | 0.64 | 2,634 / 2,584 |
| 156 | SL_Z8_R30_T100000_C5 | gas1_ics1 | 10 | +92 | 4,443 [3,874, 5,010] | 61.7 | 0.95 | no | 0.92 | 1,186 / 3,045 |
| 157 | SL_Z8_R20_T150_C2 | gas1_ics1 | 10 | +107 | 4,387 [3,832, 4,960] | 61.1 | 0.88 | no | 0.92 | 826 / 3,056 |
| 158 | SL_Z10_R20_T150_C2 | gas1_ics1 | 10 | +120 | 4,331 [3,759, 4,898] | 58.1 | 0.85 | no | 0.92 | 647 / 3,030 |
| 159 | SO_Z6_R20_T150_C2 | gas1_ics1 | 10 | +129 | 6,018 [5,505, 6,547] | 136.4 | 0.89 | no | 0.88 | 1,991 / 3,383 |
| 160 | SY_Z8_R20_T100000_C2 | gas1_ics1 | 10 | +137 | 4,837 [4,273, 5,372] | 77.4 | 0.93 | no | 0.91 | 1,394 / 3,084 |
| 161 | SL_Z4_R30_T100000_C5 | gas5_ics3 | 16 | +172 | 5,906 [5,151, 6,609] | 33.4 | 0.87 | no | 0.65 | 2,455 / 2,698 |
| 162 | SY_Z10_R30_T100000_C5 | gas1_ics1 | 10 | +173 | 4,724 [4,165, 5,272] | 73.3 | 0.94 | no | 0.91 | 1,350 / 3,084 |
| 163 | SL_Z4_R20_T100000_C5 | gas5_ics3 | 16 | +177 | 5,906 [5,152, 6,609] | 33.4 | 0.87 | no | 0.65 | 2,455 / 2,698 |
| 164 | SY_Z8_R20_T100000_C5 | gas1_ics1 | 10 | +189 | 4,781 [4,234, 5,331] | 76.5 | 0.94 | no | 0.91 | 1,380 / 3,110 |
| 165 | SL_Z8_R30_T100000_C2 | gas1_ics1 | 10 | +201 | 4,443 [3,901, 5,041] | 62.0 | 0.92 | no | 0.92 | 1,051 / 3,032 |
| 166 | SY_Z6_R30_T100000_C5 | gas1_ics1 | 10 | +201 | 4,893 [4,357, 5,444] | 81.9 | 0.94 | no | 0.91 | 1,484 / 3,136 |
| 167 | SS_Z8_R30_T100000_C5 | gas1_ics1 | 10 | +204 | 5,456 [4,954, 6,010] | 109.5 | 0.94 | no | 0.90 | 1,873 / 3,253 |
| 168 | SY_Z10_R30_T100000_C2 | gas1_ics1 | 10 | +218 | 4,724 [4,186, 5,295] | 73.4 | 0.91 | no | 0.91 | 1,245 / 3,058 |
| 169 | SL_Z4_R30_T150_C2 | gas1_ics1 | 10 | +224 | 4,724 [4,182, 5,275] | 75.3 | 0.91 | no | 0.91 | 1,140 / 3,147 |
| 170 | SY_Z6_R20_T100000_C5 | gas1_ics1 | 10 | +226 | 4,893 [4,356, 5,443] | 81.9 | 0.94 | no | 0.91 | 1,484 / 3,136 |
| 171 | SL_Z4_R20_T150_C2 | gas1_ics1 | 10 | +229 | 4,724 [4,182, 5,274] | 75.3 | 0.90 | no | 0.91 | 1,095 / 3,147 |
| 172 | SS_Z8_R30_T100000_C2 | gas1_ics1 | 10 | +232 | 5,512 [4,982, 6,036] | 110.3 | 0.94 | no | 0.90 | 1,932 / 3,240 |
| 173 | SL_Z8_R20_T100000_C2 | gas1_ics1 | 10 | +248 | 4,443 [3,901, 5,040] | 62.2 | 1.01 | no | 0.92 | 1,468 / 3,032 |
| 174 | SY_Z6_R30_T100000_C2 | gas1_ics1 | 10 | +262 | 4,949 [4,380, 5,470] | 82.3 | 0.97 | no | 0.91 | 1,678 / 3,136 |
| 175 | SS_Z8_R20_T100000_C2 | gas1_ics1 | 10 | +267 | 5,512 [4,972, 6,027] | 109.9 | 0.94 | no | 0.90 | 1,932 / 3,226 |
| 176 | SL_Z10_R30_T100000_C5 | gas1_ics1 | 10 | +269 | 4,331 [3,785, 4,933] | 58.3 | 0.95 | no | 0.92 | 1,111 / 3,019 |
| 177 | SO_Z4_R30_T100000_C5 | gas5_ics3 | 16 | +274 | 6,243 [5,756, 6,853] | 40.8 | 0.91 | no | 0.64 | 2,828 / 2,856 |
| 178 | SY_Z10_R20_T100000_C5 | gas1_ics1 | 10 | +278 | 4,724 [4,153, 5,259] | 73.0 | 0.93 | no | 0.91 | 1,320 / 3,084 |
| 179 | SO_Z4_R20_T100000_C5 | gas5_ics3 | 16 | +281 | 6,299 [5,756, 6,798] | 40.8 | 0.89 | no | 0.64 | 2,783 / 2,843 |
| 180 | SY_Z6_R20_T100000_C2 | gas1_ics1 | 10 | +287 | 4,949 [4,379, 5,468] | 82.3 | 0.96 | no | 0.91 | 1,633 / 3,122 |
| 181 | SL_Z6_R30_T100000_C5 | gas1_ics1 | 10 | +296 | 4,556 [3,977, 5,097] | 66.2 | 0.86 | no | 0.91 | 829 / 3,097 |
| 182 | SS_Z8_R20_T100000_C5 | gas1_ics1 | 10 | +303 | 5,456 [4,937, 5,991] | 108.7 | 0.94 | no | 0.90 | 1,873 / 3,240 |
| 183 | SS_Z10_R30_T100000_C2 | gas1_ics1 | 10 | +313 | 5,456 [4,903, 5,961] | 106.5 | 0.93 | no | 0.90 | 1,872 / 3,214 |
| 184 | SS_Z10_R30_T100000_C5 | gas1_ics1 | 10 | +318 | 5,399 [4,878, 5,935] | 105.8 | 0.93 | no | 0.90 | 1,813 / 3,227 |
| 185 | SL_Z8_R20_T100000_C5 | gas1_ics1 | 10 | +318 | 4,443 [3,850, 4,985] | 61.1 | 0.95 | no | 0.92 | 1,156 / 3,045 |
| 186 | SL_Z6_R20_T100000_C5 | gas1_ics1 | 10 | +320 | 4,556 [3,976, 5,096] | 66.2 | 1.00 | no | 0.91 | 1,469 / 3,084 |
| 187 | SY_Z10_R20_T100000_C2 | gas1_ics1 | 10 | +338 | 4,724 [4,172, 5,279] | 73.2 | 0.92 | no | 0.91 | 1,305 / 3,058 |
| 188 | SL_Z10_R30_T100000_C2 | gas1_ics1 | 10 | +352 | 4,387 [3,811, 4,961] | 58.5 | 1.00 | no | 0.92 | 1,393 / 3,005 |
| 189 | SO_Z4_R30_T150_C2 | gas1_ics1 | 10 | +368 | 6,131 [5,585, 6,630] | 139.6 | 0.90 | no | 0.88 | 2,110 / 3,409 |
| 190 | SO_Z4_R20_T150_C2 | gas1_ics1 | 10 | +369 | 6,131 [5,584, 6,629] | 139.6 | 0.90 | no | 0.88 | 2,110 / 3,409 |
| 191 | SL_Z10_R20_T100000_C5 | gas1_ics1 | 10 | +373 | 4,331 [3,778, 4,923] | 58.3 | 0.94 | no | 0.92 | 1,066 / 3,019 |
| 192 | SS_Z6_R30_T100000_C5 | gas1_ics1 | 10 | +394 | 5,568 [5,043, 6,094] | 113.6 | 0.94 | no | 0.89 | 1,962 / 3,265 |
| 193 | SS_Z10_R20_T100000_C5 | gas1_ics1 | 10 | +395 | 5,399 [4,861, 5,918] | 105.1 | 0.93 | no | 0.90 | 1,798 / 3,214 |
| 194 | SL_Z6_R30_T100000_C2 | gas1_ics1 | 10 | +399 | 4,556 [4,004, 5,126] | 66.5 | 1.00 | no | 0.91 | 1,498 / 3,070 |
| 195 | SS_Z6_R30_T100000_C2 | gas1_ics1 | 10 | +406 | 5,624 [5,072, 6,125] | 114.7 | 0.94 | no | 0.89 | 2,021 / 3,265 |
| 196 | SS_Z10_R20_T100000_C2 | gas1_ics1 | 10 | +407 | 5,399 [4,883, 5,941] | 105.7 | 0.94 | no | 0.90 | 1,857 / 3,201 |
| 197 | SY_Z4_R30_T100000_C2 | gas1_ics1 | 10 | +409 | 5,118 [4,602, 5,675] | 92.4 | 0.97 | no | 0.90 | 1,753 / 3,187 |
| 198 | SS_Z6_R20_T100000_C5 | gas1_ics1 | 10 | +413 | 5,568 [5,038, 6,090] | 113.4 | 0.94 | no | 0.89 | 1,962 / 3,265 |
| 199 | SY_Z4_R20_T100000_C2 | gas1_ics1 | 10 | +417 | 5,118 [4,601, 5,674] | 92.4 | 0.97 | no | 0.90 | 1,753 / 3,187 |
| 200 | SL_Z6_R20_T100000_C2 | gas1_ics1 | 10 | +421 | 4,556 [4,004, 5,125] | 66.6 | 0.97 | no | 0.91 | 1,334 / 3,070 |
| 201 | SS_Z6_R20_T100000_C2 | gas1_ics1 | 10 | +423 | 5,568 [5,068, 6,120] | 114.4 | 0.95 | no | 0.89 | 2,021 / 3,265 |
| 202 | SL_Z10_R20_T100000_C2 | gas1_ics1 | 10 | +466 | 4,387 [3,802, 4,950] | 58.5 | 1.00 | no | 0.92 | 1,393 / 2,993 |
| 203 | SS_Z4_R30_T100000_C5 | gas1_ics1 | 10 | +554 | 5,737 [5,204, 6,253] | 121.2 | 0.95 | no | 0.89 | 2,141 / 3,317 |
| 204 | SS_Z4_R30_T100000_C2 | gas1_ics1 | 10 | +580 | 5,793 [5,243, 6,292] | 122.7 | 0.95 | no | 0.89 | 2,200 / 3,303 |
| 205 | SL_Z4_R30_T100000_C2 | gas1_ics1 | 10 | +586 | 4,781 [4,239, 5,336] | 76.5 | 1.03 | no | 0.91 | 1,782 / 3,135 |
| 206 | SS_Z4_R20_T100000_C2 | gas1_ics1 | 10 | +588 | 5,793 [5,242, 6,291] | 122.7 | 0.95 | no | 0.89 | 2,200 / 3,303 |
| 207 | SL_Z4_R20_T100000_C2 | gas1_ics1 | 10 | +589 | 4,781 [4,239, 5,336] | 76.5 | 0.97 | no | 0.91 | 1,499 / 3,135 |
| 208 | SY_Z4_R30_T100000_C5 | gas1_ics1 | 10 | +620 | 5,062 [4,543, 5,616] | 90.4 | 0.96 | no | 0.90 | 1,649 / 3,187 |
| 209 | SY_Z4_R20_T100000_C5 | gas1_ics1 | 10 | +630 | 5,062 [4,542, 5,616] | 90.4 | 0.96 | no | 0.90 | 1,664 / 3,187 |
| 210 | SS_Z4_R20_T100000_C5 | gas1_ics1 | 10 | +745 | 5,681 [5,181, 6,232] | 120.0 | 0.96 | no | 0.89 | 2,126 / 3,317 |
| 211 | SO_Z6_R20_T100000_C5 | gas1_ics1 | 10 | +764 | 5,962 [5,448, 6,494] | 133.2 | 0.91 | no | 0.88 | 2,067 / 3,356 |
| 212 | SO_Z6_R20_T100000_C2 | gas1_ics1 | 10 | +777 | 6,018 [5,491, 6,535] | 135.1 | 0.91 | no | 0.88 | 2,126 / 3,343 |
| 213 | SL_Z4_R30_T100000_C5 | gas1_ics1 | 10 | +789 | 4,724 [4,176, 5,273] | 74.7 | 0.91 | no | 0.91 | 1,186 / 3,135 |
| 214 | SL_Z4_R20_T100000_C5 | gas1_ics1 | 10 | +796 | 4,724 [4,176, 5,273] | 74.7 | 0.93 | no | 0.91 | 1,276 / 3,135 |
| 215 | SO_Z4_R30_T100000_C2 | gas1_ics1 | 10 | +1,013 | 6,131 [5,610, 6,655] | 140.4 | 0.92 | no | 0.88 | 2,275 / 3,381 |
| 216 | SO_Z4_R20_T100000_C2 | gas1_ics1 | 10 | +1,016 | 6,131 [5,609, 6,654] | 140.4 | 0.92 | no | 0.88 | 2,275 / 3,381 |
| 217 | SO_Z4_R30_T100000_C5 | gas1_ics1 | 10 | +1,180 | 6,074 [5,537, 6,584] | 136.7 | 0.92 | no | 0.88 | 2,186 / 3,381 |
| 218 | SO_Z4_R20_T100000_C5 | gas1_ics1 | 10 | +1,185 | 6,074 [5,536, 6,584] | 136.7 | 0.92 | no | 0.88 | 2,186 / 3,381 |
