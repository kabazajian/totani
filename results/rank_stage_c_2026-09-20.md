# Stage C — 2026-09-20

Consensus selection (median Δ(−2lnL) across 217 model rows; BIC threshold p·ln 1080 = 63 per template): **bubble_totani → loopI_shells → sun_moon_ic**

| step | chosen | median Δ | min | max | threshold | runner-up (median Δ) |
|---|---|---|---|---|---|---|
| 1 | bubble_totani | 28,819 | 20,615 | 39,682 | 126 | bubble_flat (19,194) |
| 2 | loopI_shells | 7,042 | 4,891 | 8,492 | 126 | loopI_shells+nps (7,035) |
| 3 | sun_moon_ic | 81 | 12 | 146 | 63 | erosita_outer_NS (71) |
| 4 | (stop) erosita_outer_NS | 57 | -95 | 129 | 126 | loopII (42) |

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

## Ranking with the consensus set (bubble_totani, loopI_shells, sun_moon_ic)

| rank | model | variant | n | Δ(−2lnL) vs Totani baseline | halo @ 20.8 | TS @ 20.8 | quad Σ/joint | any neg | inj min | N / S |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | a2020 | - | 21 | -1,573 | 8,549 [8,056, 9,054] | 43.2 | 0.71 | no | 0.98 | 3,268 / 2,798 |
| 2 | SS_Z10_R30_T150_C2 | gas5_ics3 | 15 | -1,466 | 6,074 [5,175, 6,988] | 17.9 | 0.14 | yes | 0.97 | -4 / 824 |
| 3 | SS_Z10_R20_T150_C2 | gas5_ics3 | 15 | -1,462 | 6,131 [5,227, 7,037] | 18.3 | 0.19 | yes | 0.99 | 71 / 1,096 |
| 4 | SY_Z10_R30_T150_C2 | gas5_ics3 | 15 | -1,454 | 6,187 [5,252, 7,132] | 18.6 | 0.32 | yes | 1.00 | 522 / 1,472 |
| 5 | SS_Z8_R30_T150_C2 | gas5_ics3 | 15 | -1,453 | 6,187 [5,271, 7,078] | 20.3 | 0.34 | yes | 0.99 | 673 / 1,403 |
| 6 | SY_Z10_R20_T150_C2 | gas5_ics3 | 15 | -1,451 | 6,187 [5,268, 7,145] | 18.6 | 0.31 | yes | 0.99 | 597 / 1,313 |
| 7 | SS_Z8_R20_T150_C2 | gas5_ics3 | 15 | -1,450 | 6,243 [5,326, 7,129] | 20.6 | 0.36 | yes | 0.99 | 763 / 1,469 |
| 8 | SY_Z8_R30_T150_C2 | gas5_ics3 | 15 | -1,436 | 6,299 [5,358, 7,222] | 20.9 | 0.48 | yes | 0.98 | 1,168 / 1,845 |
| 9 | SS_Z10_R30_T150_C5 | gas5_ics3 | 15 | -1,433 | 6,412 [5,536, 7,412] | 27.7 | 0.50 | yes | 0.99 | 923 / 2,291 |
| 10 | SY_Z8_R20_T150_C2 | gas5_ics3 | 15 | -1,432 | 6,299 [5,396, 7,256] | 21.0 | 0.46 | yes | 0.98 | 1,243 / 1,683 |
| 11 | SS_Z10_R20_T150_C5 | gas5_ics3 | 15 | -1,430 | 6,749 [5,538, 7,466] | 23.4 | 0.49 | yes | 0.97 | 998 / 2,335 |
| 12 | SS_Z8_R30_T150_C5 | gas5_ics3 | 15 | -1,429 | 6,580 [5,653, 7,494] | 25.6 | 0.65 | yes | 0.99 | 1,764 / 2,539 |
| 13 | SL_Z10_R30_T150_C2 | gas5_ics3 | 15 | -1,428 | 6,299 [5,340, 7,237] | 19.6 | 0.44 | yes | 0.98 | 866 / 1,894 |
| 14 | SS_Z8_R20_T150_C5 | gas5_ics3 | 15 | -1,425 | 6,580 [5,672, 7,506] | 25.8 | 0.63 | yes | 0.98 | 1,809 / 2,313 |
| 15 | SS_Z6_R30_T150_C2 | gas5_ics3 | 15 | -1,424 | 6,412 [5,515, 7,299] | 26.0 | 0.58 | yes | 0.98 | 1,860 / 1,847 |
| 16 | SS_Z6_R20_T150_C2 | gas5_ics3 | 15 | -1,423 | 6,412 [5,523, 7,307] | 26.1 | 0.60 | yes | 0.98 | 1,875 / 1,987 |
| 17 | SL_Z10_R20_T150_C2 | gas5_ics3 | 15 | -1,421 | 6,299 [5,357, 7,252] | 19.7 | 0.46 | yes | 0.98 | 911 / 1,974 |
| 18 | SY_Z10_R30_T150_C5 | gas5_ics3 | 15 | -1,410 | 6,412 [5,478, 7,334] | 22.7 | 0.54 | yes | 0.99 | 1,373 / 2,089 |
| 19 | SL_Z8_R30_T150_C2 | gas5_ics3 | 15 | -1,408 | 6,412 [5,466, 7,344] | 22.1 | 0.54 | yes | 0.99 | 1,422 / 2,025 |
| 20 | SY_Z10_R20_T150_C5 | gas5_ics3 | 15 | -1,405 | 6,412 [5,505, 7,359] | 22.8 | 0.57 | yes | 0.67 | 1,492 / 2,157 |
| 21 | SS_Z6_R30_T150_C5 | gas5_ics3 | 15 | -1,405 | 6,637 [5,774, 7,497] | 30.0 | 0.85 | no | 0.97 | 2,965 / 2,658 |
| 22 | SS_Z6_R20_T150_C5 | gas5_ics3 | 15 | -1,404 | 6,637 [5,780, 7,506] | 30.0 | 0.86 | no | 0.99 | 2,965 / 2,723 |
| 23 | SL_Z8_R20_T150_C2 | gas5_ics3 | 15 | -1,403 | 6,412 [5,499, 7,375] | 22.1 | 0.54 | yes | 1.00 | 1,466 / 2,015 |
| 24 | SL_Z10_R30_T150_C5 | gas5_ics3 | 15 | -1,401 | 6,580 [5,624, 7,493] | 22.3 | 0.79 | no | 0.99 | 2,879 / 2,334 |
| 25 | SY_Z8_R30_T150_C5 | gas5_ics3 | 15 | -1,401 | 6,524 [5,581, 7,412] | 21.8 | 0.71 | yes | 0.98 | 2,228 / 2,372 |
| 26 | SY_Z8_R20_T150_C5 | gas5_ics3 | 15 | -1,396 | 6,524 [5,595, 7,425] | 21.8 | 0.68 | yes | 0.98 | 2,198 / 2,237 |
| 27 | SL_Z10_R20_T150_C5 | gas5_ics3 | 15 | -1,391 | 6,524 [5,645, 7,514] | 22.7 | 0.61 | yes | 0.98 | 1,703 / 2,308 |
| 28 | SL_Z8_R30_T150_C5 | gas5_ics3 | 15 | -1,376 | 6,862 [5,718, 7,598] | 24.7 | 0.64 | yes | 0.66 | 2,259 / 2,117 |
| 29 | SL_Z8_R20_T150_C5 | gas5_ics3 | 15 | -1,370 | 6,918 [5,718, 7,617] | 24.8 | 0.64 | yes | 0.98 | 2,258 / 2,184 |
| 30 | SO_Z6_R20_T150_C5 | gas5_ics3 | 15 | -1,369 | 6,299 [5,626, 7,344] | 35.4 | 0.78 | yes | 0.98 | 2,379 / 2,553 |
| 31 | SY_Z10_R30_T100000_C2 | gas5_ics3 | 15 | -1,369 | 6,918 [5,988, 7,835] | 18.9 | 0.51 | yes | 0.98 | 1,644 / 1,915 |
| 32 | SY_Z6_R20_T150_C5 | gas5_ics3 | 15 | -1,368 | 6,637 [5,708, 7,472] | 26.2 | 0.80 | no | 0.97 | 2,920 / 2,418 |
| 33 | SL_Z6_R30_T150_C2 | gas5_ics3 | 15 | -1,367 | 6,580 [5,661, 7,488] | 27.5 | 0.71 | yes | 0.96 | 2,368 / 2,288 |
| 34 | SO_Z6_R20_T150_C2 | gas5_ics3 | 15 | -1,367 | 6,299 [5,587, 7,367] | 35.8 | 0.57 | yes | 0.97 | 1,573 / 2,020 |
| 35 | SY_Z6_R30_T150_C5 | gas5_ics3 | 15 | -1,366 | 6,580 [5,704, 7,466] | 26.2 | 0.82 | no | 0.99 | 2,905 / 2,459 |
| 36 | SL_Z6_R20_T150_C2 | gas5_ics3 | 15 | -1,363 | 6,580 [5,667, 7,493] | 27.5 | 0.71 | yes | 1.00 | 2,368 / 2,301 |
| 37 | SS_Z4_R20_T150_C2 | gas5_ics3 | 15 | -1,361 | 6,862 [6,012, 7,538] | 41.2 | 0.83 | no | 0.99 | 2,968 / 2,754 |
| 38 | SS_Z4_R30_T150_C2 | gas5_ics3 | 15 | -1,361 | 6,862 [6,009, 7,534] | 41.2 | 0.84 | no | 1.00 | 2,983 / 2,754 |
| 39 | SY_Z8_R30_T100000_C2 | gas5_ics3 | 15 | -1,360 | 6,862 [5,959, 7,792] | 20.6 | 0.60 | yes | 0.98 | 2,169 / 1,933 |
| 40 | SY_Z10_R20_T100000_C2 | gas5_ics3 | 15 | -1,355 | 6,974 [6,045, 7,887] | 19.2 | 0.54 | yes | 0.98 | 1,793 / 1,967 |
| 41 | SY_Z6_R30_T150_C2 | gas5_ics3 | 15 | -1,354 | 6,637 [5,700, 7,511] | 27.4 | 0.67 | yes | 0.99 | 2,188 / 2,244 |
| 42 | SY_Z8_R20_T100000_C2 | gas5_ics3 | 15 | -1,353 | 6,918 [5,983, 7,813] | 20.7 | 0.58 | yes | 0.97 | 2,154 / 1,870 |
| 43 | SY_Z6_R20_T150_C2 | gas5_ics3 | 15 | -1,353 | 6,637 [5,704, 7,516] | 27.4 | 0.69 | yes | 1.00 | 2,218 / 2,333 |
| 44 | SS_Z8_R30_T100000_C2 | gas5_ics3 | 15 | -1,349 | 6,862 [5,927, 7,799] | 22.8 | 0.53 | yes | 0.98 | 1,690 / 1,962 |
| 45 | SS_Z10_R30_T100000_C2 | gas5_ics3 | 15 | -1,346 | 6,862 [5,933, 7,834] | 21.0 | 0.48 | yes | 0.98 | 1,133 / 2,144 |
| 46 | SS_Z8_R20_T100000_C2 | gas5_ics3 | 15 | -1,342 | 7,199 [6,170, 8,040] | 27.3 | 0.60 | yes | 1.00 | 1,900 / 2,455 |
| 47 | SS_Z10_R20_T100000_C2 | gas5_ics3 | 15 | -1,336 | 6,974 [5,999, 7,891] | 21.6 | 0.51 | yes | 0.98 | 1,297 / 2,253 |
| 48 | SL_Z10_R30_T100000_C2 | gas5_ics3 | 15 | -1,334 | 7,087 [6,138, 8,000] | 19.5 | 0.63 | yes | 0.98 | 2,272 / 2,173 |
| 49 | SS_Z4_R20_T150_C5 | gas5_ics3 | 15 | -1,329 | 6,974 [6,110, 7,588] | 43.0 | 0.91 | no | 0.99 | 3,476 / 2,873 |
| 50 | SS_Z4_R30_T150_C5 | gas5_ics3 | 15 | -1,328 | 6,918 [6,107, 7,585] | 43.0 | 0.92 | no | 0.99 | 3,402 / 2,975 |
| 51 | SL_Z6_R30_T150_C5 | gas5_ics3 | 15 | -1,328 | 6,974 [5,850, 7,651] | 30.2 | 0.75 | no | 0.66 | 2,905 / 2,311 |
| 52 | SS_Z6_R30_T100000_C2 | gas5_ics3 | 15 | -1,327 | 6,974 [6,091, 7,823] | 27.5 | 0.71 | yes | 0.99 | 2,784 / 2,201 |
| 53 | SL_Z6_R20_T150_C5 | gas5_ics3 | 15 | -1,326 | 6,749 [5,879, 7,641] | 30.3 | 0.77 | no | 1.00 | 2,905 / 2,269 |
| 54 | SS_Z6_R20_T100000_C2 | gas5_ics3 | 15 | -1,324 | 6,974 [6,104, 7,836] | 27.6 | 0.75 | yes | 0.99 | 2,814 / 2,417 |
| 55 | SL_Z10_R20_T100000_C2 | gas5_ics3 | 15 | -1,321 | 7,087 [6,176, 8,034] | 19.6 | 0.55 | yes | 0.67 | 2,049 / 1,849 |
| 56 | SY_Z6_R30_T100000_C2 | gas5_ics3 | 15 | -1,319 | 6,974 [6,090, 7,878] | 24.3 | 0.70 | yes | 1.00 | 2,783 / 2,120 |
| 57 | SL_Z8_R30_T100000_C2 | gas5_ics3 | 15 | -1,316 | 7,087 [6,139, 7,986] | 21.8 | 0.57 | yes | 0.66 | 2,051 / 1,981 |
| 58 | SY_Z6_R20_T100000_C2 | gas5_ics3 | 15 | -1,316 | 6,974 [6,099, 7,886] | 24.3 | 0.68 | yes | 0.98 | 2,783 / 1,968 |
| 59 | SY_Z4_R20_T150_C5 | gas5_ics3 | 15 | -1,314 | 6,862 [6,024, 7,537] | 37.2 | 0.84 | no | 0.98 | 3,281 / 2,469 |
| 60 | SY_Z4_R30_T150_C5 | gas5_ics3 | 15 | -1,313 | 6,862 [6,022, 7,535] | 37.2 | 0.86 | no | 0.98 | 3,432 / 2,469 |
| 61 | SO_Z4_R20_T150_C5 | gas5_ics3 | 15 | -1,312 | 6,693 [6,138, 7,211] | 45.4 | 0.89 | no | 0.99 | 3,297 / 2,666 |
| 62 | SL_Z8_R20_T100000_C2 | gas5_ics3 | 15 | -1,309 | 7,087 [6,154, 7,998] | 21.9 | 0.61 | yes | 0.97 | 2,394 / 1,943 |
| 63 | SO_Z4_R30_T150_C5 | gas5_ics3 | 15 | -1,309 | 6,693 [6,135, 7,208] | 45.4 | 0.88 | no | 1.00 | 3,297 / 2,589 |
| 64 | SY_Z4_R20_T150_C2 | gas5_ics3 | 15 | -1,304 | 6,918 [6,060, 7,590] | 39.9 | 0.79 | no | 0.98 | 2,908 / 2,538 |
| 65 | SY_Z4_R30_T150_C2 | gas5_ics3 | 15 | -1,302 | 6,805 [6,051, 7,593] | 39.9 | 0.82 | no | 0.98 | 3,028 / 2,551 |
| 66 | SY_Z10_R30_T100000_C5 | gas5_ics3 | 15 | -1,298 | 7,199 [6,275, 8,097] | 25.6 | 0.61 | yes | 0.98 | 2,333 / 2,077 |
| 67 | SY_Z8_R30_T100000_C5 | gas5_ics3 | 15 | -1,292 | 7,143 [6,222, 8,029] | 26.8 | 0.78 | no | 0.66 | 3,011 / 2,564 |
| 68 | SO_Z4_R30_T150_C2 | gas5_ics3 | 15 | -1,289 | 6,749 [6,207, 7,279] | 48.0 | 0.84 | no | 0.99 | 3,027 / 2,627 |
| 69 | SO_Z4_R20_T150_C2 | gas5_ics3 | 15 | -1,288 | 6,749 [6,205, 7,277] | 48.0 | 0.84 | no | 0.99 | 3,057 / 2,627 |
| 70 | SS_Z8_R30_T100000_C5 | gas5_ics3 | 15 | -1,285 | 7,368 [6,186, 8,090] | 29.1 | 0.68 | yes | 0.98 | 2,304 / 2,701 |
| 71 | SS_Z10_R30_T100000_C5 | gas5_ics3 | 15 | -1,280 | 7,199 [6,272, 8,135] | 27.7 | 0.58 | yes | 0.99 | 1,510 / 2,651 |
| 72 | SL_Z10_R30_T100000_C5 | gas5_ics3 | 15 | -1,279 | 7,143 [6,246, 8,086] | 22.9 | 0.71 | no | 0.66 | 2,966 / 2,079 |
| 73 | SY_Z10_R20_T100000_C5 | gas5_ics3 | 15 | -1,275 | 7,030 [6,129, 7,949] | 23.7 | 0.61 | yes | 0.66 | 2,153 / 2,110 |
| 74 | SL_Z4_R30_T150_C2 | gas5_ics3 | 15 | -1,275 | 6,862 [6,089, 7,615] | 41.0 | 0.81 | no | 0.99 | 3,118 / 2,416 |
| 75 | SL_Z4_R20_T150_C2 | gas5_ics3 | 15 | -1,274 | 6,918 [6,098, 7,617] | 41.0 | 0.80 | no | 0.99 | 3,088 / 2,416 |
| 76 | SS_Z8_R20_T100000_C5 | gas5_ics3 | 15 | -1,274 | 7,255 [6,268, 8,111] | 29.6 | 0.69 | yes | 0.99 | 2,334 / 2,688 |
| 77 | SS_Z10_R20_T100000_C5 | gas5_ics3 | 15 | -1,272 | 7,255 [6,317, 8,173] | 28.3 | 0.60 | yes | 0.99 | 1,629 / 2,742 |
| 78 | SY_Z8_R20_T100000_C5 | gas5_ics3 | 15 | -1,266 | 7,030 [6,144, 7,938] | 23.4 | 0.70 | yes | 1.00 | 2,439 / 2,500 |
| 79 | SL_Z10_R20_T100000_C5 | gas5_ics3 | 15 | -1,266 | 7,199 [6,266, 8,101] | 22.9 | 0.61 | yes | 0.97 | 2,438 / 1,951 |
| 80 | SL_Z8_R30_T100000_C5 | gas5_ics3 | 15 | -1,266 | 7,143 [6,211, 8,033] | 24.6 | 0.65 | yes | 0.99 | 2,649 / 1,972 |
| 81 | SS_Z4_R20_T100000_C2 | gas5_ics3 | 15 | -1,261 | 7,143 [6,310, 7,873] | 40.0 | 0.76 | yes | 0.99 | 3,008 / 2,406 |
| 82 | SS_Z6_R30_T100000_C5 | gas5_ics3 | 15 | -1,260 | 7,199 [6,307, 8,069] | 33.0 | 0.79 | no | 1.00 | 2,995 / 2,724 |
| 83 | SS_Z4_R30_T100000_C2 | gas5_ics3 | 15 | -1,260 | 7,143 [6,306, 7,867] | 40.0 | 0.83 | no | 1.00 | 3,279 / 2,634 |
| 84 | SS_Z6_R20_T100000_C5 | gas5_ics3 | 15 | -1,259 | 7,199 [6,318, 8,078] | 33.1 | 0.80 | no | 0.94 | 3,010 / 2,749 |
| 85 | SL_Z6_R30_T100000_C2 | gas5_ics3 | 15 | -1,257 | 7,199 [6,328, 8,114] | 27.9 | 0.78 | no | 1.00 | 3,234 / 2,348 |
| 86 | SL_Z6_R20_T100000_C2 | gas5_ics3 | 15 | -1,257 | 7,199 [6,310, 8,102] | 26.4 | 0.67 | yes | 1.00 | 2,828 / 2,018 |
| 87 | SO_Z6_R20_T100000_C2 | gas5_ics3 | 15 | -1,255 | 6,918 [6,073, 7,721] | 34.2 | 0.68 | yes | 0.99 | 2,275 / 2,435 |
| 88 | SO_Z6_R20_T100000_C5 | gas5_ics3 | 15 | -1,253 | 6,862 [6,020, 7,617] | 33.4 | 0.75 | yes | 1.00 | 2,425 / 2,706 |
| 89 | SL_Z8_R20_T100000_C5 | gas5_ics3 | 15 | -1,247 | 7,199 [6,309, 8,116] | 25.0 | 0.65 | yes | 0.98 | 2,664 / 1,987 |
| 90 | SY_Z6_R20_T100000_C5 | gas5_ics3 | 15 | -1,247 | 7,312 [6,352, 8,088] | 30.9 | 0.72 | yes | 0.66 | 2,934 / 2,323 |
| 91 | SO_Z4_R20_T100000_C2 | gas5_ics3 | 15 | -1,245 | 7,030 [6,530, 7,576] | 43.3 | 0.83 | no | 0.99 | 3,220 / 2,587 |
| 92 | SY_Z6_R30_T100000_C5 | gas5_ics3 | 15 | -1,240 | 7,087 [6,226, 7,947] | 27.1 | 0.68 | yes | 0.97 | 2,769 / 2,064 |
| 93 | SL_Z4_R20_T150_C5 | gas5_ics3 | 15 | -1,239 | 6,974 [6,190, 7,668] | 40.1 | 0.82 | no | 1.00 | 3,341 / 2,397 |
| 94 | SL_Z4_R30_T150_C5 | gas5_ics3 | 15 | -1,233 | 6,637 [6,092, 7,666] | 37.8 | 0.88 | no | 0.98 | 3,461 / 2,396 |
| 95 | SY_Z4_R20_T100000_C2 | gas5_ics3 | 15 | -1,218 | 7,143 [6,281, 7,878] | 35.7 | 0.74 | yes | 1.00 | 3,053 / 2,220 |
| 96 | SY_Z4_R30_T100000_C2 | gas5_ics3 | 15 | -1,218 | 7,143 [6,279, 7,876] | 35.8 | 0.74 | yes | 1.00 | 2,918 / 2,334 |
| 97 | SO_Z4_R30_T100000_C2 | gas5_ics3 | 15 | -1,209 | 7,030 [6,491, 7,538] | 43.8 | 0.81 | no | 1.00 | 3,205 / 2,522 |
| 98 | SL_Z6_R30_T100000_C5 | gas5_ics3 | 15 | -1,207 | 7,255 [6,364, 8,097] | 29.0 | 0.75 | no | 0.66 | 3,175 / 2,284 |
| 99 | SL_Z6_R20_T100000_C5 | gas5_ics3 | 15 | -1,206 | 7,255 [6,364, 8,097] | 28.9 | 0.65 | yes | 0.96 | 2,754 / 1,993 |
| 100 | SS_Z4_R30_T100000_C5 | gas5_ics3 | 15 | -1,196 | 7,143 [6,343, 7,871] | 41.3 | 0.87 | no | 0.97 | 3,399 / 2,816 |
| 101 | SS_Z4_R20_T100000_C5 | gas5_ics3 | 15 | -1,195 | 7,255 [6,405, 7,920] | 41.2 | 0.84 | no | 0.99 | 3,385 / 2,679 |
| 102 | SO_Z4_R20_T100000_C5 | gas5_ics3 | 15 | -1,194 | 6,974 [6,444, 7,488] | 41.6 | 0.82 | no | 0.99 | 3,207 / 2,524 |
| 103 | SO_Z4_R30_T100000_C5 | gas5_ics3 | 15 | -1,193 | 6,974 [6,443, 7,487] | 41.6 | 0.81 | no | 0.99 | 3,237 / 2,410 |
| 104 | SL_Z4_R30_T100000_C2 | gas5_ics3 | 15 | -1,163 | 7,312 [6,480, 8,084] | 38.4 | 0.73 | no | 0.99 | 3,159 / 2,200 |
| 105 | SL_Z4_R20_T100000_C2 | gas5_ics3 | 15 | -1,163 | 7,312 [6,482, 8,085] | 38.4 | 0.73 | no | 0.99 | 3,144 / 2,200 |
| 106 | SY_Z4_R20_T100000_C5 | gas5_ics3 | 15 | -1,154 | 7,312 [6,451, 7,988] | 37.5 | 0.76 | no | 0.94 | 3,264 / 2,289 |
| 107 | SY_Z4_R30_T100000_C5 | gas5_ics3 | 15 | -1,153 | 7,199 [6,444, 7,987] | 37.5 | 0.77 | no | 1.00 | 3,264 / 2,289 |
| 108 | SL_Z4_R20_T100000_C5 | gas5_ics3 | 15 | -1,104 | 7,368 [6,520, 8,034] | 38.3 | 0.84 | no | 0.99 | 3,595 / 2,604 |
| 109 | SL_Z4_R30_T100000_C5 | gas5_ics3 | 15 | -1,103 | 7,480 [6,495, 8,031] | 36.7 | 0.73 | no | 0.65 | 3,279 / 2,152 |
| 110 | SS_Z4_R30_T100000_C5 | gas1_ics1 | 9 | -383 | 6,018 [5,529, 6,514] | 155.5 | 0.74 | no | 0.98 | 1,669 / 2,805 |
| 111 | SS_Z4_R20_T100000_C5 | gas1_ics1 | 9 | -374 | 6,018 [5,553, 6,534] | 158.1 | 0.76 | no | 0.98 | 1,759 / 2,806 |
| 112 | SY_Z4_R30_T100000_C5 | gas1_ics1 | 9 | -362 | 6,131 [5,657, 6,633] | 162.5 | 0.73 | no | 1.00 | 1,773 / 2,685 |
| 113 | SY_Z4_R20_T100000_C5 | gas1_ics1 | 9 | -358 | 6,131 [5,653, 6,629] | 162.3 | 0.71 | no | 1.00 | 1,698 / 2,685 |
| 114 | SS_Z6_R30_T100000_C5 | gas1_ics1 | 9 | -318 | 6,074 [5,588, 6,568] | 159.4 | 0.75 | no | 1.00 | 1,669 / 2,875 |
| 115 | SL_Z4_R30_T100000_C5 | gas1_ics1 | 9 | -317 | 6,187 [5,695, 6,672] | 162.1 | 0.74 | no | 0.99 | 2,027 / 2,533 |
| 116 | SL_Z4_R20_T100000_C5 | gas1_ics1 | 9 | -316 | 6,187 [5,689, 6,666] | 161.8 | 0.74 | no | 0.99 | 2,057 / 2,519 |
| 117 | SS_Z4_R30_T150_C5 | gas1_ics1 | 9 | -306 | 5,962 [5,477, 6,463] | 153.1 | 0.74 | yes | 1.00 | 1,399 / 2,990 |
| 118 | SS_Z4_R20_T150_C5 | gas1_ics1 | 9 | -301 | 5,962 [5,474, 6,461] | 153.1 | 0.74 | yes | 1.00 | 1,399 / 2,990 |
| 119 | SY_Z4_R30_T150_C5 | gas1_ics1 | 9 | -298 | 6,018 [5,526, 6,508] | 154.0 | 0.74 | no | 1.00 | 1,639 / 2,842 |
| 120 | SY_Z4_R30_T100000_C2 | gas1_ics1 | 9 | -296 | 6,243 [5,754, 6,732] | 167.0 | 0.72 | yes | 1.00 | 1,698 / 2,797 |
| 121 | SS_Z6_R20_T100000_C5 | gas1_ics1 | 9 | -295 | 6,074 [5,592, 6,569] | 160.2 | 0.75 | no | 1.00 | 1,683 / 2,902 |
| 122 | SY_Z4_R20_T150_C5 | gas1_ics1 | 9 | -295 | 6,018 [5,523, 6,505] | 153.9 | 0.74 | no | 0.99 | 1,639 / 2,842 |
| 123 | SY_Z4_R20_T100000_C2 | gas1_ics1 | 9 | -292 | 6,243 [5,751, 6,729] | 167.0 | 0.72 | yes | 1.00 | 1,698 / 2,797 |
| 124 | SO_Z4_R20_T100000_C5 | gas1_ics1 | 9 | -292 | 5,737 [5,257, 6,230] | 145.1 | 0.68 | no | 0.99 | 1,550 / 2,373 |
| 125 | SS_Z4_R30_T100000_C2 | gas1_ics1 | 9 | -282 | 6,074 [5,598, 6,580] | 160.0 | 0.72 | yes | 0.98 | 1,549 / 2,840 |
| 126 | SS_Z4_R20_T100000_C2 | gas1_ics1 | 9 | -275 | 6,074 [5,597, 6,579] | 160.1 | 0.72 | yes | 0.98 | 1,549 / 2,840 |
| 127 | SO_Z4_R30_T100000_C5 | gas1_ics1 | 9 | -265 | 5,737 [5,261, 6,234] | 145.3 | 0.68 | no | 0.95 | 1,550 / 2,373 |
| 128 | SY_Z6_R30_T100000_C5 | gas1_ics1 | 9 | -264 | 6,299 [5,781, 6,766] | 165.9 | 0.74 | no | 0.99 | 1,877 / 2,804 |
| 129 | SL_Z4_R30_T100000_C2 | gas1_ics1 | 9 | -259 | 6,356 [5,839, 6,818] | 169.7 | 0.74 | no | 1.00 | 2,042 / 2,658 |
| 130 | SL_Z4_R20_T100000_C2 | gas1_ics1 | 9 | -259 | 6,299 [5,833, 6,812] | 169.3 | 0.74 | no | 1.00 | 2,027 / 2,645 |
| 131 | SL_Z6_R30_T100000_C5 | gas1_ics1 | 9 | -248 | 6,356 [5,870, 6,863] | 165.6 | 0.77 | no | 1.00 | 2,252 / 2,665 |
| 132 | SY_Z6_R20_T100000_C5 | gas1_ics1 | 9 | -247 | 6,299 [5,781, 6,765] | 166.3 | 0.76 | no | 0.99 | 1,952 / 2,832 |
| 133 | SS_Z6_R30_T150_C5 | gas1_ics1 | 9 | -240 | 5,962 [5,472, 6,455] | 152.8 | 0.76 | no | 1.00 | 1,609 / 2,950 |
| 134 | SL_Z4_R30_T150_C5 | gas1_ics1 | 9 | -238 | 6,074 [5,585, 6,565] | 155.7 | 0.77 | no | 0.99 | 2,013 / 2,651 |
| 135 | SL_Z4_R20_T150_C5 | gas1_ics1 | 9 | -238 | 6,074 [5,580, 6,561] | 155.4 | 0.77 | no | 0.99 | 2,013 / 2,651 |
| 136 | SL_Z6_R20_T100000_C5 | gas1_ics1 | 9 | -236 | 6,356 [5,868, 6,860] | 165.6 | 0.78 | no | 0.99 | 2,252 / 2,679 |
| 137 | SS_Z8_R30_T100000_C5 | gas1_ics1 | 9 | -235 | 6,131 [5,619, 6,601] | 159.9 | 0.74 | no | 0.82 | 1,698 / 2,836 |
| 138 | SS_Z6_R30_T100000_C2 | gas1_ics1 | 9 | -223 | 6,131 [5,663, 6,640] | 164.4 | 0.74 | yes | 0.98 | 1,683 / 2,856 |
| 139 | SS_Z6_R20_T150_C5 | gas1_ics1 | 9 | -220 | 5,962 [5,476, 6,457] | 153.3 | 0.77 | no | 1.00 | 1,624 / 2,950 |
| 140 | SL_Z8_R30_T100000_C5 | gas1_ics1 | 9 | -205 | 6,468 [5,976, 6,982] | 166.5 | 0.79 | no | 0.99 | 2,387 / 2,719 |
| 141 | SS_Z4_R30_T150_C2 | gas1_ics1 | 9 | -200 | 6,018 [5,506, 6,494] | 158.6 | 0.73 | yes | 0.96 | 1,370 / 3,040 |
| 142 | SY_Z8_R30_T100000_C5 | gas1_ics1 | 9 | -200 | 6,356 [5,859, 6,851] | 167.0 | 0.75 | no | 1.00 | 1,907 / 2,859 |
| 143 | SS_Z6_R20_T100000_C2 | gas1_ics1 | 9 | -199 | 6,131 [5,669, 6,644] | 165.3 | 0.74 | yes | 1.00 | 1,683 / 2,856 |
| 144 | SY_Z4_R30_T150_C2 | gas1_ics1 | 9 | -198 | 6,131 [5,616, 6,600] | 159.1 | 0.75 | no | 1.00 | 1,609 / 2,983 |
| 145 | SS_Z10_R30_T100000_C5 | gas1_ics1 | 9 | -196 | 6,187 [5,692, 6,669] | 164.6 | 0.74 | no | 0.81 | 1,743 / 2,811 |
| 146 | SY_Z4_R20_T150_C2 | gas1_ics1 | 9 | -195 | 6,131 [5,613, 6,597] | 159.0 | 0.75 | no | 1.00 | 1,609 / 2,983 |
| 147 | SS_Z8_R30_T150_C5 | gas1_ics1 | 9 | -195 | 6,018 [5,555, 6,535] | 157.6 | 0.78 | no | 0.97 | 1,818 / 2,870 |
| 148 | SS_Z8_R20_T100000_C5 | gas1_ics1 | 9 | -193 | 6,131 [5,645, 6,621] | 163.1 | 0.75 | no | 0.94 | 1,728 / 2,862 |
| 149 | SS_Z4_R20_T150_C2 | gas1_ics1 | 9 | -190 | 6,018 [5,504, 6,491] | 158.5 | 0.71 | yes | 0.97 | 1,385 / 2,862 |
| 150 | SL_Z6_R30_T100000_C2 | gas1_ics1 | 9 | -168 | 6,524 [6,053, 7,045] | 176.3 | 0.79 | no | 0.99 | 2,387 / 2,753 |
| 151 | SO_Z6_R20_T100000_C5 | gas1_ics1 | 9 | -162 | 5,793 [5,334, 6,307] | 149.3 | 0.72 | no | 0.98 | 1,624 / 2,550 |
| 152 | SL_Z8_R20_T100000_C5 | gas1_ics1 | 9 | -160 | 6,524 [6,009, 7,010] | 170.0 | 0.80 | no | 0.99 | 2,461 / 2,733 |
| 153 | SY_Z6_R30_T100000_C2 | gas1_ics1 | 9 | -159 | 6,412 [5,915, 6,896] | 175.3 | 0.75 | no | 1.00 | 1,922 / 2,905 |
| 154 | SS_Z8_R20_T150_C5 | gas1_ics1 | 9 | -157 | 6,074 [5,571, 6,546] | 159.4 | 0.78 | no | 1.00 | 1,848 / 2,884 |
| 155 | SL_Z6_R20_T100000_C2 | gas1_ics1 | 9 | -155 | 6,524 [6,054, 7,044] | 176.6 | 0.79 | no | 0.99 | 2,402 / 2,753 |
| 156 | SY_Z6_R30_T150_C5 | gas1_ics1 | 9 | -155 | 6,131 [5,647, 6,632] | 159.8 | 0.80 | no | 0.99 | 1,938 / 2,960 |
| 157 | SY_Z8_R20_T100000_C5 | gas1_ics1 | 9 | -151 | 6,356 [5,890, 6,877] | 170.8 | 0.78 | no | 1.00 | 2,101 / 2,872 |
| 158 | SS_Z10_R30_T150_C5 | gas1_ics1 | 9 | -147 | 6,074 [5,598, 6,576] | 159.3 | 0.80 | no | 0.97 | 1,878 / 2,994 |
| 159 | SL_Z10_R30_T100000_C5 | gas1_ics1 | 9 | -146 | 6,580 [6,098, 7,109] | 171.6 | 0.81 | no | 0.98 | 2,551 / 2,748 |
| 160 | SL_Z4_R30_T150_C2 | gas1_ics1 | 9 | -144 | 6,243 [5,724, 6,707] | 163.0 | 0.76 | no | 0.99 | 1,954 / 2,818 |
| 161 | SL_Z4_R20_T150_C2 | gas1_ics1 | 9 | -143 | 6,187 [5,721, 6,703] | 162.8 | 0.77 | no | 0.99 | 1,954 / 2,818 |
| 162 | SY_Z10_R30_T100000_C5 | gas1_ics1 | 9 | -142 | 6,468 [5,967, 6,959] | 172.2 | 0.78 | no | 1.00 | 2,191 / 2,875 |
| 163 | SO_Z4_R20_T150_C5 | gas1_ics1 | 9 | -142 | 5,793 [5,292, 6,271] | 145.0 | 0.68 | yes | 0.99 | 1,444 / 2,490 |
| 164 | SY_Z6_R20_T100000_C2 | gas1_ics1 | 9 | -140 | 6,412 [5,917, 6,897] | 176.0 | 0.76 | no | 1.00 | 1,937 / 2,905 |
| 165 | SY_Z6_R20_T150_C5 | gas1_ics1 | 9 | -139 | 6,131 [5,648, 6,631] | 160.1 | 0.80 | no | 0.99 | 1,938 / 2,960 |
| 166 | SO_Z4_R30_T150_C5 | gas1_ics1 | 9 | -139 | 5,793 [5,295, 6,274] | 145.2 | 0.68 | yes | 1.00 | 1,444 / 2,490 |
| 167 | SS_Z10_R20_T100000_C5 | gas1_ics1 | 9 | -139 | 6,187 [5,715, 6,687] | 167.8 | 0.75 | no | 0.85 | 1,788 / 2,836 |
| 168 | SS_Z8_R30_T100000_C2 | gas1_ics1 | 9 | -139 | 6,187 [5,700, 6,679] | 165.1 | 0.76 | no | 1.00 | 1,758 / 2,939 |
| 169 | SS_Z6_R30_T150_C2 | gas1_ics1 | 9 | -134 | 6,018 [5,503, 6,486] | 155.1 | 0.76 | no | 0.96 | 1,609 / 2,987 |
| 170 | SO_Z4_R20_T100000_C2 | gas1_ics1 | 9 | -121 | 5,793 [5,306, 6,281] | 147.1 | 0.66 | yes | 1.00 | 1,385 / 2,431 |
| 171 | SY_Z8_R30_T100000_C2 | gas1_ics1 | 9 | -119 | 6,468 [5,990, 6,981] | 175.0 | 0.76 | no | 1.00 | 1,998 / 2,920 |
| 172 | SO_Z4_R30_T100000_C2 | gas1_ics1 | 9 | -119 | 5,793 [5,310, 6,285] | 147.3 | 0.66 | yes | 1.00 | 1,385 / 2,431 |
| 173 | SO_Z6_R20_T100000_C2 | gas1_ics1 | 9 | -113 | 5,906 [5,400, 6,370] | 153.3 | 0.68 | no | 0.98 | 1,490 / 2,510 |
| 174 | SL_Z8_R30_T100000_C2 | gas1_ics1 | 9 | -113 | 6,637 [6,155, 7,160] | 176.9 | 0.81 | no | 0.99 | 2,522 / 2,835 |
| 175 | SS_Z6_R20_T150_C2 | gas1_ics1 | 9 | -112 | 6,018 [5,509, 6,490] | 155.8 | 0.77 | no | 0.97 | 1,639 / 3,000 |
| 176 | SS_Z10_R30_T100000_C2 | gas1_ics1 | 9 | -104 | 6,243 [5,778, 6,755] | 169.9 | 0.77 | no | 0.98 | 1,878 / 2,914 |
| 177 | SL_Z10_R20_T100000_C5 | gas1_ics1 | 9 | -102 | 6,637 [6,107, 7,114] | 173.3 | 0.80 | no | 0.98 | 2,581 / 2,761 |
| 178 | SS_Z8_R20_T100000_C2 | gas1_ics1 | 9 | -97 | 6,187 [5,717, 6,692] | 167.3 | 0.77 | no | 1.00 | 1,788 / 2,966 |
| 179 | SY_Z10_R20_T100000_C5 | gas1_ics1 | 9 | -94 | 6,468 [5,981, 6,968] | 174.7 | 0.80 | no | 0.99 | 2,251 / 2,901 |
| 180 | SS_Z10_R20_T150_C5 | gas1_ics1 | 9 | -93 | 6,131 [5,630, 6,604] | 162.5 | 0.79 | no | 0.99 | 1,983 / 2,844 |
| 181 | SL_Z6_R30_T150_C5 | gas1_ics1 | 9 | -90 | 6,243 [5,774, 6,767] | 161.1 | 0.82 | no | 1.00 | 2,373 / 2,729 |
| 182 | SL_Z8_R20_T100000_C2 | gas1_ics1 | 9 | -85 | 6,637 [6,160, 7,163] | 177.9 | 0.81 | no | 0.99 | 2,567 / 2,834 |
| 183 | SS_Z8_R30_T150_C2 | gas1_ics1 | 9 | -83 | 6,074 [5,568, 6,549] | 158.1 | 0.75 | no | 0.99 | 1,595 / 2,975 |
| 184 | SO_Z6_R20_T150_C5 | gas1_ics1 | 9 | -82 | 5,849 [5,346, 6,322] | 148.8 | 0.69 | no | 0.98 | 1,460 / 2,592 |
| 185 | SY_Z8_R20_T100000_C2 | gas1_ics1 | 9 | -82 | 6,468 [6,000, 6,988] | 176.7 | 0.79 | no | 0.99 | 2,222 / 2,919 |
| 186 | SL_Z6_R20_T150_C5 | gas1_ics1 | 9 | -80 | 6,243 [5,775, 6,767] | 161.4 | 0.82 | no | 0.99 | 2,373 / 2,770 |
| 187 | SS_Z10_R20_T100000_C2 | gas1_ics1 | 9 | -72 | 6,299 [5,818, 6,786] | 174.3 | 0.78 | no | 0.99 | 1,937 / 2,953 |
| 188 | SY_Z10_R30_T100000_C2 | gas1_ics1 | 9 | -67 | 6,580 [6,104, 7,095] | 180.2 | 0.79 | no | 0.99 | 2,222 / 2,978 |
| 189 | SL_Z10_R30_T100000_C2 | gas1_ics1 | 9 | -55 | 6,805 [6,285, 7,293] | 182.3 | 0.82 | no | 0.99 | 2,702 / 2,864 |
| 190 | SS_Z10_R30_T150_C2 | gas1_ics1 | 9 | -50 | 6,131 [5,633, 6,612] | 161.6 | 0.77 | no | 0.99 | 1,790 / 2,936 |
| 191 | SS_Z8_R20_T150_C2 | gas1_ics1 | 9 | -48 | 6,074 [5,601, 6,578] | 161.2 | 0.76 | no | 1.00 | 1,625 / 2,988 |
| 192 | SY_Z8_R30_T150_C5 | gas1_ics1 | 9 | -46 | 6,299 [5,803, 6,788] | 168.0 | 0.84 | no | 0.99 | 2,223 / 3,056 |
| 193 | SY_Z6_R30_T150_C2 | gas1_ics1 | 9 | -30 | 6,243 [5,734, 6,720] | 168.5 | 0.80 | no | 0.99 | 1,894 / 3,075 |
| 194 | SY_Z8_R20_T150_C5 | gas1_ics1 | 9 | -14 | 6,299 [5,812, 6,795] | 169.2 | 0.84 | no | 0.99 | 2,253 / 3,015 |
| 195 | SY_Z6_R20_T150_C2 | gas1_ics1 | 9 | -13 | 6,243 [5,736, 6,721] | 168.9 | 0.80 | no | 0.99 | 1,909 / 3,075 |
| 196 | SL_Z10_R20_T100000_C2 | gas1_ics1 | 9 | -9 | 6,805 [6,306, 7,309] | 184.8 | 0.83 | no | 0.99 | 2,762 / 2,877 |
| 197 | SY_Z10_R20_T100000_C2 | gas1_ics1 | 9 | -8 | 6,637 [6,133, 7,118] | 183.9 | 0.79 | no | 0.99 | 2,238 / 3,003 |
| 198 | SL_Z8_R30_T150_C5 | gas1_ics1 | 9 | -0 | 6,468 [5,965, 6,963] | 168.4 | 0.86 | no | 0.98 | 2,703 / 2,840 |
| 199 | SL_Z6_R30_T150_C2 | gas1_ics1 | 9 | +0 | 6,412 [5,904, 6,898] | 167.9 | 0.81 | no | 0.99 | 2,299 / 2,910 |
| 200 | SS_Z10_R20_T150_C2 | gas1_ics1 | 9 | +9 | 6,187 [5,674, 6,650] | 165.5 | 0.78 | no | 1.00 | 1,894 / 2,948 |
| 201 | SL_Z6_R20_T150_C2 | gas1_ics1 | 9 | +11 | 6,412 [5,904, 6,897] | 168.0 | 0.81 | no | 1.00 | 2,314 / 2,897 |
| 202 | SY_Z10_R30_T150_C5 | gas1_ics1 | 9 | +13 | 6,356 [5,880, 6,869] | 169.9 | 0.85 | no | 0.99 | 2,373 / 3,044 |
| 203 | SO_Z4_R20_T150_C2 | gas1_ics1 | 9 | +24 | 5,793 [5,329, 6,310] | 146.9 | 0.68 | yes | 0.99 | 1,415 / 2,539 |
| 204 | SL_Z8_R20_T150_C5 | gas1_ics1 | 9 | +25 | 6,468 [5,970, 6,966] | 169.0 | 0.86 | no | 0.98 | 2,732 / 2,852 |
| 205 | SO_Z4_R30_T150_C2 | gas1_ics1 | 9 | +27 | 5,849 [5,333, 6,314] | 147.1 | 0.67 | yes | 0.99 | 1,400 / 2,539 |
| 206 | SY_Z8_R30_T150_C2 | gas1_ics1 | 9 | +38 | 6,356 [5,861, 6,850] | 172.3 | 0.83 | no | 0.99 | 2,179 / 3,078 |
| 207 | SO_Z6_R20_T150_C2 | gas1_ics1 | 9 | +42 | 5,849 [5,363, 6,339] | 150.1 | 0.71 | no | 0.99 | 1,535 / 2,603 |
| 208 | SL_Z10_R30_T150_C5 | gas1_ics1 | 9 | +57 | 6,580 [6,050, 7,057] | 170.0 | 0.86 | no | 0.98 | 2,838 / 2,853 |
| 209 | SY_Z10_R20_T150_C5 | gas1_ics1 | 9 | +61 | 6,412 [5,896, 6,882] | 172.1 | 0.85 | no | 0.99 | 2,403 / 3,044 |
| 210 | SY_Z8_R20_T150_C2 | gas1_ics1 | 9 | +73 | 6,412 [5,895, 6,880] | 175.8 | 0.83 | no | 0.99 | 2,239 / 3,104 |
| 211 | SL_Z8_R30_T150_C2 | gas1_ics1 | 9 | +83 | 6,580 [6,074, 7,074] | 174.5 | 0.85 | no | 0.99 | 2,584 / 2,980 |
| 212 | SY_Z10_R30_T150_C2 | gas1_ics1 | 9 | +86 | 6,468 [5,958, 6,950] | 175.6 | 0.84 | no | 0.99 | 2,329 / 3,120 |
| 213 | SL_Z10_R20_T150_C5 | gas1_ics1 | 9 | +93 | 6,580 [6,058, 7,062] | 171.2 | 0.87 | no | 0.98 | 2,868 / 2,866 |
| 214 | SL_Z8_R20_T150_C2 | gas1_ics1 | 9 | +108 | 6,580 [6,100, 7,098] | 176.5 | 0.86 | no | 0.98 | 2,644 / 3,007 |
| 215 | SY_Z10_R20_T150_C2 | gas1_ics1 | 9 | +138 | 6,468 [5,986, 6,973] | 178.5 | 0.85 | no | 0.99 | 2,389 / 3,120 |
| 216 | SL_Z10_R30_T150_C2 | gas1_ics1 | 9 | +152 | 6,693 [6,179, 7,187] | 178.1 | 0.86 | no | 1.00 | 2,749 / 3,021 |
| 217 | SL_Z10_R20_T150_C2 | gas1_ics1 | 9 | +185 | 6,693 [6,197, 7,202] | 179.7 | 0.87 | no | 0.98 | 2,779 / 3,021 |
| 218 | gll_iem_v07 | - | 6 | +850 | 3,150 [2,614, 3,663] | 36.4 | 0.50 | yes | 0.97 | 110 / 1,473 |
