# Structure selection and scoring — cells1_psfuc_resp (2026-09-29)

Bins used for selection AND inference: [3, 4, 5, 6, 7, 8] (K = 6). Fraction-of-rows rule, cutoff 0.8; BIC p ln n with n = cells × K. Voters: 16 stratified 2012 models + A2020 per decomposition; Galp21 rows scored, not voting. 'favor' = TS > 4 at 20.8 GeV, positive, all four quadrants positive. Q⁺ over the K bins with the binomial χ² mixture; Z₂₀₁₈ = Φ⁻¹(1 − p/2).

## gas1_ics1, residual

Selected: **bubble_totani → erosita_bright_NS → loopI_shells+nps**

- 2012 subset: favor 17/17; halo20 median 6,907 [6,515–8,910]; Z₂₀₁₈ median 28.8 [21.9–30.3]; ⟨σv⟩ median +6.9, 2σ+ 17/17
- Galp21: favor 1/1; halo20 median 7,220 [7,220–7,220]; Z₂₀₁₈ median 30.0 [30.0–30.0]; ⟨σv⟩ median +7.0, 2σ+ 1/1
- own paths, bubble chosen: {'2012': {'bubble_totani': 17}, 'G21': {'bubble_totani': 1}}

- step 1: clears {'bubble_flat': 1.0, 'bubble_structured': 1.0, 'bubble_totani': 1.0}; flat − structured improvement per row median -188 (flat larger in 2/17 rows)

## gas1_ics1, no residual

Selected: **bubble_structured → erosita_unif_NS → loopI_haslam_raw**

- 2012 subset: favor 17/17; halo20 median 4,454 [3,597–6,194]; Z₂₀₁₈ median 14.8 [9.5–19.6]; ⟨σv⟩ median +3.5, 2σ+ 17/17
- Galp21: favor 1/1; halo20 median 4,478 [4,478–4,478]; Z₂₀₁₈ median 14.0 [14.0–14.0]; ⟨σv⟩ median +3.4, 2σ+ 1/1
- own paths, bubble chosen: {'2012': {'bubble_structured': 15, 'bubble_flat': 2}, 'G21': {'bubble_structured': 1}}

- step 1: clears {'bubble_flat': 1.0, 'bubble_structured': 1.0}; flat − structured improvement per row median -188 (flat larger in 2/17 rows)

## gas5_ics3, residual

Selected: **bubble_totani → loopI_shells+nps**

- 2012 subset: favor 17/17; halo20 median 7,478 [7,118–8,910]; Z₂₀₁₈ median 17.1 [12.4–24.5]; ⟨σv⟩ median +7.8, 2σ+ 17/17
- Galp21: —
- own paths, bubble chosen: {'2012': {'bubble_totani': 17}}

- step 1: clears {'bubble_flat': 1.0, 'bubble_structured': 1.0, 'bubble_totani': 1.0}; flat − structured improvement per row median +234 (flat larger in 16/17 rows)

## gas5_ics3, no residual

Selected: **bubble_flat → erosita_outer_NS**

- 2012 subset: favor 17/17; halo20 median 6,281 [5,516–7,793]; Z₂₀₁₈ median 12.0 [7.0–22.0]; ⟨σv⟩ median +6.2, 2σ+ 17/17
- Galp21: —
- own paths, bubble chosen: {'2012': {'bubble_flat': 16, 'bubble_structured': 1}}

- step 1: clears {'bubble_flat': 1.0, 'bubble_structured': 1.0}; flat − structured improvement per row median +234 (flat larger in 16/17 rows)

