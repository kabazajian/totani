# Structure selection and scoring — cells2_psfuc_resp (2026-09-29)

Bins used for selection AND inference: [3, 4, 5, 6, 7, 8] (K = 6). Fraction-of-rows rule, cutoff 0.8; BIC p ln n with n = cells × K. Voters: 16 stratified 2012 models + A2020 per decomposition; Galp21 rows scored, not voting. 'favor' = TS > 4 at 20.8 GeV, positive, all four quadrants positive. Q⁺ over the K bins with the binomial χ² mixture; Z₂₀₁₈ = Φ⁻¹(1 − p/2).

## gas1_ics1, residual

Selected: **bubble_totani → loopI_shells+nps → erosita_bright_NS**

- 2012 subset: favor 17/17; halo20 median 6,631 [6,299–8,704]; Z₂₀₁₈ median 27.3 [19.8–28.2]; ⟨σv⟩ median +6.6, 2σ+ 17/17
- Galp21: favor 12/12; halo20 median 7,111 [6,567–7,535]; Z₂₀₁₈ median 31.9 [27.2–33.4]; ⟨σv⟩ median +7.3, 2σ+ 12/12
- own paths, bubble chosen: {'2012': {'bubble_totani': 17}, 'G21': {'bubble_totani': 12}}

- step 1: clears {'bubble_flat': 1.0, 'bubble_structured': 1.0, 'bubble_totani': 1.0}; flat − structured improvement per row median -165 (flat larger in 2/17 rows)

## gas1_ics1, no residual

Selected: **bubble_structured → erosita_unif_NS → loopI_haslam_raw**

- 2012 subset: favor 17/17; halo20 median 4,097 [3,184–5,810]; Z₂₀₁₈ median 12.6 [7.5–17.4]; ⟨σv⟩ median +3.1, 2σ+ 17/17
- Galp21: favor 12/12; halo20 median 4,053 [3,639–4,420]; Z₂₀₁₈ median 12.1 [9.1–14.9]; ⟨σv⟩ median +3.0, 2σ+ 12/12
- own paths, bubble chosen: {'2012': {'bubble_structured': 15, 'bubble_flat': 2}, 'G21': {'bubble_structured': 12}}

- step 1: clears {'bubble_flat': 1.0, 'bubble_structured': 1.0}; flat − structured improvement per row median -165 (flat larger in 2/17 rows)

## gas5_ics3, residual

Selected: **bubble_totani → loopI_shells+nps**

- 2012 subset: favor 17/17; halo20 median 7,440 [7,082–8,704]; Z₂₀₁₈ median 15.6 [11.1–22.0]; ⟨σv⟩ median +7.5, 2σ+ 17/17
- Galp21: favor 12/12; halo20 median 7,527 [7,015–8,029]; Z₂₀₁₈ median 15.2 [13.4–16.2]; ⟨σv⟩ median +7.7, 2σ+ 12/12
- own paths, bubble chosen: {'2012': {'bubble_totani': 17}, 'G21': {'bubble_totani': 12}}

- step 1: clears {'bubble_flat': 1.0, 'bubble_structured': 1.0, 'bubble_totani': 1.0}; flat − structured improvement per row median +252 (flat larger in 16/17 rows)

## gas5_ics3, no residual

Selected: **bubble_flat → erosita_outer_NS**

- 2012 subset: favor 17/17; halo20 median 6,154 [5,248–7,588]; Z₂₀₁₈ median 11.2 [6.1–21.3]; ⟨σv⟩ median +6.1, 2σ+ 17/17
- Galp21: favor 12/12; halo20 median 6,740 [6,153–7,281]; Z₂₀₁₈ median 15.1 [12.4–18.1]; ⟨σv⟩ median +6.8, 2σ+ 12/12
- own paths, bubble chosen: {'2012': {'bubble_flat': 16, 'bubble_structured': 1}, 'G21': {'bubble_flat': 12}}

- step 1: clears {'bubble_flat': 1.0, 'bubble_structured': 1.0}; flat − structured improvement per row median +252 (flat larger in 16/17 rows)

