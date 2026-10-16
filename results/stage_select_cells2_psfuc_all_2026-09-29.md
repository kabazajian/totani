# Structure selection and scoring — cells2_psfuc_all (2026-09-29)

Bins used for selection AND inference: [0, 1, 2, 3, 4, 5, 6, 7, 8] (K = 9). Fraction-of-rows rule, cutoff 0.8; BIC p ln n with n = cells × K. Voters: 16 stratified 2012 models + A2020 per decomposition; Galp21 rows scored, not voting. 'favor' = TS > 4 at 20.8 GeV, positive, all four quadrants positive. Q⁺ over the K bins with the binomial χ² mixture; Z₂₀₁₈ = Φ⁻¹(1 − p/2).

## gas1_ics1, residual

Selected: **bubble_totani → loopI_shells+nps → erosita_bright_NS**

- 2012 subset: favor 17/17; halo20 median 6,631 [6,299–8,704]; Z₂₀₁₈ median 28.0 [22.9–35.6]; ⟨σv⟩ median +5.6, 2σ+ 17/17
- Galp21: favor 12/12; halo20 median 7,111 [6,567–7,535]; Z₂₀₁₈ median 36.8 [28.3–40.6]; ⟨σv⟩ median +7.2, 2σ+ 12/12
- own paths, bubble chosen: {'2012': {'bubble_totani': 17}, 'G21': {'bubble_totani': 12}}

- step 1: clears {'bubble_flat': 1.0, 'bubble_structured': 1.0, 'bubble_totani': 1.0}; flat − structured improvement per row median +1129 (flat larger in 15/17 rows)

## gas1_ics1, no residual

Selected: **bubble_flat → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic**

- 2012 subset: favor 17/17; halo20 median 5,735 [4,763–7,404]; Z₂₀₁₈ median 28.7 [16.2–44.5]; ⟨σv⟩ median +5.5, 2σ+ 17/17
- Galp21: favor 12/12; halo20 median 5,688 [5,238–6,100]; Z₂₀₁₈ median 28.6 [19.9–37.5]; ⟨σv⟩ median +5.7, 2σ+ 12/12
- own paths, bubble chosen: {'2012': {'bubble_flat': 15, 'bubble_structured': 2}, 'G21': {'bubble_flat': 12}}

- step 1: clears {'bubble_flat': 1.0, 'bubble_structured': 1.0}; flat − structured improvement per row median +1129 (flat larger in 15/17 rows)

## gas5_ics3, residual

Selected: **bubble_totani → loopI_shells+nps → erosita_unif_NS**

- 2012 subset: favor 17/17; halo20 median 7,436 [7,022–8,704]; Z₂₀₁₈ median 14.1 [9.8–24.9]; ⟨σv⟩ median +7.4, 2σ+ 17/17
- Galp21: favor 12/12; halo20 median 7,468 [7,015–7,925]; Z₂₀₁₈ median 15.0 [13.8–18.7]; ⟨σv⟩ median +7.6, 2σ+ 12/12
- own paths, bubble chosen: {'2012': {'bubble_totani': 17}, 'G21': {'bubble_totani': 12}}

- step 1: clears {'bubble_flat': 1.0, 'bubble_structured': 1.0, 'bubble_totani': 1.0}; flat − structured improvement per row median +1910 (flat larger in 17/17 rows)

## gas5_ics3, no residual

Selected: **bubble_flat → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic**

- 2012 subset: favor 10/17; halo20 median 6,205 [5,363–7,404]; Z₂₀₁₈ median 11.5 [5.1–33.5]; ⟨σv⟩ median +6.0, 2σ+ 17/17
- Galp21: favor 12/12; halo20 median 6,695 [6,092–6,882]; Z₂₀₁₈ median 18.9 [13.8–23.3]; ⟨σv⟩ median +7.2, 2σ+ 12/12
- own paths, bubble chosen: {'2012': {'bubble_flat': 17}, 'G21': {'bubble_flat': 12}}

- step 1: clears {'bubble_flat': 1.0, 'bubble_structured': 1.0}; flat − structured improvement per row median +1910 (flat larger in 17/17 rows)

