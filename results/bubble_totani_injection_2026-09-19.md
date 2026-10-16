# Injection test of Totani's bubble construction — 2026-09-19

Kev: "Do the injection test." `scripts/bubble_totani_injection.py` →
`fits/bubble_totani_injection_2026-09-19.json`; 148 s; peak RSS 6.0 GB.
His configuration throughout (GALPROP SLZ6R30T150C2 gas + ICS, two Loop I
shells, isotropic, point sources, GCE in the construction fit).

## The question

The standard gate holds the templates fixed and asks whether the fitter
recovers a planted halo — it does (0.99–1.00 in every bin,
`bubble_totani_2026-09-19`). But his positive and negative templates are
*built from the data at 4.3 GeV*: a halo present there goes into them, and
thereafter a template with a free normalisation per bin carries the halo's
morphology. So: plant a halo, **rebuild the templates as he builds them**,
fit, and ask what fraction comes back, per bin.

Injections at 0, 1 and 2 times the fitted 20.8 GeV amplitude, with two
spectra: **flat** (the same amplitude in every bin) and **Totani-like** (his
own fitted spectrum, positive part, zero at 1.5–2.6 GeV, peaking at 20 GeV —
does his construction let his own signal through?). Recovered/injected
slope in covariance form; alongside, the same injections with the templates
frozen at their a = 0 form.

## Result: recovered / injected slope

| | 1.5 | 2.6 | 4.3 | 7.3 | 12.3 | 20.8 | 35 | 59 | 100 GeV |
|---|---|---|---|---|---|---|---|---|---|
| flat spectrum, templates **rebuilt** | **−5.2** | **−0.8** | **0.45** | 0.79 | 0.95 | 1.00 | 0.99 | 1.00 | 1.00 |
| flat spectrum, templates fixed | −0.03 | 0.74 | 0.89 | 0.97 | 1.00 | 1.01 | 0.99 | 1.00 | 1.00 |
| Totani-like, templates **rebuilt** | — | — | **0.41** | 0.82 | 0.90 | **0.93** | 0.91 | 0.92 | 0.87 |
| Totani-like, templates fixed | — | — | 0.97 | 1.01 | 1.00 | 1.01 | 0.99 | 1.01 | 1.00 |

(— : nothing injected in that bin.)

## Reading

1. **At the construction energy the templates eat the halo.** 55–59% of a
   planted halo at 4.3 GeV disappears into the positive template, whatever
   the spectrum. Expected from his own words ("even if a halo component is
   present near this energy bin, it will be absorbed"), now measured.

2. **The 20 GeV peak survives — 93% of it.** A Totani-like halo comes back
   at 0.90–0.93 from 12 to 59 GeV. His construction does not destroy his
   signal at the energies where he reports it. That is the fair statement
   in his favour.

3. **But the construction manufactures his spectral shape.** A halo with a
   *flat* amplitude in every bin is recovered at 0.45 at 4.3 GeV, 0.79 at
   7.3, 0.95 at 12.3 and 1.00 from 20 GeV up — so a truly flat halo would
   come out *rising* from 4 to 20 GeV. And at 1.5–2.6 GeV the recovered
   amplitude goes the wrong way: slope **−5.2** and −0.8 — planting a
   positive halo makes the fitted low-energy halo *more negative*, because
   the rebuilt templates absorb it at 4.3 GeV and then, scaled freely at
   1.5 GeV, over-subtract. The "negative below 3 GeV, rising to a 20 GeV
   peak" shape is what this construction does to any halo, including a flat
   one. The shape is not evidence about the halo's spectrum.

4. **The 1.5 GeV bin is not a measurement in this configuration at all**:
   even with the templates frozen, the slope there is −0.03 — the halo
   template is fully degenerate with the residual templates at that energy
   (they were built from the model's 4.3 GeV residual, which at 1.5 GeV is
   dominated by the GALPROP gas over-prediction). Foster's N/S asymmetry
   headline is built on exactly the low-energy bins of this kind of fit;
   under his 2012-set fits the north median is negative below 5 GeV. This
   test says those bins do not carry halo information.

5. **Where the number is trustworthy**: 12–100 GeV, to 7–13%, *given* the
   templates. What the templates are built from — the model's own residual
   — is the question, not the fitter.

## What it changes

- His 6,412 ± 500 at 20.8 GeV should be read as 6,412 / 0.93 ≈ 6,900 of
  "true" halo-shaped counts *if* the halo is real, or as what his model
  leaves at 20.8 GeV that is not at 4.3 GeV *if* it is not. The test cannot
  tell those apart; it only says his construction passes 93% of the former.
- His spectral shape below 12 GeV, and every argument built on it — his
  "consistent with zero below 2 GeV", the 20 GeV *peak* as opposed to a
  plateau, and (in the same class of fit) Foster's low-energy N/S
  asymmetry — is a property of the construction.
- For the ranking programme: the residual-template construction should be
  a separate study, not a candidate in the structure selection; it wins any
  likelihood contest at 4.3 GeV by construction and reshapes the halo
  spectrum by construction.

Analysis and text: Claude Fable 5.1 (`claude-fable-5-1[1m]`), under Kev's
direction; 165 tests passing.
