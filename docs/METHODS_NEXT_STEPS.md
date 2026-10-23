# DKFZ Boundary-ID: next-stage analysis protocol (23 October 2026)

## Status and provenance
- Inputs: genuine public NSCLC Visium samples P5 (GSM8855707) and P12 (GSM8855714), previously processed and aligned to tissue coordinates. The uploaded original files are still available separately.
- Latest outputs: `edge_summary.csv`, `P5_edge_diagnostic.csv.gz`, and `P12_edge_diagnostic.csv.gz` are **observational edge-level diagnostics**, calculated from the expression-only K-means baseline and spatial adjacency graph.
- `identifiability_counterexample.png` is a **synthetic theoretical illustration**, not a patient-data output.
- No claims of validated biological/technical boundary classification, calibrated false discovery rate, or causal attribution.

## Immediate finding from genuine data
- P5: 15,700 adjacency edges, boundary fraction 23.465%; boundary frequency 15.64% in the lowest quartile of local log-UMI jump and 37.55% in the highest quartile.
- P12: 5,789 adjacency edges, boundary fraction 12.610%; lowest/highest quartile boundary frequencies 4.63% / 27.21%.
- These descriptive results are compatible with depth-associated cluster boundaries but **cannot identify the reason**: local RNA content, tissue composition, capture efficiency, and image morphology may vary together.
- Spots/edges are correlated within tissue; do not apply independent-edge statistical significance tests.

## Scientific question and estimand
Let G=(V,E) be the physically grounded spot-adjacency graph, Y_ig the raw gene-count matrix, s_i technical capture/sequencing exposure, mu_ig biological mean abundance, and z_i latent biological molecular state.

Poisson/NB expectation: E[Y_ig | s_i,mu_ig] = s_i * mu_ig.
For an edge e=(i,j), the *biological boundary* target is B*_e = I(z_i != z_j).
This latent target is **not generically identifiable from Y alone** if s_i and all-gene biology change by reciprocally compensating multiplicative factors. Explicit restrictions, anchors, orthogonal modalities, replicates, or interventions are needed.

The readily estimable quantity from data is B_e(A,m,r), an ALGORITHM-INDUCED crossing indicator for analysis pipeline A, nuisance-handling option m, and randomized perturbation r. We should report its empirical stability as an algorithm-sensitivity profile, NOT as P(B*_e=1|Y) or a formal identification bound.

## Proposed phases (and required data)
1. **Already done**: integrity QC, RNA-only filter, spatial graph, baseline transcriptomic clusters, global thinning, and depth association; latest edge-level diagnostic supplied.
2. **Next computation with current files only**: edge persistence from 10-20 independent perturbations for each of 3 count-thinning regimes (50%, 25%, equalized); compare 2 normalization choices and multiple clustering resolutions. Export per-edge q_e(m) and min/max over predeclared settings as **sensitivity ranges**, not confidence intervals. Use cached sparse matrix, fixed HVGs, and reproducible RNG seeds; chunk work by specimen/condition.
3. **Method development, current files + synthetic simulation**: explicitly specify technical capture field, biological smooth-gradient field, true biological discontinuity, noise and cell mixing. Include six scenarios: discrete boundary, pure technical gradient, coupled biology/depth, smooth biological field, observational equivalence, and rare biological state. Derive what is and is not identifiable.
4. **Spatially valid testing**: develop constrained null preserving target spatial covariance and nuisance depth profiles. Verify false-positive control in simulations *before* using p-values on tissue specimens.
5. **Independent validation (additional inputs needed)**: pathology annotations made without consulting clusters, independent molecular markers/DNA profiles, and at least several more patient specimens; masked-label or held-out-patient validation.

## Precautions
- A spatial graph produced from nearest neighbors is not an independent sample of edges.
- Count thinning decreases measurement information; a decreased ARI cannot be interpreted as artifact proportion.
- Regressing out RNA depth can remove genuine biology.
- Method benchmarking must avoid reusing held-out pathology labels to tune model configurations.
- Two patients (P5 and P12) are insufficient to quantify clinical generalizability.

## What user needs to do now
Nothing to start step 2 or step 3. The originals and previous analysis outputs are already available. If the user wants to move to phase 5 later, upload additional independent pathology annotations/DNA spatial data and more patients.
