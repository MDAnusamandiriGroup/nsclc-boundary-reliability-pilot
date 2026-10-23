# DKFZ PRELUNG-compatible methodological proof of concept (exploratory)

**Question.** Do expression-derived molecular-state boundaries remain stable as sequencing depth is experimentally perturbed, or are the inferred states closely related to depth variation and spatial structure?

**Data.** Two public NSCLC Visium FFPE tissue specimens from GEO GSE292299: P5, 5,388 spots; P12, 2,131 spots. Original spatial coordinates and histology were matched to all barcodes. Each count matrix includes 18,085 RNA features; 35 antibody-capture features were excluded.

**Method.** Separate RNA count-normalized PCA/K-means clustering (k=4) per specimen; 3 stochastic repetitions each for 50% and 25% count retention and spot-specific equal-depth thinning; adjusted Rand index (ARI) against original clusters; spatial spot graph; descriptive PC1–depth correlation; edge-based depth sensitivity; aggressive depth-regressed representation used as stress test. Spatial coordinates are used for descriptive maps and spatial graph diagnostics, not as clustering inputs.

**Results.** PC1 remained strongly correlated with log RNA depth (Spearman ρ=0.868 in P5, −0.951 in P12). Mean ARI was 0.834/0.875 under half-depth thinning, 0.705/0.701 under quarter-depth thinning, and 0.721/0.531 under spot-specific equal-depth thinning (P5/P12). Labels from depth-regressed expression had ARI 0.412/0.368 relative to the original representation, showing strong sensitivity to a potential nuisance adjustment. Boundary edges differed more in local RNA depth than non-boundary edges (descriptive depth-jump AUROC 0.635/0.728). RNA depth itself was spatially autocorrelated (Moran's I 0.833/0.858).

**Interpretation.** These observations expose a potentially important limitation of conventional molecular-state detection: apparently structured and spatially coherent expression clusters can be highly sensitive to the handling of RNA depth. This is **not evidence that the clusters are artifacts**; depth can reflect real tissue biology, and neither selected clustering nor simple permutation nulls constitute validated biological inference.

**Novelty target, not yet achieved.** A calibrated conditional *boundary-reliability* framework that distinguishes identifiable biological transitions from depth-dependent and spatially correlated alternatives; to be tested in simulation and validated using multiple biological specimens and independent morphology. Current output is a methodological pilot, not a completed algorithm or biomarker claim.

**Files.** `README.md` describes the full analysis, caveats and reproduction; PNG figures, CSV results and Python scripts accompany this summary.
