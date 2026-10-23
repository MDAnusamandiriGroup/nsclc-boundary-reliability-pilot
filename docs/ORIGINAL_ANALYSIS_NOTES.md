# DKFZ NSCLC spatial transcriptomics — depth/parameter robustness proof of concept

**Status:** exploratory computational pilot on real public Visium data. **Not** a new validated method, clinical diagnostic, histological classification, or inferential proof that observed domains are biologically real.

## Data and provenance

Public GEO dataset GSE292299 (NSCLC Visium FFPE). Two separate lung NSCLC specimens:

- P5: `GSM8855707_NSCLC_P5_filtered_feature_bc_matrix.h5` and `GSM8855707_NSCLC_P5_spatial.tar.gz`
- P12: `GSM8855714_NSCLC_P12_filtered_feature_bc_matrix.h5` and `GSM8855714_NSCLC_P12_spatial.tar.gz`

GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE292299

**Raw data intentionally excluded** from the analysis ZIP. Obtain original files via GEO or provide the four originals; do not silently substitute artificial data. All analyses below were run on the supplied original H5 and spatial archives. Archive paths are read in memory and are not extracted into working directories.

## Data integrity checks

| Field | P5 | P12 |
| --- | ---: | ---: |
| Matched in-tissue barcodes | 5,388 | 2,131 |
| Gene Expression RNA features | 18,085 | 18,085 |
| Antibody Capture features excluded | 35 | 35 |
| Spots matched uniquely to coordinates | All | All |
| Median RNA UMI counts / spot | 1,985 | 3,231 |
| Median detected RNA genes / spot | 1,453 | 2,115 |
| Local tissue-neighbour edges | 15,700 | 5,789 |

Histological image, spatial coordinates, and scale factors are read from each original spatial archive.

## Computational procedure

1. Read HDF5 count matrix as sparse spot × feature; **retain Gene Expression only**, removing 35 `Antibody Capture` features (the former must not be conflated with RNA).
2. Align in-tissue barcodes and spatial spot coordinates, and verify exact set equality.
3. QC: spot RNA counts, detected genes, fraction of MT- transcripts; no spot filtering (explicit exploratory baseline).
4. Library size normalization to 10,000 RNA counts per spot, `log1p` transform, select top 1,200 variable RNA genes (detected in ≥2% of spots).
5. PCA (20 randomized PCs) and **expression-only** K-means clustering on first 15 PCs, k=4 (not spatial clustering). **Spatial coordinates are not used to select clusters**; they are used to project clusters and evaluate neighbourhood structure.
6. Construct a tissue graph from the six nearest neighbours by image pixel distance, excluding distances >1.3× the median nearest neighbour distance and de-duplicating undirected edges.
7. Count-thinning sensitivity: Binomial sampling of counts at 50%, 25%, and 25th-percentile target depth per spot. Three independent stochastic runs per condition. Re-normalize, refit PCA/K-means using a fixed baseline HVG set. Compare spot cluster assignment to baseline using Adjusted Rand Index (ARI).
8. Neighbour-disagreement sensitivity: proportion of perturbed runs in which two connected spots receive different cluster assignments. This is a descriptive *boundary frequency*, not a calibrated posterior probability.
9. Spatial structure diagnostics: graph Moran's I of log RNA counts and observed fraction of neighbours assigned the same cluster. Two illustrative permutation nulls: globally shuffled labels and labels shuffled within 4×4 spatial tiles. **Neither null preserves the exact spatial autocorrelation function; their p-values are exploratory and are NOT calibrated evidence of a biological boundary.**
10. Additional depth nuisance sensitivity: regress log RNA depth from normalized individual features, refit PCA and k=4 clustering, and compare ARI. **This sensitivity test is deliberately aggressive and could remove true biology**; it is not an improved normalization or ground truth.
11. Parameter sensitivity: compare k=3 and k=5 against k=4 at fixed PCA representation.
12. Descriptive boundary-versus-depth comparison: AUROC of absolute *local* log count difference for identifying edges crossing an expression cluster boundary. The AUROC is descriptive (edges are spatially dependent), not evidence that technical depth caused the boundary.
13. Exploratory RNA maps for EPCAM, KRT19, and PTPRC as plausibility illustrations; no biomarker discovery or histological validation claimed.

## Main results

| Metric | P5 | P12 |
| --- | ---: | ---: |
| Spearman correlation PC1 vs log RNA counts | +0.868 | −0.951 |
| Mean ARI after half-depth thinning (3 runs) | 0.834 | 0.875 |
| Mean ARI after quarter-depth thinning (3 runs) | 0.705 | 0.701 |
| Mean ARI after 25th percentile equal-depth thinning (3 runs) | 0.721 | 0.531 |
| ARI vs depth-regressed alternative | 0.412 | 0.368 |
| Descriptive boundary depth-jump AUROC | 0.635 | 0.728 |
| Moran's I for log RNA counts on tissue graph | 0.833 | 0.858 |
| Same-cluster edge fraction (k=4) | 0.765 | 0.874 |

Strong depth dependence *after conventional library size normalization* is a useful hypothesis-generating finding. It **does not prove depth is a technical artifact**, because RNA content, cellular density, and spatial tissue morphology may produce genuine spatial gradients in RNA counts.

## Key limitations

- Only two distinct patients, one spot-based Visium section each: no valid cohort-level contrast or progression inference. Repeated spots, neighbouring edges and perturbation replicates are **not independent patient-level biological replicates**.
- No independent pathology segmentation labels or molecular ground truth available within these supplied files. Coloured expression clusters are **not** established tumor, stroma, normal lung, AIS, or MIA.
- Depth residualization may remove biologically meaningful variation. Thinning both removes reads and changes measurement precision, so lower ARI cannot be equated to a measured fraction of false clusters.
- `k=4` is selected for a demonstration; k=3/5 checks are sensitivity evaluations, not a proof of optimal state number.
- Permutation nulls are intentionally simple and do **not** maintain conditional spatial autocorrelation; a properly designed conditional spatial null must be validated with dedicated simulation before offering FDR/p-value claims.
- No independent replication beyond these two observational specimens, and no genomic/histological/radiomic multimodal integration yet.
- The methods (clustering, count thinning, stability comparisons, spatial graph statistics) are **not new methods**. Their combination is a feasibility exercise, not a novelty or journal-acceptance claim.

## Outputs

- `P5_spatial_analysis.png`, `P12_spatial_analysis.png`: counts, detected genes, expression groups, perturbation disagreement, PCA and depth correlation overlaid on original low-resolution histology
- `P5_gene_overlay.png`, `P12_gene_overlay.png`: spatial log-normalized RNA gene maps
- `cross_sample_thinning_stability.png`: depth-thinning sensitivity plot
- `sample_results.csv`: QC and graph statistics
- `perturbation_runs.csv`: all nine stochastic thinning runs per patient
- `additional_sensitivity_results.csv`: depth-residualized labels, parameter sensitivity and local-depth effect sizes
- `*_spot_level_metrics.csv.gz`: spot-level metrics, coordinates, original expression k=4 assignment, stability diagnostics
- `*_selected_HVGs.csv`: selected RNA features
- `run.log`, `extended.log`: original run logs
- `run_analysis.py`, `extended_checks.py`: Python reproduction scripts

## How to reproduce

```bash
python -m pip install numpy pandas scipy scikit-learn matplotlib pillow h5py

# Place the original four input files in /path/to/input
export OPENBLAS_NUM_THREADS=2
export OMP_NUM_THREADS=2
export MKL_NUM_THREADS=2
python run_analysis.py --input-dir /path/to/input --output-dir /path/to/output
# Extended checks use fixed paths /mnt/data and are intended for the original container setup.
python extended_checks.py
```

Seed is fixed at 26023. Floating-point and multithreading differences may cause small differences across Python/scikit-learn versions. The results have been executed locally in this response and are not hypothetical expected results.

## Hypothesis for further work

**Conditional reliability of spatial molecular-state boundaries:** can the evidence for each boundary survive changes in technical count depth, clustering parameters and plausible spatially correlated background structures, with calibrated uncertainty when depth is biologically confounded?

A publishable methodological contribution requires a formally defined estimand and null hypothesis, simulation under known biological and technical depth-generating mechanisms, spatially valid inference, comparison to existing spatial domain and uncertainty methods, and multi-patient blinded histology or orthogonal evidence. None of these novel inferential guarantees is claimed here.
