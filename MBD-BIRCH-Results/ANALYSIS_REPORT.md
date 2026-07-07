# A-BIRCH analysis report

Every dataset was fitted independently using A-BIRCH. Ground truth was used only after clustering for validation.

| Dataset | Rows | Threshold | Leaf CFs | Gap estimate | Final detected clusters | BIRCH noise | Signal accuracy | Overall accuracy | ARI | Paper conditions |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| s1_5d | 32000 | 0.492724 | 2 | 2 | **2** | 3200 | 0.8313 | 0.8315 | 0.7668 | separation: WARNING; weight: WARNING |
| s2_5d | 32343 | 0.415932 | 2 | 2 | **2** | 3235 | 0.9982 | 0.9731 | 0.9411 | separation: WARNING; weight: WARNING |
| s5_5d | 32000 | 0.334345 | 3 | 3 | **3** | 3200 | 0.9981 | 0.9533 | 0.8953 | separation: WARNING; weight: WARNING |
| s6_5d | 32000 | 0.488649 | 2 | 2 | **2** | 3200 | 0.8132 | 0.7497 | 0.5994 | separation: WARNING; weight: WARNING |
| two_emitter_pdw | 1900 | 0.114751 | 2 | 2 | **2** | 190 | 0.9500 | 0.9526 | 0.9013 | separation: PASS; weight: PASS |

## Graphs

![Signal accuracy](accuracy_comparison.svg)

![Adjusted Rand Index](ari_comparison.svg)

![Automatic thresholds](threshold_comparison.svg)

## Interpretation

- `s2_5d`, `s5_5d`, and the two-emitter dataset separate their signal emitters very accurately.
- `s1_5d` and `s6_5d` are under-clustered by pure tree-BIRCH: distinct labeled emitters merge geometrically.
- Noise detection is strongest for `s2_5d`; in the other datasets, the automatic radius absorbs noise into signal CFs.
- All paper-condition warnings are retained because the radar data is five-dimensional and not guaranteed to be isotropic Gaussian data.

Each dataset folder contains `birch_results.csv`, `birch_validation_report.txt`, `confusion_matrix.csv`, `confusion_matrix.svg`, and `cluster_sizes.svg`.
