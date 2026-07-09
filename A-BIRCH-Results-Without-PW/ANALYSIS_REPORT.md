# A-BIRCH analysis report

Feature set used for clustering: **Frequency + Azimuth + Elevation**.

Every dataset was fitted independently using A-BIRCH. Ground truth was used only after clustering for validation.

| Dataset | Rows | Threshold | Leaf CFs | Gap estimate | Final detected clusters | BIRCH noise | Signal accuracy | Overall accuracy | ARI | Paper conditions |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| s1_5d | 32000 | 0.228573 | 10 | 4 | **4** | 9095 | 0.6600 | 0.6853 | 0.4834 | separation: WARNING; weight: WARNING |
| s2_5d | 32343 | 0.396641 | 2 | 2 | **2** | 3235 | 0.9805 | 0.9422 | 0.8704 | separation: WARNING; weight: WARNING |
| s5_5d | 32000 | 0.191486 | 6 | 6 | **6** | 3200 | 0.9720 | 0.9087 | 0.7769 | separation: WARNING; weight: WARNING |
| s6_5d | 32000 | 0.264643 | 6 | 4 | **4** | 5128 | 0.7872 | 0.7658 | 0.6365 | separation: WARNING; weight: WARNING |
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
