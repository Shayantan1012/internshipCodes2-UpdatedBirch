# BIRCH analysis report

Generated from the supplied CSV files with standardized features, a CF compression threshold of 0.65, branching/memory parameter 32, and silhouette-based final cluster selection. Ground-truth labels were not used for fitting or choosing `k`.

| Dataset | Rows | Ground-truth signal classes | CF micro-clusters | Discovered clusters | Silhouette | Purity | ARI |
|---|---:|---:|---:|---:|---:|---:|---:|
| s1_5d | 32,000 | 3 | 338 | **4** | 0.639 | 0.833 | 0.742 |
| s2_5d | 32,343 | 2 | 604 | **3** | 0.646 | 1.000 | 1.000 |
| s5_5d | 32,000 | 2 | 481 | **3** | 0.611 | 1.000 | 0.999 |
| s6_5d | 32,000 | 3 | 242 | **2** | 0.631 | 0.815 | 0.753 |
| two_emitter_pdw_labeled | 1,900 | 2 | 70 | **9** | 0.970 | 1.000 | 0.503 |

## What the clusters mean

- **s1_5d — 4 clusters.** One compact low-PW/high-angle group (14,728 rows), two high-PW/low-angle groups (4,631 and 10,359), and one smaller high-amplitude group (2,282). The model over-segments relative to the three signal emitter labels and has moderate label agreement.
- **s2_5d — 3 clusters.** Two signal signatures are recovered perfectly. The additional 3,144-row group is dominated by the known noise/outlier structure; this explains why three geometric clusters coexist with two signal classes.
- **s5_5d — 3 clusters.** Two signal emitters are recovered almost perfectly. A 2,841-row high-amplitude/intermediate-PW group separates from the two main emitter profiles and largely represents outlier/noise behavior.
- **s6_5d — 2 clusters.** The unsupervised geometry favors a high-PW/low-angle group and a low-PW/high-angle group. Three labeled emitters therefore merge into two broader regimes. This is the main under-clustering case and should be rerun with `--clusters 3` when emitter enumeration is the goal.
- **two_emitter_pdw_labeled — 9 clusters, but still 2 emitters.** Four large clusters split the two emitters by TOA region: Emitter 1 appears as groups of 572 and 442; Emitter 2 as 445 and 375. Five tiny clusters (10–16 points each) capture noise/outliers. Purity is 1.0, while ARI is only 0.503 because ARI correctly penalizes this over-segmentation. Thus nine is the number of compact geometric groups, not the number of physical emitters.

## Important interpretation rule

Silhouette selects compact geometric groups; it does not know what an “emitter” is. For physical emitter counting, inspect `cluster_profiles.csv`, noise populations, and ARI/confusion results together. To test a known emitter hypothesis, rerun with `--clusters K`; the default remains fully unsupervised and label-independent.
