# SD Index for BIRCH Clusters

This folder contains a standalone C++17 implementation of the SD clustering
validity index from:

Halkidi, Vazirgiannis, and Batistakis, "Quality Scheme Assessment in the
Clustering Process."

The implementation evaluates the cluster labels already produced by BIRCH.
It does not run k-means or create a new partition.

## Formula

```text
SD(c) = alpha * Scat(c) + Dis(c)
alpha = Dis(c_max)
```

`Scat` measures average within-cluster scattering. `Dis` measures separation
between cluster centers. Lower SD is preferred when comparing partitions
evaluated on the same data and with the same `c_max` configuration.

The current workflow evaluates one fixed two-cluster BIRCH partition, so
`c_max = c = 2` and `alpha = Dis(2)`.

## Build on Ubuntu or WSL

From the project root:

```bash
cd SD


g++ -std=c++17 -O2 -Wall -Wextra -pedantic sd.cpp -o sd
```

## Run

```bash
./sd ../Birch-Implementation/birch_results.csv ../metrics/birch_sd_metrics.txt
```

The program:

1. Reads the original five feature columns and `Predicted_Cluster`.
2. Normalizes all features using the full BIRCH output range.
3. Excludes rows predicted as `Noise` by LOF.
4. Evaluates the remaining BIRCH emitter clusters.
5. Prints and saves `Scat`, `Dis`, `Alpha`, and `SD`.

Expected output:

```text
Points evaluated: 1811
BIRCH clusters: 2
Scat: 0.874785
Dis: 3.798511
Alpha: 3.798511
SD: 7.121390
```
