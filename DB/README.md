# Davies-Bouldin Index for BIRCH Clusters

This folder contains a standalone C++17 implementation of the cluster
separation measure from:

David L. Davies and Donald W. Bouldin, "A Cluster Separation Measure," IEEE
Transactions on Pattern Analysis and Machine Intelligence, 1979.

The implementation evaluates the emitter labels already produced by BIRCH.
It does not run k-means or generate a new partition.

## Formula

For each cluster `i`:

```text
S_i  = average distance of points in cluster i to centroid i
M_ij = distance between centroids i and j
R_ij = (S_i + S_j) / M_ij
R_i  = max(R_ij), j != i
DB   = average(R_i)
```

Lower DB values indicate clusters that are more compact and better separated.

## Build on Ubuntu or WSL

From the project root:

```bash
cd DB


g++ -std=c++17 -O2 -Wall -Wextra -pedantic db.cpp -o db
```

## Run

```bash
./db ../Birch-Implementation/birch_results.csv ../metrics/birch_db_metrics.txt
```

The program:

1. Reads the five original features and `Predicted_Cluster`.
2. Normalizes features using the full BIRCH output range.
3. Excludes rows predicted as `Noise` by LOF.
4. Evaluates the remaining BIRCH emitter clusters.
5. Saves cluster scatter, worst similarity, and the DB index.

