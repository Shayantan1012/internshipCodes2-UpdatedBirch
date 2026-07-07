# A-BIRCH: BIRCH-only implementation

This implementation combines the distance-based kNN outlier method of Dang, Ngan and Liu with A-BIRCH threshold estimation and MBD-BIRCH multiple-branch descent. It contains no LOF preprocessing and no final global k-means phase.

Pipeline:

1. Min-max normalize features.
   A strictly monotonic leading TOA/sequence column is retained in output but excluded from Euclidean clustering distance because it represents acquisition order rather than emitter identity.
2. Calculate every point's distance to its k-th nearest neighbor and remove the configured top fraction of k-distances as outliers before BIRCH. Defaults follow the paper's `k=3`; the cutoff fraction is configurable because the paper selects `m`/the histogram threshold manually.
3. Select a deterministic representative sample.
4. Estimate its cluster structure with Gap Statistic.
5. Calculate `Rmax`, `Dmin`, and `wrmax`.
6. Apply the paper's automatic threshold formula:

   `T = Dmin / (0.3 * wrmax + 4.5) + 0.8 * Rmax`

7. Build the recursive CF tree using `(N, LS, SS)`, radius-constrained absorption, farthest-pair splits and root growth.
   At every internal node, MBD-BIRCH also explores children satisfying `d(child,p)-d(nearest,p) < s`, then inserts the point once into the best candidate leaf.
8. Use leaf CFs directly as tree-BIRCH clusters. Small CFs are also reported as noise.
9. Validate only after clustering.

The report warns when the paper's sufficient applicability conditions fail:

- `6 * Rmax <= Dmin`
- `wrmax <= 8.8 * Dmin / Rmax - 42.1`

## Easy Ubuntu run

From the repository root:

```bash
sudo apt update
sudo apt install build-essential python3
python3 run_birch.py
```

Use `python3 run_birch.py --help` for threshold override, MBD spread `s`, pilot sample size, branching factor, and minimum cluster size. Set `--mbd-spread 0` to recover ordinary tree-BIRCH descent.

## Direct build

```bash
g++ -std=c++14 -O2 -Wall -Wextra -Wpedantic \
  main.cpp Common.cpp Dataset.cpp BIRCH.cpp Validation.cpp Pipeline.cpp \
  -o birch
./birch --self-test
```
