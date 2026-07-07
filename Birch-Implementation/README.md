# A-BIRCH: BIRCH-only implementation

This implementation follows *Variations on the Clustering Algorithm BIRCH* by Lorbeer et al. It contains no LOF preprocessing and no final global k-means phase.

Pipeline:

1. Min-max normalize features.
   A strictly monotonic leading TOA/sequence column is retained in output but excluded from Euclidean clustering distance because it represents acquisition order rather than emitter identity.
2. Select a deterministic representative sample.
3. Estimate its cluster structure with Gap Statistic.
4. Calculate `Rmax`, `Dmin`, and `wrmax`.
5. Apply the paper's automatic threshold formula:

   `T = Dmin / (0.3 * wrmax + 4.5) + 0.8 * Rmax`

6. Build the recursive CF tree using `(N, LS, SS)`, radius-constrained absorption, farthest-pair splits and root growth.
7. Use leaf CFs directly as tree-BIRCH clusters. Small CFs are reported as noise; no separate outlier algorithm is used.
8. Validate only after clustering.

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

Use `python3 run_birch.py --help` for threshold override, pilot sample size, branching factor, and minimum cluster size.

## Direct build

```bash
g++ -std=c++14 -O2 -Wall -Wextra -Wpedantic \
  main.cpp Common.cpp Dataset.cpp BIRCH.cpp Validation.cpp Pipeline.cpp \
  -o birch
./birch --self-test
```
