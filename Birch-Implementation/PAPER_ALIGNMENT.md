# Alignment with “Variations on the Clustering Algorithm BIRCH”

The code implements the paper's A-BIRCH direction:

- Automatic threshold estimation from representative data.
- Gap Statistic for selecting the pilot cluster count.
- Estimation of maximum radius `Rmax`, minimum inter-cluster distance `Dmin`, and maximum neighboring weight ratio `wrmax`.
- Threshold formula `T = Dmin/(0.3*wrmax+4.5) + 0.8*Rmax`.
- Explicit warnings for both sufficient conditions given by the paper.
- Tree-BIRCH output without a required cluster-count parameter.
- No final global clustering phase.
- No LOF or other external outlier detector.
- MBD-BIRCH multiple-branch descent using Equation (11): all children within `s` of the nearest-child distance are searched, while each point is inserted only once.

The paper targets two-dimensional isotropic Gaussian clusters. The supplied radar datasets are five-dimensional and need not satisfy that assumption, so the generated condition warnings are important rather than cosmetic. Small or surplus CF entries beyond the Gap-Statistic pilot cluster estimate are labeled as noise; this is an engineering reporting rule, not LOF or global clustering.
