# Distance-based kNN noise detection + A-BIRCH + MBD-BIRCH research implementation

The canonical BIRCH-only implementation combines automatic A-BIRCH threshold estimation with MBD-BIRCH multiple-branch descent. It is in [`Birch-Implementation`](Birch-Implementation/README.md), and the canonical input files are in `Datasets/`.

## Easy Ubuntu run

```bash
sudo apt update
sudo apt install build-essential python3
python3 run_birch.py
```

The script compiles the C++ implementation and analyzes all five supplied datasets. It converts the combined `s*_5d.csv` files internally, applies distance-based kNN outlier preprocessing, estimates an independent A-BIRCH threshold, and stores each result under `A-BIRCH-Results/<dataset>/`. No LOF or global k-means phase is used.

Only four clustering features are retained consistently across every dataset: frequency, pulse width, azimuth, and elevation. TOA, row-number fields, amplitude, labels, and trailing spreadsheet columns are excluded.

See all parameter overrides with:

```bash
python3 run_birch.py --help
```
