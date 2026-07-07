#include "Common.h"
#include "Dataset.h"
#include "KNNOutlier.h"
#include "BIRCH.h"
#include "Validation.h"
#include "Pipeline.h"

int runPipeline(const std::string& featureFile, const std::string& truthFile,
                const std::string& resultFile, const std::string& reportFile,
                size_t branchingFactor, size_t minimumClusterPoints,
                size_t pilotPoints, size_t maximumPilotClusters,
                double thresholdOverride, double mbdSpread,
                size_t knnK, double knnOutlierFraction)
{
    // ========================================================================
    // STAGE 0: Validate user-supplied configuration
    // ========================================================================
    // branchingFactor controls the maximum number of entries in a CF-tree node.
    // minimumClusterPoints controls when a leaf CF is large enough to be kept as
    // a detected cluster. pilotPoints and maximumPilotClusters configure the
    // representative sample used for automatic A-BIRCH threshold estimation.
    if(branchingFactor < 2 || minimumClusterPoints < 1 || pilotPoints < 2 ||
       maximumPilotClusters < 2)
        throw std::invalid_argument("Invalid A-BIRCH parameter");

    // ========================================================================
    // STAGE 1: Load and prepare the dataset
    // ========================================================================
    // featureFile contains only observations used by clustering. truthFile is
    // loaded for validation, but its labels are never supplied to A-BIRCH.
    // normalize() scales features and excludes a monotonic leading TOA/sequence
    // field from distance calculations when such a field is detected.
    Dataset data;
    data.loadFeatures(featureFile);
    data.normalize();
    data.loadTruth(truthFile);

    // ========================================================================
    // STAGE 2: Distance-based kNN outlier preprocessing
    // ========================================================================
    // Following Dang, Ngan and Liu, each point receives the distance to its
    // k-th nearest neighbor. The largest m distances are marked as outliers,
    // where m is configured here as a fraction of the dataset. Labels are not
    // used. Only the retained points proceed to A-BIRCH/MBD-BIRCH.
    const KNNOutlierResult knn=detectKNNOutliers(data.normalized,knnK,knnOutlierFraction);
    std::vector<Point> cleanPoints;
    std::vector<size_t> cleanIndices;
    for(size_t i=0;i<data.normalized.size();++i)
        if(!knn.outlier[i]){cleanPoints.push_back(data.normalized[i]);cleanIndices.push_back(i);}
    if(cleanPoints.size()<2) throw std::runtime_error("kNN preprocessing retained too few points");

    // ========================================================================
    // STAGE 3: Estimate the A-BIRCH threshold automatically
    // ========================================================================
    // The representative sample is analyzed with Gap Statistic to estimate its
    // natural cluster count. From that structure, A-BIRCH calculates:
    //   Rmax  = maximum estimated cluster radius
    //   Dmin  = minimum distance between estimated cluster centers
    //   wrmax = maximum cluster-size ratio
    // These quantities enter the threshold formula from Lorbeer et al.
    // Ground-truth labels do not participate in this calculation.
    const ABirchEstimate estimate = estimateABirchThreshold(
        cleanPoints, pilotPoints, maximumPilotClusters);

    // A non-negative override is useful for controlled experiments. The normal
    // execution path uses the automatically estimated paper threshold.
    const double threshold = thresholdOverride >= 0.0 ? thresholdOverride : estimate.threshold;

    // ========================================================================
    // STAGE 4: Construct the BIRCH CF tree
    // ========================================================================
    // Every normalized point is inserted incrementally. A leaf CF stores the
    // sufficient statistics (N, LS, SS). Points are absorbed only when the
    // resulting radius stays below threshold; overflowing nodes split.
    // This is tree-BIRCH only: there is no LOF and no global k-means phase.
    CFTree tree(threshold, branchingFactor, mbdSpread);
    for(const Point& point : cleanPoints) tree.insert(point);

    // Collect the final leaf-level micro-clusters produced by the CF tree.
    const std::vector<CF> leafEntries = tree.leafEntries();

    // Verify the most important CF-tree invariant: the sum of all leaf weights
    // must equal the number of inserted observations. Otherwise points were lost
    // or counted twice during recursive node splitting.
    const size_t treeWeight = std::accumulate(leafEntries.begin(), leafEntries.end(), size_t{0},
        [](size_t total, const CF& entry) { return total + entry.n; });
    if(treeWeight != cleanPoints.size()) throw std::logic_error("CF-tree lost observations");

    // ========================================================================
    // STAGE 5: Convert leaf CFs into final detected clusters
    // ========================================================================
    // The largest eligible leaf CFs are retained, up to the cluster count found
    // automatically by Gap Statistic. Small or surplus CFs receive label -1 and
    // are reported as BIRCH noise. No ground-truth class count is used here.
    size_t clusterCount = 0;
    const std::vector<int> cleanLabels = assignTreeBirchClusters(
        cleanPoints, leafEntries, minimumClusterPoints, estimate.pilotClusters, clusterCount);
    std::vector<int> labels(data.normalized.size(),-1);
    for(size_t i=0;i<cleanIndices.size();++i) labels[cleanIndices[i]]=cleanLabels[i];

    // Numeric cluster IDs are anonymous. emitterNames() gives stable readable
    // names using increasing mean feature values; this does not alter membership.
    const std::vector<std::string> anonymousNames = emitterNames(data.raw, labels, clusterCount);
    // Cluster IDs are arbitrary. Ground truth is used only now, after fitting,
    // to align display names for interpretable accuracy/confusion matrices.
    const std::vector<std::string> predicted =
        alignClusterNamesForEvaluation(data.truth, anonymousNames);

    // ========================================================================
    // STAGE 6: Collect cluster and noise populations
    // ========================================================================
    std::vector<size_t> noiseIndices;
    for(size_t i = 0; i < labels.size(); ++i)
        if(labels[i] < 0) noiseIndices.push_back(i);

    std::vector<size_t> clusterPopulations(clusterCount, 0);
    for(int label : labels)
        if(label >= 0 && static_cast<size_t>(label) < clusterCount)
            ++clusterPopulations[static_cast<size_t>(label)];

    // ========================================================================
    // STAGE 7: Save row-level predictions and create the validation report
    // ========================================================================
    // Ground truth is consulted only from this point onward. saveResults writes
    // one prediction per original observation. The report records threshold
    // diagnostics, detected clusters, populations, confusion matrix and metrics.
    saveResults(resultFile, data, predicted);
    std::ofstream report(reportFile);
    if(!report) throw std::runtime_error("Cannot create validation report");
    report << "A-BIRCH / TREE-BIRCH VALIDATION REPORT\n"
           << "=====================================\n\n"
           << "No LOF and no global k-means phase were used.\n\n"
           << "DISTANCE-BASED kNN OUTLIER PREPROCESSING\n"
           << "k=" << knnK << " Outlier fraction=" << knnOutlierFraction
           << "\nk-distance cutoff=" << knn.cutoff
           << "\nPreprocessing outliers removed=" << knn.count
           << "\nPoints retained for BIRCH=" << cleanPoints.size() << "\n\n"
           << "MBD-BIRCH multiple-branch spread s=" << mbdSpread << "\n"
           << "MBD-BIRCH enabled=" << (mbdSpread>0.0?"YES":"NO") << "\n\n"
           << "Monotonic leading sequence/TOA excluded from distance: "
           << (data.ignoredLeadingSequence?"YES":"NO") << "\n\n"
           << "AUTOMATIC THRESHOLD ESTIMATION\n"
           << "Pilot clusters (Gap Statistic)=" << estimate.pilotClusters
           << "\nRmax=" << estimate.maximumRadius
           << " Dmin=" << estimate.minimumDistance
           << " wrmax=" << estimate.maximumWeightRatio
           << "\nFormula threshold=" << estimate.threshold
           << " Applied threshold=" << threshold
           << "\nCondition 6*Rmax <= Dmin: " << (estimate.separationCondition?"PASS":"WARNING")
           << "\nWeight-ratio condition: " << (estimate.weightCondition?"PASS":"WARNING")
           << "\n\nTREE RESULTS\nBranching factor=" << branchingFactor
           << " Minimum cluster points=" << minimumClusterPoints
           << "\nLeaf CFs=" << leafEntries.size()
           << "\nAutomatically estimated clusters (Gap Statistic)=" << estimate.pilotClusters
           << "\nFinal detected BIRCH clusters=" << clusterCount
           << "\nRetained clusters=" << clusterCount
           << " Noise rows from small CFs=" << noiseIndices.size()
           << " CF weight=" << treeWeight << "\n\n"
           << "DETECTED CLUSTER POPULATIONS\n";
    for(size_t cluster = 0; cluster < clusterPopulations.size(); ++cluster)
        report << "Internal cluster " << cluster << "="
               << clusterPopulations[cluster] << " rows\n";
    report << "Noise=" << noiseIndices.size() << " rows\n"
           << "Total assigned=" << labels.size() << " rows\n\n";
    writeMetrics(report, data.truth, predicted, noiseIndices);
    writeValidityIndices(report, data.normalized, predicted);

    // Print the same metric summary to the terminal for immediate feedback.
    writeMetrics(std::cout, data.truth, predicted, noiseIndices);
    std::cout << "\nRows=" << data.raw.size() << " LeafCFs=" << leafEntries.size()
              << " Clusters=" << clusterCount << " Noise=" << noiseIndices.size()
              << " Threshold=" << threshold << "\nCreated " << resultFile
              << " and " << reportFile << '\n';
    return 0;
}
