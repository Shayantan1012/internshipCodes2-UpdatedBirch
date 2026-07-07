#ifndef BIRCH_PIPELINE_H
#define BIRCH_PIPELINE_H

int runPipeline(const std::string& featureFile = "../Datasets/two_emitter_pdw_dataset.csv",
                const std::string& truthFile = "../Datasets/two_emitter_pdw_labeled.csv",
                const std::string& resultFile = "birch_results.csv",
                const std::string& reportFile = "birch_validation_report.txt",
                size_t branchingFactor = 32, size_t minimumClusterPoints = 5,
                size_t pilotPoints = 1000, size_t maximumPilotClusters = 10,
                double thresholdOverride = -1.0, double mbdSpread = 0.05,
                size_t knnK = 3, double knnOutlierFraction = 0.10);
void runSelfTests();

#endif
