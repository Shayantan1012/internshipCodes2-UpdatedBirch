#ifndef BIRCH_VALIDATION_H
#define BIRCH_VALIDATION_H

#include "Common.h"
#include "Dataset.h"

struct Score
{
    double precision;
    double recall;
    double f1;
};

Score classScore(const std::vector<std::string>& truth, const std::vector<std::string>& predicted,
                 const std::string& name, bool excludeTruthNoise);
void writeMetrics(std::ostream& output, const std::vector<std::string>& truth,
                  const std::vector<std::string>& predicted, const std::vector<size_t>& noiseIndices);
void saveResults(const std::string& filename, const Dataset& data, const std::vector<std::string>& predicted);
void writeValidityIndices(std::ostream& output, const std::vector<Point>& normalizedPoints,
                          const std::vector<std::string>& predicted);

#endif
