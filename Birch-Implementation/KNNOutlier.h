#ifndef KNN_OUTLIER_H
#define KNN_OUTLIER_H

#include "Common.h"

struct KNNOutlierResult
{
    std::vector<double> kDistances;
    std::vector<bool> outlier;
    double cutoff = 0.0;
    size_t count = 0;
};

KNNOutlierResult detectKNNOutliers(const std::vector<Point>& points,
                                   size_t k, double outlierFraction);

#endif
