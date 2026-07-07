#ifndef BIRCH_DATASET_H
#define BIRCH_DATASET_H

#include "Common.h"

struct Dataset
{
    std::vector<Point> raw;
    std::vector<Point> normalized;
    std::vector<std::string> truth;
    bool ignoredLeadingSequence = false;

    void loadFeatures(const std::string& filename);
    void normalize();
    void loadTruth(const std::string& filename);
};

#endif
