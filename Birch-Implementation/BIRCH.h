#ifndef BIRCH_TREE_H
#define BIRCH_TREE_H

#include "Common.h"

struct CF
{
    size_t n = 0;
    Point linearSum;
    Point squareSum;

    explicit CF(size_t dimensions);

    void add(const Point& point);
    void merge(const CF& other);
    Point centroid() const;
    double radius() const;
};

struct Node
{
    bool leaf;
    std::vector<CF> entries;
    std::vector<std::unique_ptr<Node>> children;
    explicit Node(bool isLeaf);
};

struct ABirchEstimate
{
    double threshold = 0.0;
    double maximumRadius = 0.0;
    double minimumDistance = 0.0;
    double maximumWeightRatio = 1.0;
    size_t pilotClusters = 0;
    bool separationCondition = false;
    bool weightCondition = false;
};

class CFTree
{
    double threshold;
    double mbdSpread;
    size_t branchingFactor;
    std::unique_ptr<Node> root;

    static CF summary(const Node& node);
    static size_t closest(const std::vector<CF>& entries, const Point& point);
    std::unique_ptr<Node> split(Node& node);
    std::unique_ptr<Node> insertRecursive(Node& node, const Point& point);
    void collectMbdLeaves(Node& node, const Point& point, std::vector<Node*>& leaves) const;
    bool insertIntoLeaf(Node& node, Node* target, const Point& point,
                        std::unique_ptr<Node>& sibling);
    static void collect(const Node& node, std::vector<CF>& result);

public:
    CFTree(double radiusThreshold, size_t factor, double multipleBranchSpread = 0.0);
    void insert(const Point& point);
    std::vector<CF> leafEntries() const;
};

std::vector<std::string> emitterNames(const std::vector<Point>& raw, const std::vector<int>& labels,
                                      size_t clusterCount);
ABirchEstimate estimateABirchThreshold(const std::vector<Point>& points,
                                       size_t maximumPilotPoints = 1000,
                                       size_t maximumPilotClusters = 10);
std::vector<int> assignTreeBirchClusters(const std::vector<Point>& points,
                                         const std::vector<CF>& entries,
                                         size_t minimumClusterPoints,
                                         size_t maximumRetainedClusters,
                                         size_t& clusterCount);

#endif
