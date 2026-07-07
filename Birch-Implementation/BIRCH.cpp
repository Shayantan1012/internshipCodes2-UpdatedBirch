#include "BIRCH.h"

CF::CF(size_t dimensions) : linearSum(dimensions, 0.0), squareSum(dimensions, 0.0) {}

void CF::add(const Point& point)
{
    if(point.size() != linearSum.size()) throw std::invalid_argument("CF dimension mismatch");
    ++n;
    for(size_t i=0; i<point.size(); ++i)
    {
        linearSum[i] += point[i];
        squareSum[i] += point[i] * point[i];
    }
}

void CF::merge(const CF& other)
{
    if(other.linearSum.size() != linearSum.size()) throw std::invalid_argument("CF dimension mismatch");
    n += other.n;
    for(size_t i=0; i<linearSum.size(); ++i)
    {
        linearSum[i] += other.linearSum[i];
        squareSum[i] += other.squareSum[i];
    }
}

Point CF::centroid() const
{
    Point center(linearSum.size(), 0.0);
    if(n != 0)
        for(size_t i=0; i<center.size(); ++i) center[i] = linearSum[i] / n;
    return center;
}

double CF::radius() const
{
    if(n == 0) return 0.0;
    double squared = 0.0;
    for(size_t i=0; i<linearSum.size(); ++i)
    {
        const double mean = linearSum[i] / n;
        squared += std::max(0.0, squareSum[i] / n - mean * mean);
    }
    return std::sqrt(squared);
}

Node::Node(bool isLeaf) : leaf(isLeaf) {}

CFTree::CFTree(double radiusThreshold, size_t factor, double multipleBranchSpread)
    : threshold(radiusThreshold), mbdSpread(multipleBranchSpread), branchingFactor(factor),
      root(std::make_unique<Node>(true))
{
    if(threshold < 0.0) throw std::invalid_argument("BIRCH threshold cannot be negative");
    if(mbdSpread < 0.0) throw std::invalid_argument("MBD-BIRCH spread cannot be negative");
    if(branchingFactor < 2) throw std::invalid_argument("Branching factor must be at least 2");
}

CF CFTree::summary(const Node& node)
{
    if(node.entries.empty()) throw std::logic_error("Cannot summarize empty CF node");
    CF result(node.entries.front().linearSum.size());
    for(const CF& entry : node.entries) result.merge(entry);
    return result;
}

size_t CFTree::closest(const std::vector<CF>& entries, const Point& point)
{
    size_t best = 0;
    double bestDistance = std::numeric_limits<double>::infinity();
    for(size_t i=0; i<entries.size(); ++i)
    {
        const double current = distanceSquared(entries[i].centroid(), point);
        if(current < bestDistance)
        {
            bestDistance = current;
            best = i;
        }
    }
    return best;
}

std::unique_ptr<Node> CFTree::split(Node& node)
{
    const size_t count = node.entries.size();
    size_t seed1 = 0, seed2 = 1;
    double farthest = -1.0;
    for(size_t i=0; i<count; ++i)
        for(size_t j=i+1; j<count; ++j)
        {
            const double current = distanceSquared(node.entries[i].centroid(), node.entries[j].centroid());
            if(current > farthest)
            {
                farthest = current;
                seed1 = i;
                seed2 = j;
            }
        }

    std::vector<CF> oldEntries = std::move(node.entries);
    std::vector<std::unique_ptr<Node>> oldChildren = std::move(node.children);
    auto sibling = std::make_unique<Node>(node.leaf);
    const Point center1 = oldEntries[seed1].centroid();
    const Point center2 = oldEntries[seed2].centroid();
    for(size_t i=0; i<count; ++i)
    {
        Node* destination = distanceSquared(oldEntries[i].centroid(), center1)
            <= distanceSquared(oldEntries[i].centroid(), center2) ? &node : sibling.get();
        destination->entries.push_back(std::move(oldEntries[i]));
        if(!node.leaf) destination->children.push_back(std::move(oldChildren[i]));
    }
    return sibling;
}

std::unique_ptr<Node> CFTree::insertRecursive(Node& node, const Point& point)
{
    if(node.leaf)
    {
        if(node.entries.empty())
        {
            CF entry(point.size());
            entry.add(point);
            node.entries.push_back(std::move(entry));
        }
        else
        {
            const size_t index = closest(node.entries, point);
            CF candidate = node.entries[index];
            candidate.add(point);
            if(candidate.radius() <= threshold) node.entries[index] = std::move(candidate);
            else
            {
                CF entry(point.size());
                entry.add(point);
                node.entries.push_back(std::move(entry));
            }
        }
    }
    else
    {
        const size_t index = closest(node.entries, point);
        std::unique_ptr<Node> sibling = insertRecursive(*node.children[index], point);
        node.entries[index] = summary(*node.children[index]);
        if(sibling)
        {
            node.entries.push_back(summary(*sibling));
            node.children.push_back(std::move(sibling));
        }
    }
    return node.entries.size() > branchingFactor ? split(node) : nullptr;
}

// Implements Equation (11) from the paper. At each internal node, MBD-BIRCH
// explores the nearest child and every other child whose distance is less than
// s farther away than the nearest child. The point is still inserted only once.
void CFTree::collectMbdLeaves(Node& node, const Point& point,
                              std::vector<Node*>& leaves) const
{
    if(node.leaf)
    {
        leaves.push_back(&node);
        return;
    }

    std::vector<double> distances(node.entries.size());
    double nearestDistance = std::numeric_limits<double>::infinity();
    size_t nearestIndex = 0;
    for(size_t i=0; i<node.entries.size(); ++i)
    {
        distances[i] = std::sqrt(distanceSquared(node.entries[i].centroid(), point));
        if(distances[i] < nearestDistance)
        {
            nearestDistance = distances[i];
            nearestIndex = i;
        }
    }
    for(size_t i=0; i<node.children.size(); ++i)
        if(i == nearestIndex || distances[i] - nearestDistance < mbdSpread)
            collectMbdLeaves(*node.children[i], point, leaves);
}

// Inserts into a leaf selected during the MBD search and updates/splits every
// node on that one root-to-leaf path. The boolean indicates that target was
// found below this node; sibling carries an overflow split back to the parent.
bool CFTree::insertIntoLeaf(Node& node, Node* target, const Point& point,
                            std::unique_ptr<Node>& sibling)
{
    if(node.leaf)
    {
        if(&node != target) return false;
        if(node.entries.empty())
        {
            CF entry(point.size());
            entry.add(point);
            node.entries.push_back(std::move(entry));
        }
        else
        {
            const size_t index = closest(node.entries, point);
            CF candidate = node.entries[index];
            candidate.add(point);
            if(candidate.radius() <= threshold) node.entries[index] = std::move(candidate);
            else
            {
                CF entry(point.size());
                entry.add(point);
                node.entries.push_back(std::move(entry));
            }
        }
        if(node.entries.size() > branchingFactor) sibling = split(node);
        return true;
    }

    for(size_t i=0; i<node.children.size(); ++i)
    {
        std::unique_ptr<Node> childSibling;
        if(!insertIntoLeaf(*node.children[i], target, point, childSibling)) continue;
        node.entries[i] = summary(*node.children[i]);
        if(childSibling)
        {
            node.entries.push_back(summary(*childSibling));
            node.children.push_back(std::move(childSibling));
        }
        if(node.entries.size() > branchingFactor) sibling = split(node);
        return true;
    }
    return false;
}

void CFTree::collect(const Node& node, std::vector<CF>& result)
{
    if(node.leaf)
    {
        result.insert(result.end(), node.entries.begin(), node.entries.end());
        return;
    }
    for(const auto& child : node.children) collect(*child, result);
}
 
void CFTree::insert(const Point& point)
{
    std::unique_ptr<Node> sibling;
    if(mbdSpread == 0.0 || root->leaf)
        sibling = insertRecursive(*root, point);
    else
    {
        std::vector<Node*> candidateLeaves;
        collectMbdLeaves(*root, point, candidateLeaves);
        if(candidateLeaves.empty()) throw std::logic_error("MBD-BIRCH found no candidate leaf");

        Node* bestLeaf = candidateLeaves.front();
        double bestDistance = bestLeaf->entries.empty()
            ? 0.0 : distanceSquared(bestLeaf->entries[closest(bestLeaf->entries, point)].centroid(), point);
        for(Node* leaf : candidateLeaves)
        {
            const double current = leaf->entries.empty()
                ? 0.0 : distanceSquared(leaf->entries[closest(leaf->entries, point)].centroid(), point);
            if(current < bestDistance) { bestDistance = current; bestLeaf = leaf; }
        }
        if(!insertIntoLeaf(*root, bestLeaf, point, sibling))
            throw std::logic_error("MBD-BIRCH candidate leaf disappeared");
    }
    if(sibling)
    {
        auto newRoot = std::make_unique<Node>(false);
        newRoot->entries.push_back(summary(*root));
        newRoot->children.push_back(std::move(root));
        newRoot->entries.push_back(summary(*sibling));
        newRoot->children.push_back(std::move(sibling));
        root = std::move(newRoot);
    }
}
 
std::vector<CF> CFTree::leafEntries() const
{
    std::vector<CF> result;
    collect(*root, result);
    return result;
}

std::vector<std::string> emitterNames(const std::vector<Point>& raw, const std::vector<int>& labels,
                                      size_t clusterCount)
{
    std::vector<double> frequencySum(clusterCount, 0.0);
    std::vector<size_t> count(clusterCount, 0);
    for(size_t i=0; i<labels.size(); ++i)
        if(labels[i] >= 0 && static_cast<size_t>(labels[i]) < clusterCount)
        {
            // Prepared datasets place frequency in the first retained column.
            frequencySum[labels[i]] += raw[i][0];
            ++count[labels[i]];
        }
    std::vector<size_t> order(clusterCount);
    std::iota(order.begin(), order.end(), size_t{0});
    for(size_t cluster=0; cluster<clusterCount; ++cluster)
        if(count[cluster] == 0) throw std::runtime_error("A global emitter cluster is empty");
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        return frequencySum[a] / count[a] < frequencySum[b] / count[b];
    });
    std::vector<std::string> clusterName(clusterCount);
    for(size_t rank=0; rank<clusterCount; ++rank)
        clusterName[order[rank]] = "Emitter_" + std::to_string(rank + 1);
    std::vector<std::string> names;
    for(int label : labels)
        names.push_back(label < 0 ? "Noise" : clusterName.at(static_cast<size_t>(label)));
    return names;
}

namespace
{
struct PilotResult { std::vector<int> labels; std::vector<Point> centers; double dispersion=0.0; };

PilotResult pilotKMeans(const std::vector<Point>& points, size_t k)
{
    std::vector<Point> centers{points.front()};
    while(centers.size() < k)
    {
        size_t farthest=0; double farthestDistance=-1.0;
        for(size_t i=0; i<points.size(); ++i)
        {
            double nearest=std::numeric_limits<double>::infinity();
            for(const Point& center : centers) nearest=std::min(nearest, distanceSquared(points[i], center));
            if(nearest > farthestDistance) { farthestDistance=nearest; farthest=i; }
        }
        centers.push_back(points[farthest]);
    }
    std::vector<int> labels(points.size(), -1);
    for(size_t iteration=0; iteration<200; ++iteration)
    {
        bool changed=false;
        for(size_t i=0; i<points.size(); ++i)
        {
            int best=0;
            double bestDistance=distanceSquared(points[i], centers.front());
            for(size_t c=1; c<k; ++c) { const double d=distanceSquared(points[i], centers[c]);
            if(d<bestDistance){bestDistance=d;best=static_cast<int>(c);} }
            if(labels[i]!=best){labels[i]=best;changed=true;}
        }
        std::vector<Point> updated(k, Point(points.front().size(),0.0)); std::vector<size_t> counts(k,0);
        for(size_t i=0;i<points.size();++i){++counts[labels[i]];for(size_t j=0;j<points[i].size();++j)updated[labels[i]][j]+=points[i][j];}
        for(size_t c=0;c<k;++c)if(counts[c])for(double&v:updated[c])v/=counts[c];else updated[c]=centers[c];
        centers.swap(updated); if(!changed)break;
    }
    double dispersion=0.0;for(size_t i=0;i<points.size();++i)dispersion+=distanceSquared(points[i],centers[labels[i]]);
    return {labels,centers,dispersion};
}
}

ABirchEstimate estimateABirchThreshold(const std::vector<Point>& points,
                                       size_t maximumPilotPoints,
                                       size_t maximumPilotClusters)
{
    if(points.size() < 2)
        throw std::invalid_argument("A-BIRCH threshold estimation needs at least two points");

    const size_t sampleSize = std::min(maximumPilotPoints, points.size());
    std::vector<Point> sample;
    sample.reserve(sampleSize);
    for(size_t i = 0; i < sampleSize; ++i)
        sample.push_back(points[i * points.size() / sampleSize]);

    const size_t maxK = std::min({maximumPilotClusters, sample.size() - 1, size_t{12}});
    std::vector<double> gap(maxK + 1), standard(maxK + 1);
    std::vector<PilotResult> observed(maxK + 1);              
    Point low(sample.front().size(), 1.0), high(sample.front().size(), 0.0);
    for(const Point& point : sample)
        for(size_t j = 0; j < point.size(); ++j)
        {
            low[j] = std::min(low[j], point[j]);
            high[j] = std::max(high[j], point[j]);
        }

    constexpr size_t references = 5;
    constexpr unsigned gapRandomSeed = 42;
    for(size_t k = 1; k <= maxK; ++k)
    {
        observed[k] = pilotKMeans(sample, k);
        std::vector<double> logs;
        for(size_t b = 0; b < references; ++b)
        {
            std::vector<Point> reference(sampleSize, Point(low.size()));
            std::seed_seq seed{
                gapRandomSeed,
                static_cast<unsigned>(k),
                static_cast<unsigned>(b)
            };
            std::mt19937 generator(seed);
            std::uniform_real_distribution<double> uniform(0.0, 1.0);
            for(Point& point : reference)
                for(size_t j = 0; j < point.size(); ++j)
                {
                    const double randomValue = uniform(generator);
                    point[j] = low[j] + randomValue * (high[j] - low[j]);
                }
            logs.push_back(std::log(std::max(1e-15, pilotKMeans(reference, k).dispersion)));
        }
        const double mean = std::accumulate(logs.begin(), logs.end(), 0.0) / logs.size();
        double variance = 0.0;
        for(double value : logs) variance += (value - mean) * (value - mean);
        variance /= logs.size();
        gap[k] = mean - std::log(std::max(1e-15, observed[k].dispersion));
        standard[k] = std::sqrt(variance) * std::sqrt(1.0 + 1.0 / references);
    }

    size_t chosen = maxK;
    for(size_t k = 1; k < maxK; ++k)
        if(gap[k] >= gap[k + 1] - standard[k + 1]) { chosen = k; break; }
    if(chosen < 2) chosen = 2;

    const PilotResult& model = observed[chosen];
    std::vector<size_t> counts(chosen);
    std::vector<double> radiusSquared(chosen);
    for(size_t i = 0; i < sample.size(); ++i)
    {
        ++counts[model.labels[i]];
        radiusSquared[model.labels[i]] += distanceSquared(sample[i], model.centers[model.labels[i]]);
    }

    ABirchEstimate estimate;
    estimate.pilotClusters = chosen;
    for(size_t cluster = 0; cluster < chosen; ++cluster)
        estimate.maximumRadius = std::max(estimate.maximumRadius,
            std::sqrt(radiusSquared[cluster] / std::max<size_t>(1, counts[cluster])));
    estimate.minimumDistance = std::numeric_limits<double>::infinity();
    estimate.maximumWeightRatio = 1.0;
    for(size_t a = 0; a < chosen; ++a)
        for(size_t b = a + 1; b < chosen; ++b)
        {
            estimate.minimumDistance = std::min(estimate.minimumDistance,
                                                 distance(model.centers[a], model.centers[b]));
            const double ratio = double(std::max(counts[a], counts[b])) /
                                 std::max<size_t>(1, std::min(counts[a], counts[b]));
            estimate.maximumWeightRatio = std::max(estimate.maximumWeightRatio, ratio);
        }
    estimate.threshold = estimate.minimumDistance / (0.3 * estimate.maximumWeightRatio + 4.5)
                         + 0.8 * estimate.maximumRadius;
    estimate.separationCondition = 6.0 * estimate.maximumRadius <= estimate.minimumDistance;
    estimate.weightCondition = estimate.maximumWeightRatio <=
        8.8 * estimate.minimumDistance / std::max(1e-15, estimate.maximumRadius) - 42.1;
    return estimate;
}

std::vector<int> assignTreeBirchClusters(const std::vector<Point>& points,
                                         const std::vector<CF>& entries,
                                         size_t minimumClusterPoints,
                                         size_t maximumRetainedClusters,
                                         size_t& clusterCount)
{
    std::vector<size_t> order(entries.size());
    std::iota(order.begin(), order.end(), size_t{0});
    std::sort(order.begin(), order.end(),
              [&](size_t a, size_t b) { return entries[a].n > entries[b].n; });

    std::vector<int> entryLabel(entries.size(), -1);
    clusterCount = 0;
    for(size_t index : order)
        if(entries[index].n >= minimumClusterPoints && clusterCount < maximumRetainedClusters)
            entryLabel[index] = static_cast<int>(clusterCount++);

    if(clusterCount == 0)
    {
        size_t largest = 0;
        for(size_t i = 1; i < entries.size(); ++i)
            if(entries[i].n > entries[largest].n) largest = i;
        entryLabel[largest] = 0;
        clusterCount = 1;
    }

    std::vector<int> labels;
    labels.reserve(points.size());
    for(const Point& point : points)
    {
        size_t best = 0;
        double bestDistance = distanceSquared(point, entries[0].centroid());
        for(size_t i = 1; i < entries.size(); ++i)
        {
            const double current = distanceSquared(point, entries[i].centroid());
            if(current < bestDistance) { bestDistance = current; best = i; }
        }
        labels.push_back(entryLabel[best]);
    }
    return labels;
}
