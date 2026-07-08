#include "Common.h"
#include "Dataset.h"
#include "BIRCH.h"
#include "Validation.h"
#include "Pipeline.h"
#include "KNNOutlier.h"

#include <map>

Score classScore(const std::vector<std::string>& truth, const std::vector<std::string>& predicted,
                 const std::string& name, bool excludeTruthNoise)
{
    size_t tp=0, fp=0, fn=0;
    for(size_t i=0; i<truth.size(); ++i)
    {
        if(excludeTruthNoise && truth[i] == "Noise") continue;
        const bool actual = truth[i] == name, prediction = predicted[i] == name;
        if(actual && prediction) ++tp; else if(!actual && prediction) ++fp; else if(actual) ++fn;
    }
    const double precision = tp + fp == 0 ? 0.0 : static_cast<double>(tp) / (tp + fp);
    const double recall = tp + fn == 0 ? 0.0 : static_cast<double>(tp) / (tp + fn);
    return {precision, recall, precision + recall == 0.0 ? 0.0 : 2.0 * precision * recall / (precision + recall)};
}

std::vector<std::string> alignClusterNamesForEvaluation(
    const std::vector<std::string>& truth, const std::vector<std::string>& predicted)
{
    if(truth.size()!=predicted.size())
        throw std::invalid_argument("Evaluation alignment size mismatch");
    std::map<std::string,std::map<std::string,size_t>> overlap;
    for(size_t i=0;i<truth.size();++i)
        if(predicted[i]!="Noise" && truth[i]!="Noise") ++overlap[predicted[i]][truth[i]];
    std::map<std::string,std::string> mapping;
    for(const auto& cluster : overlap)
    {
        size_t bestCount=0;
        for(const auto& candidate : cluster.second)
            if(candidate.second>bestCount)
            {
                bestCount=candidate.second;
                mapping[cluster.first]=candidate.first;
            }
    }
    std::vector<std::string> aligned=predicted;
    for(std::string& label : aligned)
        if(label!="Noise" && mapping.find(label)!=mapping.end()) label=mapping[label];
    return aligned;
}

void writeMetrics(std::ostream& output, const std::vector<std::string>& truth,
                  const std::vector<std::string>& predicted, const std::vector<size_t>& noiseIndices)
{
    std::vector<bool> predictedNoise(truth.size(), false);
    for(size_t index : noiseIndices) predictedNoise[index] = true;
    size_t tp=0, fp=0, fn=0, tn=0;
    for(size_t i=0; i<truth.size(); ++i)
    {
        const bool actual = truth[i] == "Noise";
        if(actual && predictedNoise[i]) ++tp; else if(!actual && predictedNoise[i]) ++fp;
        else if(actual) ++fn; else ++tn;
    }
    const double p = tp + fp == 0 ? 0.0 : static_cast<double>(tp)/(tp+fp);
    const double r = tp + fn == 0 ? 0.0 : static_cast<double>(tp)/(tp+fn);
    output << "TREE-BIRCH SMALL-CF NOISE DETECTION\nTP="<<tp<<" FP="<<fp<<" FN="<<fn<<" TN="<<tn
           <<"\nPrecision="<<p<<" Recall="<<r<<" F1="<<(p+r==0?0:2*p*r/(p+r))<<"\n\n";

    size_t emitterCorrect=0, emitterTotal=0, overallCorrect=0;
    for(size_t i=0; i<truth.size(); ++i)
    {
        if(truth[i] != "Noise") { ++emitterTotal; if(truth[i] == predicted[i]) ++emitterCorrect; }
        if(truth[i] == predicted[i]) ++overallCorrect;
    }
    output << "BIRCH EMITTER CLUSTERING\nAccuracy="
           << (emitterTotal == 0 ? 0.0 : static_cast<double>(emitterCorrect)/emitterTotal)
           << " ("<<emitterCorrect<<'/'<<emitterTotal<<")\n";
    std::vector<std::string> emitterNames;
    for(const std::string& label : truth)
        if(label != "Noise" && std::find(emitterNames.begin(), emitterNames.end(), label) == emitterNames.end())
            emitterNames.push_back(label);
    std::sort(emitterNames.begin(), emitterNames.end());
    for(const std::string& name : emitterNames)
    {
        const Score score = classScore(truth, predicted, name, true);
        output << name<<" Precision="<<score.precision<<" Recall="<<score.recall<<" F1="<<score.f1<<'\n';
    }
    output << "\nCOMBINED PIPELINE\nAccuracy=" << static_cast<double>(overallCorrect)/truth.size()
           << " ("<<overallCorrect<<'/'<<truth.size()<<")\n";
    emitterNames.push_back("Noise");
    for(const std::string& name : emitterNames)
    {
        const Score score = classScore(truth, predicted, name, false);
        output << name<<" Precision="<<score.precision<<" Recall="<<score.recall<<" F1="<<score.f1<<'\n';
    }

    // Include every discovered cluster in the matrix. This matters when
    // tree-BIRCH finds more geometric clusters than exist in ground truth.
    for(const std::string& name : predicted)
        if(std::find(emitterNames.begin(), emitterNames.end(), name) == emitterNames.end())
            emitterNames.push_back(name);

    output << "\nCONFUSION MATRIX (rows=Actual, columns=Predicted)\nActual/Predicted";
    for(const std::string& name : emitterNames) output << ',' << name;
    output << '\n';
    std::vector<std::vector<size_t>> matrix(emitterNames.size(),
                                             std::vector<size_t>(emitterNames.size(), 0));
    for(size_t i=0; i<truth.size(); ++i)
    {
        const auto actual = std::find(emitterNames.begin(), emitterNames.end(), truth[i]);
        const auto guess = std::find(emitterNames.begin(), emitterNames.end(), predicted[i]);
        if(actual != emitterNames.end() && guess != emitterNames.end())
            ++matrix[actual-emitterNames.begin()][guess-emitterNames.begin()];
    }
    for(size_t row=0; row<emitterNames.size(); ++row)
    {
        output << emitterNames[row];
        for(size_t column=0; column<emitterNames.size(); ++column)
            output << ',' << matrix[row][column];
        output << '\n';
    }

    // Adjusted Rand index from the contingency matrix. Unlike accuracy, ARI is
    // invariant to cluster-name permutations and corrects for chance agreement.
    const auto pairs = [](double n) { return n * (n - 1.0) / 2.0; };
    double sumCells=0.0, sumRows=0.0, sumColumns=0.0;
    for(size_t row=0; row<matrix.size(); ++row)
    {
        double rowTotal=0.0;
        for(size_t column=0; column<matrix.size(); ++column)
        { sumCells += pairs(matrix[row][column]); rowTotal += matrix[row][column]; }
        sumRows += pairs(rowTotal);
    }
    for(size_t column=0; column<matrix.size(); ++column)
    {
        double columnTotal=0.0;
        for(size_t row=0; row<matrix.size(); ++row) columnTotal += matrix[row][column];
        sumColumns += pairs(columnTotal);
    }
    const double totalPairs = pairs(truth.size());
    const double expected = totalPairs == 0.0 ? 0.0 : sumRows * sumColumns / totalPairs;
    const double maximum = 0.5 * (sumRows + sumColumns);
    const double ari = maximum == expected ? 0.0 : (sumCells - expected) / (maximum - expected);
    output << "Adjusted Rand Index=" << ari << '\n';
}

void saveResults(const std::string& filename, const Dataset& data, const std::vector<std::string>& predicted)
{
    std::ofstream output(filename);
    if(!output) throw std::runtime_error("Cannot create result CSV: " + filename);
    for(const std::string& name : data.featureNames) output << name << ',';
    output << "Ground_Truth,Predicted_Cluster\n";
    output << std::setprecision(12);
    for(size_t i=0; i<data.raw.size(); ++i)
    {
        for(double value : data.raw[i]) output << value << ',';
        output << data.truth[i] << ',' << predicted[i] << '\n';
    }
}

namespace
{
Point validityMean(const std::vector<Point>& points)
{
    Point result(points.front().size(), 0.0);
    for(const Point& point : points)
        for(size_t d=0; d<point.size(); ++d) result[d] += point[d];
    for(double& value : result) value /= points.size();
    return result;
}

double validityDistance(const Point& first, const Point& second)
{
    return std::sqrt(distanceSquared(first, second));
}

Point validityVariance(const std::vector<Point>& points, const Point& center)
{
    Point result(center.size(), 0.0);
    for(const Point& point : points)
        for(size_t d=0; d<point.size(); ++d)
        {
            const double difference = point[d] - center[d];
            result[d] += difference * difference;
        }
    for(double& value : result) value /= points.size();
    return result;
}

double validityNorm(const Point& values)
{
    double sum=0.0;
    for(double value : values) sum += value*value;
    return std::sqrt(sum);
}
}

void writeValidityIndices(std::ostream& output, const std::vector<Point>& normalizedPoints,
                          const std::vector<std::string>& predicted)
{
    if(normalizedPoints.size()!=predicted.size())
        throw std::invalid_argument("Validity-index point/label size mismatch");
    std::map<std::string,std::vector<Point>> clusters;
    for(size_t i=0;i<predicted.size();++i)
        if(predicted[i]!="Noise") clusters[predicted[i]].push_back(normalizedPoints[i]);

    output << "\nINTERNAL CLUSTER VALIDITY INDICES\n";
    output << "Predicted Noise rows excluded. Lower DB and SD values are preferred.\n";
    if(clusters.size()<2)
    {
        output << "Davies-Bouldin index=undefined (requires at least two clusters)\n";
        output << "SD index=undefined (requires at least two clusters)\n";
        return;
    }

    std::map<std::string,Point> centers;
    std::map<std::string,double> scatter;
    std::vector<Point> allPoints;
    for(const auto& cluster : clusters)
    {
        centers[cluster.first]=validityMean(cluster.second);
        double sum=0.0;
        for(const Point& point : cluster.second)
            sum += validityDistance(point,centers[cluster.first]);
        scatter[cluster.first]=sum/cluster.second.size();
        allPoints.insert(allPoints.end(),cluster.second.begin(),cluster.second.end());
    }

    // Davies-Bouldin: mean over clusters of the largest pairwise similarity.
    double db=0.0;
    for(const auto& first : clusters)
    {
        double worst=-1.0;
        for(const auto& second : clusters)
            if(first.first!=second.first)
            {
                const double centerDistance=validityDistance(centers[first.first],centers[second.first]);
                if(centerDistance>0.0)
                    worst=std::max(worst,(scatter[first.first]+scatter[second.first])/centerDistance);
            }
        if(worst<0.0) throw std::runtime_error("DB undefined for coincident centroids");
        db += worst;
    }
    db /= clusters.size();

    // SD index from the supplied implementation: SD = alpha*Scat + Dis,
    // with alpha=Dis(c_max) and c_max equal to this evaluated partition size.
    const Point overallCenter=validityMean(allPoints);
    const double overallVarianceNorm=validityNorm(validityVariance(allPoints,overallCenter));
    double scat=0.0;
    if(overallVarianceNorm>0.0)
    {
        for(const auto& cluster : clusters)
            scat += validityNorm(validityVariance(cluster.second,centers[cluster.first]))/
                    overallVarianceNorm;
        scat /= clusters.size();
    }
    std::vector<double> pairDistances;
    for(auto first=centers.begin();first!=centers.end();++first)
    {
        auto second=first;++second;
        for(;second!=centers.end();++second)
        {
            const double distance=validityDistance(first->second,second->second);
            if(distance>0.0) pairDistances.push_back(distance);
        }
    }
    if(pairDistances.empty()) throw std::runtime_error("SD undefined for coincident centroids");
    const double maximumDistance=*std::max_element(pairDistances.begin(),pairDistances.end());
    const double minimumDistance=*std::min_element(pairDistances.begin(),pairDistances.end());
    double inverseDistanceSum=0.0;
    for(const auto& first : centers)
    {
        double distanceSum=0.0;
        for(const auto& second : centers)
            if(first.first!=second.first)
                distanceSum += validityDistance(first.second,second.second);
        if(distanceSum==0.0) throw std::runtime_error("SD undefined for coincident centroids");
        inverseDistanceSum += 1.0/distanceSum;
    }
    const double dis=(maximumDistance/minimumDistance)*inverseDistanceSum;
    const double alpha=dis;
    const double sd=alpha*scat+dis;

    output << std::fixed << std::setprecision(6);
    output << "Points evaluated=" << allPoints.size() << '\n';
    output << "Clusters evaluated=" << clusters.size() << '\n';
    output << "Davies-Bouldin index=" << db << '\n';
    output << "SD Scat=" << scat << '\n';
    output << "SD Dis=" << dis << '\n';
    output << "SD Alpha=" << alpha << '\n';
    output << "SD index=" << sd << '\n';
}

void runSelfTests()
{
    const std::vector<Point> knnPoints{{0,0},{0.01,0},{0,0.01},{10,10}};
    const KNNOutlierResult knn=detectKNNOutliers(knnPoints,1,0.25);
    if(knn.count!=1 || !knn.outlier[3])
        throw std::logic_error("Distance-based kNN outlier test failed");

    CFTree tree(0.05, 3);
    for(int i=0; i<100; ++i) tree.insert({i < 50 ? i * 0.0001 : 1.0 + i * 0.0001, 0.0});
    const std::vector<CF> entries = tree.leafEntries();
    const size_t weight = std::accumulate(entries.begin(), entries.end(), size_t{0},
        [](size_t total, const CF& entry) { return total + entry.n; });
    if(weight != 100 || entries.size() < 2) throw std::logic_error("Recursive CF-tree test failed");

    CFTree mbdTree(0.05, 3, 0.10);
    for(int i=0; i<100; ++i)
        mbdTree.insert({i < 50 ? i*0.0001 : 1.0+i*0.0001, 0.0});
    const std::vector<CF> mbdEntries=mbdTree.leafEntries();
    const size_t mbdWeight=std::accumulate(mbdEntries.begin(),mbdEntries.end(),size_t{0},
        [](size_t total,const CF& entry){return total+entry.n;});
    if(mbdWeight!=100 || mbdEntries.size()<2)
        throw std::logic_error("MBD-BIRCH multiple-branch test failed");

    const std::vector<Point> pilot{{0,0},{0.01,0},{0,0.01},{1,1},{1.01,1},{1,1.01}};
    const ABirchEstimate estimate = estimateABirchThreshold(pilot, pilot.size(), 3);
    if(!(estimate.threshold > 0.0) || estimate.pilotClusters < 2)
        throw std::logic_error("A-BIRCH threshold-estimation test failed");
}
