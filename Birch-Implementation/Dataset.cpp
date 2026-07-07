#include "Dataset.h"

void Dataset::loadFeatures(const std::string& filename)
{
    std::ifstream input(filename);
    if(!input) throw std::runtime_error("Cannot open feature CSV: " + filename);
    std::string line;
    if(!std::getline(input, line)) throw std::runtime_error("Feature CSV is empty");
    {
        std::stringstream headerParser(line);
        std::string name;
        while(std::getline(headerParser, name, ','))
        {
            if(!name.empty() && name.back()=='\r') name.pop_back();
            featureNames.push_back(name);
        }
    }
    size_t dimensions = 0;
    while(std::getline(input, line))
    {
        if(line.empty()) continue;
        std::stringstream parser(line);
        std::string value;
        Point point;
        while(std::getline(parser, value, ',')) point.push_back(std::stod(value));
        if(point.empty()) continue;
        if(dimensions == 0) dimensions = point.size();
        if(point.size() != dimensions) throw std::runtime_error("Inconsistent feature columns");
        for(double number : point)
            if(!std::isfinite(number)) throw std::runtime_error("Non-finite feature value");
        raw.push_back(std::move(point));
    }
    if(raw.empty()) throw std::runtime_error("Feature CSV contains no data rows");
    if(featureNames.size()!=dimensions)
        throw std::runtime_error("Feature header/data column-count mismatch");
}

void Dataset::normalize()
{
    const size_t dimensions = raw.front().size();
    Point minimum(dimensions, std::numeric_limits<double>::infinity());
    Point maximum(dimensions, -std::numeric_limits<double>::infinity());
    for(const Point& point : raw)
        for(size_t j=0; j<dimensions; ++j)
        {
            minimum[j] = std::min(minimum[j], point[j]);
            maximum[j] = std::max(maximum[j], point[j]);
        }
    normalized = raw;
    for(Point& point : normalized)
        for(size_t j=0; j<dimensions; ++j)
            point[j] = maximum[j] == minimum[j]
                ? 0.0
                : (point[j] - minimum[j]) / (maximum[j] - minimum[j]);

    // PDW files place TOA first. An almost strictly increasing acquisition-time
    // column describes observation order, not emitter identity; including it in
    // Euclidean BIRCH distance splits emitters into arbitrary time windows.
    if(dimensions > 1 && raw.size() > 2)
    {
        size_t increasing=0;
        for(size_t i=1;i<raw.size();++i) if(raw[i][0] > raw[i-1][0]) ++increasing;
        ignoredLeadingSequence = increasing >= static_cast<size_t>(0.98*(raw.size()-1));
        if(ignoredLeadingSequence) for(Point& point : normalized) point[0]=0.0;
    }
}

void Dataset::loadTruth(const std::string& filename)
{
    std::ifstream input(filename);
    if(!input) throw std::runtime_error("Cannot open labeled CSV: " + filename);
    std::string line;
    if(!std::getline(input, line)) throw std::runtime_error("Labeled CSV is empty");
    while(std::getline(input, line))
    {
        if(line.empty()) continue;
        std::stringstream parser(line);
        std::string value, last;
        while(std::getline(parser, value, ',')) last = value;
        if(!last.empty() && last.back() == '\r') last.pop_back();
        truth.push_back(last);
    }
    if(truth.size() != raw.size()) throw std::runtime_error("Feature/truth row-count mismatch");
}
