#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using Point = std::vector<double>;

struct Record {
    Point point;
    std::string predictedLabel;
};

class Dataset {
public:
    static Dataset fromBirchCsv(const std::string& path) {
        std::ifstream input(path.c_str());
        if (!input) {
            throw std::runtime_error("Could not open BIRCH result file: " + path);
        }

        std::string line;
        if (!std::getline(input, line)) {
            throw std::runtime_error("BIRCH result file is empty.");
        }

        const std::vector<std::string> header = split(line);
        const std::vector<std::string> required = {
            "TOA_ns", "Freq_MHz", "PW_ns", "Az_deg", "El_deg",
            "Predicted_Cluster"
        };
        std::map<std::string, std::size_t> columns;
        for (std::size_t i = 0; i < header.size(); ++i) {
            columns[trim(header[i])] = i;
        }
        for (const std::string& name : required) {
            if (columns.find(name) == columns.end()) {
                throw std::runtime_error("Missing required CSV column: " + name);
            }
        }

        std::vector<Record> records;
        std::size_t lineNumber = 1;
        while (std::getline(input, line)) {
            ++lineNumber;
            if (trim(line).empty()) {
                continue;
            }

            const std::vector<std::string> cells = split(line);
            if (cells.size() < header.size()) {
                throw std::runtime_error(
                    "Incomplete CSV row at line " + std::to_string(lineNumber));
            }

            Record record;
            for (std::size_t i = 0; i < 5; ++i) {
                record.point.push_back(parseNumber(
                    cells[columns[required[i]]], lineNumber));
            }
            record.predictedLabel = trim(cells[columns["Predicted_Cluster"]]);
            records.push_back(std::move(record));
        }

        if (records.empty()) {
            throw std::runtime_error("No data rows were loaded.");
        }
        return Dataset(std::move(records));
    }

    void normalize() {
        const std::size_t dimensions = records_.front().point.size();
        Point minimums = records_.front().point;
        Point maximums = records_.front().point;

        for (const Record& record : records_) {
            for (std::size_t d = 0; d < dimensions; ++d) {
                minimums[d] = std::min(minimums[d], record.point[d]);
                maximums[d] = std::max(maximums[d], record.point[d]);
            }
        }

        for (Record& record : records_) {
            for (std::size_t d = 0; d < dimensions; ++d) {
                const double range = maximums[d] - minimums[d];
                record.point[d] = range == 0.0
                    ? 0.0
                    : (record.point[d] - minimums[d]) / range;
            }
        }
    }

    std::map<std::string, std::vector<Point>> birchClusters() const {
        std::map<std::string, std::vector<Point>> clusters;
        for (const Record& record : records_) {
            if (record.predictedLabel != "Noise") {
                clusters[record.predictedLabel].push_back(record.point);
            }
        }
        return clusters;
    }

private:
    explicit Dataset(std::vector<Record> records)
        : records_(std::move(records)) {}

    static std::string trim(std::string value) {
        const std::string whitespace = " \t\r\n";
        const std::size_t first = value.find_first_not_of(whitespace);
        if (first == std::string::npos) {
            return "";
        }
        const std::size_t last = value.find_last_not_of(whitespace);
        return value.substr(first, last - first + 1);
    }

    static std::vector<std::string> split(const std::string& line) {
        std::vector<std::string> cells;
        std::stringstream stream(line);
        std::string cell;
        while (std::getline(stream, cell, ',')) {
            cells.push_back(cell);
        }
        return cells;
    }

    static double parseNumber(const std::string& text, std::size_t lineNumber) {
        try {
            std::size_t parsed = 0;
            const std::string cleaned = trim(text);
            const double value = std::stod(cleaned, &parsed);
            if (parsed != cleaned.size() || !std::isfinite(value)) {
                throw std::invalid_argument("not finite");
            }
            return value;
        } catch (const std::exception&) {
            throw std::runtime_error(
                "Invalid numeric value at CSV line " +
                std::to_string(lineNumber) + ": " + text);
        }
    }

    std::vector<Record> records_;
};

class VectorMath {
public:
    static Point mean(const std::vector<Point>& points) {
        Point result(points.front().size(), 0.0);
        for (const Point& point : points) {
            for (std::size_t d = 0; d < point.size(); ++d) {
                result[d] += point[d];
            }
        }
        for (double& value : result) {
            value /= static_cast<double>(points.size());
        }
        return result;
    }

    static Point variance(const std::vector<Point>& points, const Point& center) {
        Point result(center.size(), 0.0);
        for (const Point& point : points) {
            for (std::size_t d = 0; d < point.size(); ++d) {
                const double difference = point[d] - center[d];
                result[d] += difference * difference;
            }
        }
        for (double& value : result) {
            value /= static_cast<double>(points.size());
        }
        return result;
    }

    static double norm(const Point& values) {
        double sum = 0.0;
        for (double value : values) {
            sum += value * value;
        }
        return std::sqrt(sum);
    }

    static double distance(const Point& first, const Point& second) {
        double sum = 0.0;
        for (std::size_t d = 0; d < first.size(); ++d) {
            const double difference = first[d] - second[d];
            sum += difference * difference;
        }
        return std::sqrt(sum);
    }
};

struct SDScore {
    double scatter = 0.0;
    double separation = 0.0;
    double alpha = 0.0;
    double value = 0.0;
    std::size_t clusters = 0;
    std::size_t points = 0;
};

class SDIndex {
public:
    SDScore evaluate(
        const std::map<std::string, std::vector<Point>>& clusters) const {
        if (clusters.size() < 2) {
            throw std::invalid_argument("SD requires at least two clusters.");
        }

        std::vector<Point> allPoints;
        std::map<std::string, Point> centers;
        for (const auto& entry : clusters) {
            if (entry.second.empty()) {
                throw std::invalid_argument("SD is undefined for empty clusters.");
            }
            centers[entry.first] = VectorMath::mean(entry.second);
            allPoints.insert(
                allPoints.end(), entry.second.begin(), entry.second.end());
        }

        const Point overallCenter = VectorMath::mean(allPoints);
        const double overallVarianceNorm = VectorMath::norm(
            VectorMath::variance(allPoints, overallCenter));

        double scatter = 0.0;
        if (overallVarianceNorm > 0.0) {
            for (const auto& entry : clusters) {
                scatter += VectorMath::norm(VectorMath::variance(
                    entry.second, centers.at(entry.first))) /
                    overallVarianceNorm;
            }
            scatter /= static_cast<double>(clusters.size());
        }

        std::vector<double> pairDistances;
        for (auto first = centers.begin(); first != centers.end(); ++first) {
            auto second = first;
            ++second;
            for (; second != centers.end(); ++second) {
                const double distance = VectorMath::distance(
                    first->second, second->second);
                if (distance > 0.0) {
                    pairDistances.push_back(distance);
                }
            }
        }
        if (pairDistances.empty()) {
            throw std::runtime_error(
                "SD separation is undefined for coincident cluster centers.");
        }

        const double maximumDistance = *std::max_element(
            pairDistances.begin(), pairDistances.end());
        const double minimumDistance = *std::min_element(
            pairDistances.begin(), pairDistances.end());

        double inverseDistanceSum = 0.0;
        for (const auto& first : centers) {
            double distanceSum = 0.0;
            for (const auto& second : centers) {
                if (first.first != second.first) {
                    distanceSum += VectorMath::distance(
                        first.second, second.second);
                }
            }
            if (distanceSum == 0.0) {
                throw std::runtime_error(
                    "SD separation is undefined for coincident cluster centers.");
            }
            inverseDistanceSum += 1.0 / distanceSum;
        }

        SDScore score;
        score.scatter = scatter;
        score.separation =
            (maximumDistance / minimumDistance) * inverseDistanceSum;

        // The paper defines alpha = Dis(c_max). This program evaluates one
        // fixed BIRCH partition, so c_max is the current cluster count.
        score.alpha = score.separation;
        score.value = score.alpha * score.scatter + score.separation;
        score.clusters = clusters.size();
        score.points = allPoints.size();
        return score;
    }
};

class ReportWriter {
public:
    static void write(const std::string& path, const SDScore& score) {
        std::ofstream output(path.c_str());
        if (!output) {
            throw std::runtime_error("Could not create report file: " + path);
        }
        output << std::fixed << std::setprecision(6);
        output << "SD INDEX REPORT FOR BIRCH CLUSTERS\n";
        output << "==================================\n\n";
        output << "LOF-predicted Noise rows excluded.\n";
        output << "No k-means is run.\n\n";
        output << "Points evaluated: " << score.points << '\n';
        output << "BIRCH clusters: " << score.clusters << "\n\n";
        output << "Scat: " << score.scatter << '\n';
        output << "Dis: " << score.separation << '\n';
        output << "Alpha (Dis(c_max), c_max=" << score.clusters
               << "): " << score.alpha << '\n';
        output << "SD: " << score.value << '\n';
    }
};

int main(int argc, char* argv[]) {
    if (argc < 2 || argc > 3) {
        std::cerr << "Usage: " << argv[0]
                  << " <birch_results.csv> [sd_report.txt]\n";
        return 1;
    }

    try {
        Dataset dataset = Dataset::fromBirchCsv(argv[1]);
        dataset.normalize();
        const std::map<std::string, std::vector<Point>> clusters =
            dataset.birchClusters();
        const SDScore score = SDIndex().evaluate(clusters);
        const std::string reportPath = argc == 3 ? argv[2] : "sd_report.txt";

        std::cout << std::fixed << std::setprecision(6);
        std::cout << "Points evaluated: " << score.points << '\n';
        std::cout << "BIRCH clusters: " << score.clusters << '\n';
        std::cout << "Scat: " << score.scatter << '\n';
        std::cout << "Dis: " << score.separation << '\n';
        std::cout << "Alpha: " << score.alpha << '\n';
        std::cout << "SD: " << score.value << '\n';

        ReportWriter::write(reportPath, score);
        std::cout << "Report written to: " << reportPath << '\n';
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
