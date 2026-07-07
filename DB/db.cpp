#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
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
        const std::vector<std::string> features = {
            "TOA_ns", "Freq_MHz", "PW_ns", "Az_deg", "El_deg"
        };
        std::map<std::string, std::size_t> columns;
        for (std::size_t i = 0; i < header.size(); ++i) {
            columns[trim(header[i])] = i;
        }
        for (const std::string& name : features) {
            requireColumn(columns, name);
        }
        requireColumn(columns, "Predicted_Cluster");

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
            for (const std::string& feature : features) {
                record.point.push_back(parseNumber(
                    cells[columns.at(feature)], lineNumber));
            }
            record.predictedLabel = trim(
                cells[columns.at("Predicted_Cluster")]);
            records.push_back(std::move(record));
        }
        if (records.empty()) {
            throw std::runtime_error("No data rows were loaded.");
        }
        return Dataset(std::move(records));
    }

    void normalize() {
        Point minimums = records_.front().point;
        Point maximums = records_.front().point;
        for (const Record& record : records_) {
            for (std::size_t d = 0; d < record.point.size(); ++d) {
                minimums[d] = std::min(minimums[d], record.point[d]);
                maximums[d] = std::max(maximums[d], record.point[d]);
            }
        }
        for (Record& record : records_) {
            for (std::size_t d = 0; d < record.point.size(); ++d) {
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

    static void requireColumn(
        const std::map<std::string, std::size_t>& columns,
        const std::string& name) {
        if (columns.find(name) == columns.end()) {
            throw std::runtime_error("Missing required CSV column: " + name);
        }
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

class Distance {
public:
    static double euclidean(const Point& first, const Point& second) {
        double sum = 0.0;
        for (std::size_t d = 0; d < first.size(); ++d) {
            const double difference = first[d] - second[d];
            sum += difference * difference;
        }
        return std::sqrt(sum);
    }
};

struct ClusterDetails {
    Point centroid;
    double scatter = 0.0;
    double worstSimilarity = 0.0;
    std::string mostSimilarCluster;
    std::size_t size = 0;
};

struct DBScore {
    double value = 0.0;
    std::size_t points = 0;
    std::map<std::string, ClusterDetails> details;
};

class DaviesBouldinIndex {
public:
    DBScore evaluate(
        const std::map<std::string, std::vector<Point>>& clusters) const {
        if (clusters.size() < 2) {
            throw std::invalid_argument(
                "Davies-Bouldin requires at least two clusters.");
        }

        DBScore score;
        for (const auto& entry : clusters) {
            if (entry.second.empty()) {
                throw std::invalid_argument(
                    "Davies-Bouldin is undefined for empty clusters.");
            }
            ClusterDetails details;
            details.centroid = mean(entry.second);
            details.size = entry.second.size();
            details.scatter = averageDistance(
                entry.second, details.centroid);
            score.points += entry.second.size();
            score.details[entry.first] = details;
        }

        for (auto& first : score.details) {
            double maximumSimilarity = -1.0;
            std::string mostSimilar;
            for (const auto& second : score.details) {
                if (first.first == second.first) {
                    continue;
                }
                const double centroidDistance = Distance::euclidean(
                    first.second.centroid, second.second.centroid);
                if (centroidDistance == 0.0) {
                    throw std::runtime_error(
                        "Davies-Bouldin is undefined for coincident centroids.");
                }
                const double similarity =
                    (first.second.scatter + second.second.scatter) /
                    centroidDistance;
                if (similarity > maximumSimilarity) {
                    maximumSimilarity = similarity;
                    mostSimilar = second.first;
                }
            }
            first.second.worstSimilarity = maximumSimilarity;
            first.second.mostSimilarCluster = mostSimilar;
            score.value += maximumSimilarity;
        }
        score.value /= static_cast<double>(clusters.size());
        return score;
    }

private:
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

    static double averageDistance(
        const std::vector<Point>& points, const Point& centroid) {
        double sum = 0.0;
        for (const Point& point : points) {
            sum += Distance::euclidean(point, centroid);
        }
        return sum / static_cast<double>(points.size());
    }
};

class ReportWriter {
public:
    static void write(const std::string& path, const DBScore& score) {
        std::ofstream output(path.c_str());
        if (!output) {
            throw std::runtime_error("Could not create report file: " + path);
        }
        output << std::fixed << std::setprecision(6);
        output << "DAVIES-BOULDIN INDEX REPORT FOR BIRCH CLUSTERS\n";
        output << "================================================\n\n";
        output << "LOF-predicted Noise rows excluded.\n";
        output << "No k-means is run.\n";
        output << "Euclidean centroid distance and mean Euclidean "
                  "within-cluster scatter are used.\n\n";
        output << "Points evaluated: " << score.points << '\n';
        output << "BIRCH clusters: " << score.details.size() << "\n\n";
        for (const auto& entry : score.details) {
            output << "Cluster " << entry.first << ":\n";
            output << "  Size: " << entry.second.size << '\n';
            output << "  Scatter: " << entry.second.scatter << '\n';
            output << "  Most similar cluster: "
                   << entry.second.mostSimilarCluster << '\n';
            output << "  Worst similarity R_i: "
                   << entry.second.worstSimilarity << '\n';
        }
        output << "\nDavies-Bouldin index: " << score.value << '\n';
        output << "Lower values indicate better separated, more compact clusters.\n";
    }
};

int main(int argc, char* argv[]) {
    if (argc < 2 || argc > 3) {
        std::cerr << "Usage: " << argv[0]
                  << " <birch_results.csv> [db_report.txt]\n";
        return 1;
    }

    try {
        Dataset dataset = Dataset::fromBirchCsv(argv[1]);
        dataset.normalize();
        const std::map<std::string, std::vector<Point>> clusters =
            dataset.birchClusters();
        const DBScore score = DaviesBouldinIndex().evaluate(clusters);
        const std::string reportPath = argc == 3 ? argv[2] : "db_report.txt";

        std::cout << std::fixed << std::setprecision(6);
        std::cout << "Points evaluated: " << score.points << '\n';
        std::cout << "BIRCH clusters: " << score.details.size() << '\n';
        for (const auto& entry : score.details) {
            std::cout << entry.first << " scatter: "
                      << entry.second.scatter << '\n';
            std::cout << entry.first << " worst similarity: "
                      << entry.second.worstSimilarity << '\n';
        }
        std::cout << "Davies-Bouldin index: " << score.value << '\n';

        ReportWriter::write(reportPath, score);
        std::cout << "Report written to: " << reportPath << '\n';
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
