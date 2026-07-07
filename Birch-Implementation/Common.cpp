#include "Common.h"

double distanceSquared(const Point& a, const Point& b)
{
    if(a.size() != b.size()) throw std::invalid_argument("Feature dimension mismatch");
    double sum = 0.0;
    for(size_t i=0; i<a.size(); ++i)
    {
        const double difference = a[i] - b[i];
        sum += difference * difference;
    }
    return sum;
}

double distance(const Point& a, const Point& b) { return std::sqrt(distanceSquared(a, b)); }
