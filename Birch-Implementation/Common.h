#ifndef BIRCH_COMMON_H
#define BIRCH_COMMON_H

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using Point = std::vector<double>;

double distanceSquared(const Point& a, const Point& b);
double distance(const Point& a, const Point& b);

#endif
