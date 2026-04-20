#ifndef RAND_UTILS_H
#define RAND_UTILS_H

#include <cmath>
#include <cstdlib>
#include <ctime>

inline double rand_normal()
{
    static bool seeded = false;
    if (!seeded) {
        std::srand(static_cast<unsigned>(std::time(nullptr))); // seed once
        seeded = true;
    }

    static const double pi = 3.141592653589793238;
    double u = 0;
    while (u == 0) // loop to ensure u nonzero, for log
    {
        u = rand() / static_cast<double>(RAND_MAX);
    }
    double v = rand() / static_cast<double>(RAND_MAX);
    return std::sqrt(-2.0 * std::log(u)) * std::cos(2.0 * pi * v);
}

#endif