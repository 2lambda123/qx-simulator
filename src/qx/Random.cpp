#include "qx/Random.hpp"

namespace qx {
namespace random {

void seed(std::uint_fast64_t seedValue) {
    RandomNumberGenerator::getInstance().seed(seedValue);
}

double randomZeroOneDouble() {
    // std::uniform_real_distribution<double> does not give the same result
    // across platforms, so use this instead.

    double result = uniformMinMaxIntegerDistribution(
        0, UINT_FAST64_MAX, RandomNumberGenerator::getInstance()());
    assert(0. <= result && result <= 1.);
    return result;
}

std::uint_fast64_t randomInteger(std::uint_fast64_t min,
                                 std::uint_fast64_t max) {
    // std::uniform_int_distribution<std::uint_fast64_t> is not consistent
    // across platforms.
    assert(min <= max);
    assert(max - min < UINT_FAST64_MAX);

    std::uint_fast64_t numberOfBuckets = max - min + 1;
    std::uint_fast64_t bucketSize = UINT_FAST64_MAX / numberOfBuckets;
    std::uint_fast64_t limit = numberOfBuckets * bucketSize;

    assert(limit <= UINT_FAST64_MAX);

    std::uint_fast64_t r = 0;
    do {
        r = RandomNumberGenerator::getInstance()();
    } while (r >= limit);

    std::uint_fast64_t bucketIndex = r / bucketSize;

    assert(0 <= bucketIndex && bucketIndex < numberOfBuckets);

    std::size_t result = min + bucketIndex;

    assert(min <= result && result <= max);

    return result;
}

double uniformMinMaxIntegerDistribution(std::uint_fast64_t min,
                                        std::uint_fast64_t max, double x) {
    assert(min <= max);

    if (x < min) {
        return 0.;
    }

    if (x >= max) {
        return 1.;
    }

    return (std::floor(x - min) + 1) / (static_cast<double>(max - min) + 1);
}

double uniformZeroOneContinuousDistribution(double x) {
    if (x < 0) {
        return 0;
    }

    if (x > 1) {
        return 1;
    }

    return x;
}

} // namespace random
} // namespace qx