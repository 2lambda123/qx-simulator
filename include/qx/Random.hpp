#pragma once

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <random>

#include "absl/container/btree_map.h"
#include "qx/utils/FloatComparison.hpp"

namespace qx {
namespace random {

class RandomNumberGenerator {
public:
    using RandomNumberGeneratorType = std::mt19937_64;

    static RandomNumberGeneratorType &getInstance() {
        static RandomNumberGenerator instance;
        return instance.randomNumberGenerator;
    }

    RandomNumberGenerator(RandomNumberGenerator const &) = delete;

    RandomNumberGenerator& operator=(RandomNumberGenerator const &) = delete;

private:
    RandomNumberGenerator() {
        std::random_device rd;
        randomNumberGenerator.seed(rd());
    }

    RandomNumberGeneratorType randomNumberGenerator;
};

void seed(std::uint_fast64_t seedValue);

double randomZeroOneDouble();

std::uint_fast64_t randomInteger(std::uint_fast64_t min,
                                 std::uint_fast64_t max);

template <typename K>
absl::btree_map<K, std::uint64_t>
randomMultinomial(std::uint64_t N, absl::btree_map<K, double> probabilities) {
    absl::btree_map<K, double> probabilitiesDistribution;

    double acc = 0.;
    for (auto p : probabilities) {
        acc += p.second;
        probabilitiesDistribution[p.first] = acc;
    }
    assert(utils::isNull(acc - 1));

    absl::btree_map<K, std::uint64_t> result;

    for (std::uint64_t i = 0; i < N; ++i) {
        auto rand = RandomNumberGenerator::getInstance()();

        auto it = probabilitiesDistribution.begin();
        while (rand > it->second * static_cast<double>(UINT_FAST64_MAX)) {
            ++it;
        }

        ++result[it->first];
    }

    return result;
}

double uniformMinMaxIntegerDistribution(std::uint_fast64_t min,
                                        std::uint_fast64_t max, double x);

double uniformZeroOneContinuousDistribution(double x);

// TODO: add function for Bernoulli and use that in measure.

} // namespace random
} // namespace qx