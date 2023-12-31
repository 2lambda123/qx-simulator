#include "qx/core/gaussian.h"

/**
 * random
 */
qx::gaussian::random::random(double mean, double deviation)
    : mean(mean), deviation(deviation), distribution(mean, deviation) {}

/**
 * random
 */
double qx::gaussian::random::next() { return distribution(generator); }

/**
 * random
 */
void qx::gaussian::random::set_mean(double mean) { this->mean = mean; }

/**
 * random
 */
void qx::gaussian::random::set_std_deviation(double deviation) {
    this->deviation = deviation;
}
