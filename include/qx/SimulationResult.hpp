#pragma once

#include <complex>
#include <string>
#include <vector>
#include <memory>

#include "qx/utils/FloatComparison.hpp"
#include "qx/Common.hpp"

namespace qx {

namespace core {
class MixedStateBase;
}

struct SimulationResult {
    using MeasurementRegisterStatistics = std::vector<std::pair<std::string, double>>;
    using State = std::vector<std::pair<std::string, std::complex<double>>>;
    
    std::uint64_t iterations = 0;
    MeasurementRegisterStatistics measurementRegisterStatistics;
    State state;
};

std::ostream &operator<<(std::ostream &os, SimulationResult const &r);

SimulationResult generateSimulationResult(std::uint64_t nIterations, std::unique_ptr<core::MixedStateBase> quantumState);

} // namespace qx