#include "qx/SimulationResult.hpp"

#include "qx/Core.hpp"
#include "qx/Random.hpp"
#include <iomanip>
#include <iostream>
#include <variant>

namespace qx {

std::ostream &operator<<(std::ostream &os, SimulationResult const &r) {
    os << std::setprecision(config::OUTPUT_DECIMALS) << std::fixed;
    os << "\33[2K\r------------------------------------------------------------"
          "--------------\n";

    os << "Final quantum state\n";

    if (r.state.empty()) {
        os << " < mixed state >\n";
    }

    for (auto const &[bitString, amplitude] : r.state) {
        os << fmt::format("{}       {} + {}*i    (norm = {})\n",
                   bitString, amplitude.real(), amplitude.imag(), std::norm(amplitude));
    }

    os << "\nMeasurement register averaging\n";

    for (const auto &[bitString, probability] :
         r.measurementRegisterStatistics) {
        os << bitString << "       " << probability << "\n";
    }
    return os;
}

SimulationResult
generateSimulationResult(std::uint64_t nIterations,
                         std::unique_ptr<core::MixedStateBase> quantumState) {
    SimulationResult simulationResult;

    simulationResult.iterations = nIterations;

    core::AssociationVectorStringMap<double>
        wrappedMeasurementRegisterStatistics(
            simulationResult.measurementRegisterStatistics);
    quantumState->getMeasurementRegisterStatistics(
        wrappedMeasurementRegisterStatistics);

    core::AssociationVectorStringMap<std::complex<double>> wrappedState(
        simulationResult.state);
    quantumState->getPureState(wrappedState);

    return simulationResult;
}
} // namespace qx