#include "qx/SimulationResult.hpp"

#include "qx/Core.hpp"
#include "qx/Random.hpp"
#include <iomanip>
#include <iostream>
#include <variant>

namespace qx {

std::ostream &operator<<(std::ostream &os, SimulationResult const &r) {
    os << std::setprecision(config::OUTPUT_DECIMALS) << std::fixed;
    os << "\33[2K\r--------------------------------------------------------------------------\n";

    os << "Final quantum state\n";

    if (r.state.empty()) {
        os << " < mixed state >\n";
    }

    for (auto const &kv : r.state) {
        auto const &bitString = kv.first;
        auto const &amplitude = kv.second;
        os << bitString << "       " << amplitude.real() << " + "
           << amplitude.imag() << "*i   "
           << " (norm = " << std::norm(amplitude) << ")\n";
    }

    os << std::endl << "Measurement register averaging\n";

    for (const auto &[bitString, probability] : r.measurementRegisterStatistics) {
        os << bitString << "       " << probability << ")\n";
    }
    return os;
}

SimulationResult generateSimulationResult(std::uint64_t nIterations, std::unique_ptr<core::MixedStateBase> quantumState) {
    SimulationResult simulationResult;

    simulationResult.iterations = nIterations;

    core::AssociationVectorStringMap<double> wrappedMeasurementRegisterStatistics(simulationResult.measurementRegisterStatistics);
    quantumState->getMeasurementRegisterStatistics(wrappedMeasurementRegisterStatistics);

    core::AssociationVectorStringMap<std::complex<double>> wrappedState(simulationResult.state);
    quantumState->getPureState(wrappedState);

    return simulationResult;
}
} // namespace qx