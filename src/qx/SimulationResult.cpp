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
        // Can you write this as?
        auto const &[bitString, amplitude] = kv;

        os << bitString << "       " << amplitude.real() << " + "
           << amplitude.imag() << "*i   "
           << " (norm = " << std::norm(amplitude) << ")\n";
        // The fmt library is very convenient for these things
        // Actually, once you start using it, you don't want to use iostreams any more
        // It also supports Unicode, date/time formatting...
        fmt::print(os, "{}       {} + {}*i    (norm = {})\n",
                   bitString, amplitude.real(), amplitude.imag(), std::norm(amplitude));
        // Actually, std::print should be available with C++23
        // https://en.cppreference.com/w/cpp/io/print
        // So, same code as above, but with std::print
        // But I'm not sure it is yet as powerful as fmt
        // E.g. Unicode, date/time formatting...
    }

    os << "\nMeasurement register averaging\n";

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