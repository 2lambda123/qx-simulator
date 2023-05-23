#include "qx/Simulator.hpp"

#include "qx/Circuit.hpp"
#include "qx/LibqasmInterface.hpp"
#include "qx/Random.hpp"
#include "qx/SimulationResult.hpp"
#include "cqasm.hpp"

#include <fstream>
#include <iostream>
#include <vector>

namespace qx {

namespace {
using Program = cqasm::v1::ast::One<cqasm::v1::semantic::Program>;
using AnalysisResult = cqasm::v1::analyzer::AnalysisResult;

auto createAnalyzer() {
    cqasm::v1::analyzer::Analyzer analyzer{"1.2"};
    analyzer.register_default_functions_and_mappings();
    return analyzer;
}

std::string getErrors(AnalysisResult const& analysisResult) {
    assert(!analysisResult.errors.empty());

    std::stringstream errorsString;
    bool first = true;
    for (auto const& error: analysisResult.errors) {
        if (!first) {
            errorsString << "\n";
        } else {
            first = false;
        }
        errorsString << error;
    }

    return errorsString.str();
}

Program parseFile(std::string const &filePath) {
    auto res = createAnalyzer().analyze(filePath);

    if (!res.errors.empty()) {
        throw std::runtime_error(getErrors(res));
    }

    return res.root;
}

Program parseString(std::string const &s) {
    // Possible alternative code
    /*
    if (auto res = createAnalyzer().analyze_string(s); res.errors.empty()) {
        return res.root;
    } else {
        throw std::runtime_error(getErrors(res));
    }
    */
    auto res = createAnalyzer().analyze_string(s);

    if (!res.errors.empty()) {
        throw std::runtime_error(getErrors(res));
    }

    return res.root;
}

Program unwrap(auto parseResult, std::stringstream& errors) { // throws when analysis fails
    return parseResult.unwrap(errors);
}

std::unique_ptr<qx::core::MixedStateBase>
execute(Program program, Operations const& operations, std::size_t iterations,
                   std::optional<std::uint_fast64_t> seed) {
    if (program.empty()) {
        throw std::runtime_error("Program is empty");
    }
    if (!program->error_model.empty()) {
        throw std::runtime_error("Probabilistic error models are not supported");
    }
    if (seed) {
        random::seed(*seed);
    }

    std::vector<qx::Circuit> circuits;
    std::uint64_t qubitCount = program->num_qubits;
    for (auto const &subcircuit : program->subcircuits) {
        circuits.push_back(loadCqasmCode(*subcircuit, operations, qubitCount));
    }

    std::unique_ptr<qx::core::MixedStateBase> quantumState;
    if (qubitCount <= 64) {
        quantumState = std::make_unique<qx::core::MixedState<64>>(qubitCount);
    } else if (qubitCount <= 128) {
        quantumState = std::make_unique<qx::core::MixedState<128>>(qubitCount);
    } else if (qubitCount <= 256) {
        quantumState = std::make_unique<qx::core::MixedState<256>>(qubitCount);
    } else if (qubitCount <= 512) {
        quantumState = std::make_unique<qx::core::MixedState<512>>(qubitCount);
    } else {
        throw std::runtime_error("Cannot handle that many qubits in this version of QX-simulator");
    }

    for (auto &circuit : circuits) {
        circuit.execute(*quantumState);
    }
    return quantumState;
}
}

std::unique_ptr<qx::core::MixedStateBase>
executeStringImpl(std::string const &s, Operations const& operations, std::size_t iterations,
                  std::optional<std::uint_fast64_t> seed) {
    auto program = parseString(s);
    return execute(program, operations, iterations, seed);
}

std::variant<SimulationResult, SimulationError>
executeString(std::string const &s, Operations operations, std::size_t iterations,
              std::optional<std::uint_fast64_t> seed) {
    try {
        auto quantumState = executeStringImpl(s, operations, iterations, seed);
        return generateSimulationResult( iterations, std::move(quantumState));
    } catch (std::exception const& e) {
        return SimulationError{e.what()};
    }

}

std::unique_ptr<qx::core::MixedStateBase>
executeFileImpl(std::string const &filePath, Operations operations, std::size_t iterations,
                std::optional<std::uint_fast64_t> seed) {
    auto program = parseFile(filePath);
    return execute(std::move(program), operations, iterations, seed);
}

std::variant<SimulationResult, SimulationError>
executeFile(std::string const &filePath, Operations operations, std::size_t iterations,
            std::optional<std::uint_fast64_t> seed) {
    try {
        auto quantumState = executeFileImpl(filePath, operations, iterations, seed);
        return generateSimulationResult(iterations, std::move(quantumState));
    } catch (std::exception const& e) {
        return SimulationError{e.what()};
    }
}

} // namespace qx