#pragma once

#include "qx/Core.hpp"
#include "qx/DefaultOperations.hpp"
#include "qx/SimulationResult.hpp"

#include <optional>
#include <string>

namespace qx {

struct SimulationError {
    std::string message = "Simulation failed";
};

std::unique_ptr<qx::core::MixedStateBase>
executeStringImpl(std::string const &s, Operations const &operations,
                  std::size_t iterations,
                  std::optional<std::uint_fast64_t> seed);

std::variant<SimulationResult, SimulationError>
executeString(std::string const &s,
              Operations operations = default_operations::defaultOperations,
              std::size_t iterations = 1,
              std::optional<std::uint_fast64_t> seed = std::nullopt);

std::unique_ptr<qx::core::MixedStateBase>
executeFileImpl(std::string const &filePath, Operations operations,
                std::size_t iterations, std::optional<std::uint_fast64_t> seed);

std::variant<SimulationResult, SimulationError>
executeFile(std::string const &filePath,
            Operations operations = default_operations::defaultOperations,
            std::size_t iterations = 1,
            std::optional<std::uint_fast64_t> seed = std::nullopt);

} // namespace qx