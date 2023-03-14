#pragma once

#include "qx/Simulator.hpp"
#include "qx/Operations.hpp"

namespace qxelarator {

std::unique_ptr<qx::core::MixedStateBase>
execute_string(std::string const &s, qx::Operations operations = qx::default_operations::defaultOperations, std::size_t iterations = 1,
               std::optional<std::uint_fast64_t> seed = std::nullopt) {
    return qx::executeStringImpl(s, operations, iterations, seed);
}

std::unique_ptr<qx::core::MixedStateBase>
execute_file(std::string const &filePath, qx::Operations operations = qx::default_operations::defaultOperations, std::size_t iterations = 1,
             std::optional<std::uint_fast64_t> seed = std::nullopt) {
    return qx::executeFileImpl(filePath, operations, iterations, seed);
}

} // namespace qxelarator