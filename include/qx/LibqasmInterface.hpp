#pragma once

#include "qx/Operations.hpp"
#include <memory>

namespace cqasm::v1::semantic {
class Subcircuit;
}

namespace qx {

class Circuit;

qx::Circuit loadCqasmCode(cqasm::v1::semantic::Subcircuit const &subcircuit, Operations const& operations, std::uint64_t qubitCount);

} // namespace qx