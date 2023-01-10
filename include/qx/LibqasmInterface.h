#pragma once

#include <memory>

namespace cqasm::v1::semantic {
    class Subcircuit;
}

namespace qx {

class Circuit;

std::shared_ptr<qx::Circuit> loadCqasmCode(uint64_t qubits_count,
                                            cqasm::v1::semantic::Subcircuit const &subcircuit);

}