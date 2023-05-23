#pragma once

#include "qx/Core.hpp"

#include <optional>
#include <string>
#include <vector>
#include <variant>

namespace qx {
class Circuit {
public:
    explicit Circuit(std::uint64_t n, std::string name = "", std::uint64_t iterations = 1)
        : numberOfQubits(n), name(name), iterations(iterations) {}

    void addInstruction(CircuitInstruction const& circuitInstruction) {
        instructions.push_back(circuitInstruction);
    }

    void execute(core::MixedStateBase &quantumState) const;

    std::string getName() const { return name; }

    std::uint64_t getNumberOfQubits() const { return numberOfQubits; }

private:
    std::vector<CircuitInstruction> instructions;
    std::uint64_t numberOfQubits = 0;
    // Consider not using const members, as they can be problematic (e.g. for move operations)
    // For this case, if 'name' and 'iterations' are private, and you don't provide set accessors,
    // they cannot be modified anyway, can they?
    std::string const name;
    std::uint64_t const iterations = 1;
};
} // namespace qx