#include "qx/Circuit.hpp"

#include <algorithm>

namespace qx {

void Circuit::execute(core::MixedStateBase &quantumState) const {
    std::uint64_t it = iterations;
    while (it-- > 0) {
        for (std::uint64_t instructionIndex = 0; instructionIndex < instructions.size(); ++instructionIndex) {
            quantumState(instructions[instructionIndex]);
        }
    }
}

} // namespace qx