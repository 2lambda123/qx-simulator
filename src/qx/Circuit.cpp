#include "qx/Circuit.hpp"

#include <algorithm>

namespace qx {

void Circuit::execute(core::MixedStateBase &quantumState) const {
    std::uint64_t it = iterations;
    while (it-- > 0) {
        for (auto const& instruction: instructions) {
            quantumState(instruction);
        }
    }
}

} // namespace qx