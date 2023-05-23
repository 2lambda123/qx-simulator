#include "qx/Circuit.hpp"

#include <algorithm>
#include <numeric>  // iota

namespace qx {

void Circuit::execute(core::MixedStateBase &quantumState) const {
    std::range::for_each(
        std::views::iota(0, iterations),
        [](auto& i) { quantumState(instructions[i]); }
    );
}

} // namespace qx