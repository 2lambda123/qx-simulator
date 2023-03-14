#include "qx/Core.hpp"

namespace qx {
namespace core {

template <std::uint64_t MaxNumberOfQubits>
bool MixedState<MaxNumberOfQubits>::isConsistent() const {
    if (data.empty()) {
        return false;
    }

    double accumulator = 0.;

    for (auto const& [key, factor]: data) {
        accumulator += std::norm(factor);
    }

    return utils::isNull(accumulator - 1.);
}

template bool MixedState<64>::isConsistent() const;
template bool MixedState<128>::isConsistent() const;
template bool MixedState<256>::isConsistent() const;
template bool MixedState<512>::isConsistent() const;

} // namespace core
} // namespace qx