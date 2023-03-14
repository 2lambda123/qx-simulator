#pragma once

#include "qx/utils/BasisVector.hpp"
#include "qx/utils/StrongTypes.hpp"

namespace qx {

template <std::uint64_t MaxNumberOfQubits = 64> struct Key {
    using BasisVector = utils::BasisVector<MaxNumberOfQubits>;

    EnsembleIndex ensembleIndex;
    BasisVector measurementRegister;
    BasisVector basisVector;

    auto operator<=>(Key const &) const = default;
};

} // namespace qx