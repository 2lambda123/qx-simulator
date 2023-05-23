#pragma once

#include "qx/utils/StrongTypes.hpp"
#include "qx/utils/BasisVector.hpp"

namespace qx {

template <std::uint64_t MaxNumberOfQubits = 64>
struct Key {
    using BasisVector = utils::BasisVector<MaxNumberOfQubits>;

    EnsembleIndex ensembleIndex;
    BasisVector measurementRegister;
    BasisVector basisVector;

    // Maybe worth a look if these two functions can just be replaced by a default spaceship operator
    // auto operator<=>(Key const &other) = default;
    // I guess it only needs EnsembleIndex and BasisVector to implement that operator (not all STL types yet implement it)

    bool operator==(Key const& other) const {
        return ensembleIndex == other.ensembleIndex &&
               measurementRegister == other.measurementRegister &&
               basisVector == other.basisVector;
    }

    bool operator<(Key const& other) const {
        if (ensembleIndex == other.ensembleIndex) {
            if (measurementRegister == other.measurementRegister) {
                return basisVector < other.basisVector;
            }
            return measurementRegister < other.measurementRegister;
        }

        return ensembleIndex < other.ensembleIndex;
    }
};

}