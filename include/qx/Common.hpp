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

    bool operator==(Key const& other) const {
        return ensembleIndex == other.ensembleIndex && measurementRegister == other.measurementRegister && basisVector == other.basisVector;
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