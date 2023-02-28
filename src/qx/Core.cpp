#include "qx/Core.hpp"

namespace qx {
namespace core {
namespace {
template <std::size_t NumberOfOperands>
void applyImpl(DenseUnitaryMatrix<1 << NumberOfOperands> const &matrix,
               std::array<QubitIndex, NumberOfOperands> const &operands,
               BasisVector index, std::complex<double> value,
               SparseArray::Map &storage) {
    utils::Bitset<NumberOfOperands> reducedIndex;
    for (std::size_t i = 0; i < NumberOfOperands; ++i) {
        reducedIndex.set(i,
                         index.test(operands[NumberOfOperands - i - 1].value));
    }

    for (std::size_t i = 0; i < (1 << NumberOfOperands); ++i) {
        std::complex<double> addedValue =
            value * matrix.at(i, reducedIndex.toSizeT());

        if (isNotNull(addedValue)) {
            auto newIndex = index;

            for (std::size_t k = 0; k < NumberOfOperands; ++k) {
                newIndex.set(operands[NumberOfOperands - k - 1].value,
                             utils::getBit(i, k));
            }

            auto it = storage.try_emplace(newIndex, 0);
            auto newValue = it.first->second + addedValue;
            it.first->second = newValue;
        }
    }
}
} // namespace

void SparseArray::set(BasisVector index, std::complex<double> value) {
#ifndef NDEBUG
    if (index.toSizeT() >= size) {
        throw std::runtime_error("SparseArray::set index out of bounds");
    }
#endif

    if (std::abs(value) < config::EPS) {
        return;
    }

    data.try_emplace(index, value);
}

void SparseArray::cleanupZeros() {
    absl::erase_if(data, [](auto const &kv) { return !isNotNull(kv.second); });
    zeroCounter = 0;
}

void QuantumState::testInitialize(
    std::initializer_list<std::pair<std::string, std::complex<double>>>
        values) {
    data.clear();
    double norm = 0;
    for (auto const &kv : values) {
        BasisVector index(kv.first);
        data.set(index, kv.second);
        norm += std::norm(kv.second);
    }
    assert(!isNotNull(norm - 1));
};

template <std::size_t NumberOfOperands>
QuantumState &
QuantumState::apply(DenseUnitaryMatrix<1 << NumberOfOperands> const &m,
                    std::array<QubitIndex, NumberOfOperands> const &operands) {
    assert(NumberOfOperands <= numberOfQubits &&
           "Quantum gate has more operands than the number of qubits in this "
           "quantum state");
    assert(std::find_if(operands.begin(), operands.end(),
                        [this](auto qubitIndex) {
                            return qubitIndex.value >= numberOfQubits;
                        }) == operands.end() &&
           "Operand refers to a non-existing qubit");

    data.applyLinear(std::bind(&applyImpl<NumberOfOperands>, m, operands,
                               std::placeholders::_1, std::placeholders::_2,
                               std::placeholders::_3));

    return *this;
}

// Explicit instantiation for use in Circuit::execute, otherwise linking error.

template QuantumState &
QuantumState::apply<1>(DenseUnitaryMatrix<1 << 1> const &m,
                       std::array<QubitIndex, 1> const &operands);

template QuantumState &
QuantumState::apply<2>(DenseUnitaryMatrix<1 << 2> const &m,
                       std::array<QubitIndex, 2> const &operands);

template QuantumState &
QuantumState::apply<3>(DenseUnitaryMatrix<1 << 3> const &m,
                       std::array<QubitIndex, 3> const &operands);


void StatisticalEnsemble::shrink() {
    assert(isConsistent());

    absl::btree_set<BasisVector> usedKets;
    for (auto const& [key, complexAmplitude]: data) {
        usedKets.insert(key.basisVector);
    }
    auto numberOfUsedKets = usedKets.size();

    std::size_t const& dimension1 = numberOfUsedKets;
    std::size_t dimension2 = std::min(ensembleSize, numberOfUsedKets);
    std::size_t maxSize = (2 * dimension1 - dimension2) * (dimension2 + 1) / 2;

    if (data.size() <= maxSize) {
        return;
    }

    std::size_t i = 0;
    while(i < std::min(ensembleSize, numberOfUsedKets)) {
        shrinkImpl(i, usedKets);

        std::size_t newEnsembleSize = 0;
        absl::erase_if(data, [&](auto const& kv) {
            if (isNotNull(kv.second)) {
                newEnsembleSize = std::max(newEnsembleSize, kv.first.ensembleIndex + 1);
                return false;
            }
            return true;
        });

        ensembleSize = newEnsembleSize;
        ++i;
    }

    assert(data.size() <= maxSize);
}

bool StatisticalEnsemble::isConsistent() {
    double acc = 0.;
    
    for (auto const& [key, complexAmplitude]: data) {
        acc += std::norm(complexAmplitude);
    }

    return !isNotNull(acc - 1.);
}

void StatisticalEnsemble::shrinkImpl(std::size_t minimumEnsembleIndex, absl::btree_set<BasisVector> &possibleKets) {
    assert(!data.empty());

    auto isValidKey = [&] (auto const& k) {
        return k.ensembleIndex >= minimumEnsembleIndex &&
            possibleKets.count(k.basisVector) != 0;
    };

    if (possibleKets.empty()) {
        return;
    }

    auto candidate = *possibleKets.begin();

    double xSquaredNormWithoutFirstElement = 0.;
    
    for (auto const& [key, complexAmplitude]: data) {
        if (!isValidKey(key)) {
            continue;
        }

        if (key.ensembleIndex != minimumEnsembleIndex && key.basisVector == candidate) {
            xSquaredNormWithoutFirstElement += std::norm(complexAmplitude);
        }
    }

    auto firstElementIt = data.find({candidate, minimumEnsembleIndex});
    std::complex<double> xFirstElement = firstElementIt == data.end() ? 0. : firstElementIt->second;

    double alpha = std::sqrt(xSquaredNormWithoutFirstElement + std::norm(xFirstElement)); // FIXME: set the complex argument of alpha to improve floating point errors.  

    auto u = [&](std::size_t i) -> std::complex<double> {
        assert(i >= minimumEnsembleIndex);
        auto dataIt = data.find({candidate, i});
        if (i != minimumEnsembleIndex && dataIt == data.end()) {
            return 0.;
        } else if (i != minimumEnsembleIndex) {
            return dataIt->second;
        } else {
            assert(i == minimumEnsembleIndex);
            return xFirstElement - alpha;
        }
    };

    double uInvSquaredNorm = 1 / (xSquaredNormWithoutFirstElement + std::norm(xFirstElement - alpha));

    auto Q = [&](std::size_t i, std::size_t j) -> std::complex<double> {
        assert(i >= minimumEnsembleIndex && j >= minimumEnsembleIndex);
        std::complex<double> result = i == j ? 1. : 0.;

        result -= 2. * u(i) * std::conj(u(j)) * uInvSquaredNorm;

        return result;
    };

    UnderlyingT newData;
    newData[{ .basisVector = candidate, .ensembleIndex = minimumEnsembleIndex}] = alpha;

    for (auto const& [key, complexAmplitude]: data) {
        if (!isValidKey(key)) {
            newData[key] = complexAmplitude; // FIXME: a lot of copy. Maybe erase+merge?
            continue;
        }

        auto const& [basisVector, ensembleIndex] = key;

        if (basisVector == candidate) {
            continue;
        }

        for (std::size_t i = minimumEnsembleIndex; i < ensembleSize; ++i) {
            newData[{basisVector, i}] += Q(i, ensembleIndex) * complexAmplitude;
        }
    }

    possibleKets.erase(candidate);
    data.swap(newData);
}


} // namespace core
} // namespace qx