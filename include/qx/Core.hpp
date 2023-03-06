#pragma once

#include "absl/container/flat_hash_map.h"
#include "absl/container/btree_set.h"
#include "absl/container/inlined_vector.h"
#include <cassert>
#include <complex>
#include <limits>

#include "qx/Common.hpp"
#include "qx/CompileTimeConfiguration.hpp"

namespace qx {
namespace core {

inline constexpr bool isNotNull(std::complex<double> c) {
#if defined(_MSC_VER)
    return c.real() > config::EPS || -c.real() > config::EPS ||
           c.imag() > config::EPS || -c.imag() > config::EPS;
#else
    return std::abs(c.real()) > config::EPS || std::abs(c.imag()) > config::EPS;
#endif
}

struct QubitIndex {
    std::size_t value;
};

template <std::size_t N> class DenseMatrix {
public:
    using UnderlyingT = std::array<std::array<std::complex<double>, N>, N>;

    static constexpr DenseMatrix<N> identity() {
        UnderlyingT m;
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = 0; j < N; ++j) {
                m[i][j] = (i == j) ? 1 : 0;
            }
        }

        return DenseMatrix(m);
    }

    constexpr DenseMatrix() : matrix() {};

    explicit constexpr DenseMatrix(UnderlyingT const &m)
        : matrix(m) {}

    inline constexpr std::complex<double> at(std::size_t i,
                                             std::size_t j) const {
        return matrix[i][j];
    }

    inline constexpr std::complex<double>& at(std::size_t i, std::size_t j) {
        return matrix[i][j];
    }

    constexpr DenseMatrix<N> dagger() const {
        UnderlyingT m;
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = 0; j < N; ++j) {
                m[i][j] = std::conj(at(j, i));
            }
        }

        return DenseMatrix(m);
    }

    constexpr bool operator==(DenseMatrix<N> const &other) const {
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = 0; j < N; ++j) {
                if (isNotNull(at(i, j) - other.at(i, j))) {
                    return false;
                }
            }
        }
        return true;
    }

    constexpr DenseMatrix<N>
    operator*(DenseMatrix<N> const &other) const {
        UnderlyingT m;
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = 0; j < N; ++j) {
                m[i][j] = 0;
                for (std::size_t k = 0; k < N; ++k) {
                    m[i][j] += at(i, k) * other.at(k, j);
                }
            }
        }

        return DenseMatrix(m);
    }

    constexpr void
    operator+=(DenseMatrix<N> const &other) {
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = 0; j < N; ++j) {
                at(i, j) += other.at(i, j);
            }
        }
    }

private:
    std::array<std::array<std::complex<double>, N>, N> matrix;
};

template <std::size_t N>
void operator<<(std::ostream &os, DenseMatrix<N> const &m) {
    for (std::size_t i = 0; i < N; ++i) {
        bool first = true;
        for (std::size_t j = 0; j < N; ++j) {
            if (!first) {
                os << "  ";
            } else {
                first = false;
            }

            os << m.at(i, j);
        }
        os << std::endl;
    }
}

template <std::size_t N>
constexpr DenseMatrix<N> operator*(double d, DenseMatrix<N> m) {
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < N; ++j) {
            m.at(i, j) *= d;
        }
    }

    return m;
}

template <std::size_t N> class DenseUnitaryMatrix : public DenseMatrix<N> {
public:
    explicit constexpr DenseUnitaryMatrix(typename DenseMatrix<N>::UnderlyingT const &m)
        : DenseMatrix<N>(m) {
            if (!(*this * DenseMatrix<N>::dagger() == DenseMatrix<N>::identity())) {
                throw std::runtime_error("Matrix is not unitary");
            };
    }

    constexpr DenseUnitaryMatrix(DenseMatrix<N> const &m)
        : DenseMatrix<N>(m) {
            if (!(*this * DenseMatrix<N>::dagger() == DenseMatrix<N>::identity())) {
                throw std::runtime_error("Matrix is not unitary");
            };
    }
};

class QuantumState;

class SparseArray {
public:
    using Map = absl::flat_hash_map<BasisVector, std::complex<double>>;
    using Iterator = Map::const_iterator;

    SparseArray() = delete;

    explicit SparseArray(std::size_t s) : size(s){};

    std::size_t getSize() const { return size; }

    std::vector<std::complex<double>> testToVector() const {
        std::vector<std::complex<double>> result(getSize(), 0);

        for (auto const &kv : *this) {
            result[kv.first.toSizeT()] = kv.second;
        }

        return result;
    }

    Iterator begin() const { return data.cbegin(); }

    Iterator end() const { return data.cend(); }

    void set(BasisVector index, std::complex<double> value);

    void clear() { data.clear(); }

    SparseArray &operator*=(double d) {
        std::for_each(data.begin(), data.end(),
                      [d](auto &kv) { kv.second *= d; });
        return *this;
    }

    template <typename F> void forEach(F &&f) {
        cleanupZeros();
        std::for_each(data.begin(), data.end(), f);
    }

    template <typename F> void forEachSorted(F &&f) {
        cleanupZeros();
        std::vector<std::pair<BasisVector, std::complex<double>>> sorted(
            data.begin(), data.end());
        std::sort(sorted.begin(), sorted.end(),
                  [](auto const &left, auto const &right) {
                      return left.first < right.first;
                  });
        std::for_each(sorted.begin(), sorted.end(), f);
    }

    template <typename F> void eraseIf(F &&pred) { absl::erase_if(data, pred); }

private:
    friend QuantumState;

    // Let f build a new SparseArray to replace *this, assuming f is linear.
    template <typename F> void applyLinear(F &&f) {
        // Every ZERO_CYCLE_SIZE gates, cleanup the 0s
        if (zeroCounter >= config::ZERO_CYCLE_SIZE) {
            cleanupZeros();
        }
        ++zeroCounter;

        Map result;

        for (auto const &kv : data) {
            f(kv.first, kv.second, result);
        }

        data.swap(result);
    }

    void cleanupZeros();

    std::size_t const size = 0;
    std::uint64_t zeroCounter = 0;
    Map data;
};

class QuantumState {
public:
    explicit QuantumState(std::size_t n)
        : numberOfQubits(n), data(1 << numberOfQubits) {
        assert(numberOfQubits > 0 && "QuantumState needs at least one qubit");
        assert(numberOfQubits <= config::MAX_QUBIT_NUMBER &&
               "QuantumState currently cannot support that many qubits with "
               "this version of OpenQL");
        data.set(BasisVector{}, 1); // Start initialized in state 00...000
    };

    std::size_t getNumberOfQubits() const { return numberOfQubits; }

    void reset() {
        data.clear();
        data.set(BasisVector{}, 1); // Start initialized in state 00...000
        measurementRegister.reset();
    }

    void testInitialize(
        std::initializer_list<std::pair<std::string, std::complex<double>>>
            values);

    template <std::size_t NumberOfOperands>
    QuantumState &
    apply(DenseUnitaryMatrix<1 << NumberOfOperands> const &m,
          std::array<QubitIndex, NumberOfOperands> const &operands);

    template <typename F> void forEach(F &&f) { data.forEachSorted(f); }

    BasisVector getMeasurementRegister() const { return measurementRegister; }

    BasisVector &getMeasurementRegister() { return measurementRegister; }

    template <typename F>
    void measure(QubitIndex qubitIndex, F &&randomGenerator) {
        auto rand = randomGenerator();
        double probabilityOfMeasuringOne = 0.;

        data.forEach([qubitIndex, &probabilityOfMeasuringOne](auto const &kv) {
            if (kv.first.test(qubitIndex.value)) {
                probabilityOfMeasuringOne += std::norm(kv.second);
            }
        });

        if (rand < probabilityOfMeasuringOne) {
            data.eraseIf([qubitIndex](auto const &kv) {
                return !kv.first.test(qubitIndex.value);
            });
            data *= std::sqrt(1 / probabilityOfMeasuringOne);
            measurementRegister.set(qubitIndex.value, true);
        } else {
            data.eraseIf([qubitIndex](auto const &kv) {
                return kv.first.test(qubitIndex.value);
            });
            data *= std::sqrt(1 / (1 - probabilityOfMeasuringOne));
            measurementRegister.set(qubitIndex.value, false);
        }
    }

    template <typename F> void measureAll(F &&randomGenerator) {
        auto rand = randomGenerator();
        double probability = 0.;

        auto measuredState = std::invoke([this, &probability, rand] {
            for (auto const &kv :
                 data) { // Does this work with non-ordered iteration?
                probability += std::norm(kv.second);
                if (probability > rand) {
                    return kv;
                }
            }
            throw std::runtime_error(
                "Vector was not normalized at measurement location (a bug)");
        });

        data.clear();
        data.set(measuredState.first,
                 measuredState.second / std::abs(measuredState.second));
        measurementRegister = measuredState.first;
    }

    template <typename F>
    void prep(QubitIndex qubitIndex, F &&randomGenerator) {
        // Measure + conditional X, and reset the measurement register.
        auto rand = randomGenerator();
        double probabilityOfMeasuringOne = 0.;

        data.forEach([qubitIndex, &probabilityOfMeasuringOne](auto const &kv) {
            if (kv.first.test(qubitIndex.value)) {
                probabilityOfMeasuringOne += std::norm(kv.second);
            }
        });

        if (rand < probabilityOfMeasuringOne) {
            data.eraseIf([qubitIndex](auto const &kv) {
                return !kv.first.test(qubitIndex.value);
            });
            SparseArray::Map newData;
            for (auto kv : data.data) {
                auto newKey = kv.first;
                newKey.set(qubitIndex.value, false);
                newData.insert(std::make_pair(
                    newKey,
                    kv.second * std::sqrt(1 / probabilityOfMeasuringOne)));
            }
            data.data = newData; // Could fix the interface
        } else {
            data.eraseIf([qubitIndex](auto const &kv) {
                return kv.first.test(qubitIndex.value);
            });
            data *= std::sqrt(1 / (1 - probabilityOfMeasuringOne));
        }
        measurementRegister.set(qubitIndex.value, false);
    };

private:
    std::size_t const numberOfQubits = 1;
    SparseArray data;
    BasisVector measurementRegister{};
};

class StatisticalEnsemble {
public:
    template <std::size_t NumberOfOperands>
    static bool areValidKrausOperators(std::initializer_list<DenseMatrix<1 << NumberOfOperands>> krausOperators) {
        // Check completeness condition.
        DenseMatrix<1 << NumberOfOperands> sum;
        for (auto const& m: krausOperators) {
            sum += m.dagger() * m;
        }
        return sum == DenseMatrix<1 << NumberOfOperands>::identity();
    }

    StatisticalEnsemble(std::size_t n) : numberOfQubits(n), ensembleSize(1) {
        data.insert({KeyT{.basisVector = BasisVector(), .ensembleIndex = 0}, 1.});
    }

    constexpr bool operator==(StatisticalEnsemble const &other) const {
        return numberOfQubits == other.numberOfQubits && ensembleSize == other.ensembleSize && data == other.data; // FIXME
    }

    template <std::size_t NumberOfOperands>
    void applyKrausOperators(std::initializer_list<DenseMatrix<1 << NumberOfOperands>> krausOperators, std::array<QubitIndex, NumberOfOperands> const &operands) {
        assert(isConsistent());

        assert(areValidKrausOperators<NumberOfOperands>(krausOperators) && "Kraus operators don't satisfy completeness constraint");

        UnderlyingT newData;

        for (auto const& [key, amplitude]: data) {
            std::size_t operatorIndex = 0;
            for (auto krausOperatorIt = krausOperators.begin(); krausOperatorIt != krausOperators.end(); ++krausOperatorIt) {
                auto const& [basisVector, ensembleIndex] = key;

                utils::Bitset<NumberOfOperands> reducedBasisVector;
                for (std::size_t i = 0; i < NumberOfOperands; ++i) {
                    reducedBasisVector.set(i, basisVector.test(operands[NumberOfOperands - i - 1].value));
                }

                for (std::size_t i = 0; i < (1 << NumberOfOperands); ++i) {
                    auto const& krausValue = krausOperatorIt->at(i, reducedBasisVector.toSizeT());

                    if (krausValue == 0.) { // No fancy double comparison here.
                        continue;
                    }
                    
                    std::complex<double> addedValue = amplitude * krausValue;

                    auto modifiedBasisVector = basisVector;

                    for (std::size_t k = 0; k < NumberOfOperands; ++k) {
                        modifiedBasisVector.set(operands[NumberOfOperands - k - 1].value,
                                    utils::getBit(i, k));
                    }

                    auto it = newData.try_emplace({ .basisVector = modifiedBasisVector, .ensembleIndex = ensembleIndex * krausOperators.size() + operatorIndex }, 0);
                    it.first->second += addedValue;
                }

                ++operatorIndex;
            }
        }

        ensembleSize = ensembleSize * krausOperators.size();
        data.swap(newData);
    }

    template <std::size_t NumberOfOperands>
    absl::InlinedVector<double, 2> applyKrausOperatorsGetOutcomeProbabilities(std::initializer_list<DenseMatrix<1 << NumberOfOperands>> krausOperators, std::array<QubitIndex, NumberOfOperands> const &operands) {
        applyKrausOperators(krausOperators, operands);

        absl::InlinedVector<double, 2> result(krausOperators.size(), 0.);

        for (auto const& [key, complexAmplitude]: data) { // FIXME: this relies on the assignation of ensemble indices by applyKrausOperators...
            auto const& ensembleIndex = key.ensembleIndex;

            result[ensembleIndex % krausOperators.size()] += std::norm(complexAmplitude);
        }

        return result;
    }

    std::size_t getEnsembleSize() const {
        return ensembleSize;
    }

    void shrink();

    template <std::size_t N>
    DenseMatrix<N> testToDensityMatrix() const {
        if (1 << numberOfQubits != N) {
            throw std::runtime_error("Density matrix size must match number of qubits");
        }

        DenseMatrix<N> result;

        for (std::size_t k = 0; k < ensembleSize; ++k) {
            for (std::size_t i = 0; i < N; ++i) {
                auto row = data.find(KeyT{BasisVector(i), k});
                if (row == data.end()) {
                    continue;
                }

                for (std::size_t j = 0; j < N; ++j) {
                    auto column = data.find(KeyT{BasisVector(j), k});
                    if (column == data.end()) {
                        continue;
                    }

                    result.at(i, j) += row->second * std::conj(column->second);
                }
            }
        }

        return result;
    }

    template <std::size_t N>
    std::vector<std::array<std::complex<double>, N>> testToMatrix() const {
        if (1 << numberOfQubits != N) {
            throw std::runtime_error("Density matrix size must match number of qubits");
        }

        std::vector<std::array<std::complex<double>, N>> result(ensembleSize, std::array<std::complex<double>, N>());
        
        for (auto const& [key, complexAmplitude]: data) {
            auto const& [basisVector, ensembleIndex] = key;
            assert(ensembleIndex < ensembleSize);
            result[ensembleIndex][basisVector.toSizeT()] = complexAmplitude;
        }

        return result;
    }

private:
    struct KeyT {
        BasisVector basisVector;
        std::size_t ensembleIndex = 0;
        
        template <typename H> friend H AbslHashValue(H h, KeyT const &key) {
            return H::combine(std::move(h), key.basisVector, key.ensembleIndex);
        }

        bool operator==(KeyT const& other) const {
            return ensembleIndex == other.ensembleIndex && basisVector == other.basisVector;
        }
    };

    using UnderlyingT = absl::flat_hash_map<KeyT, std::complex<double>>;

    bool isConsistent();
    
    void shrinkImpl(std::size_t minimumEnsembleIndex, absl::btree_set<BasisVector> &possibleKets);

    UnderlyingT data;
    std::size_t const numberOfQubits = 1;
    std::size_t ensembleSize = 1;
};

} // namespace core
} // namespace qx