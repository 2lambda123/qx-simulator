#pragma once

#include "absl/container/btree_map.h"
#include "absl/container/inlined_vector.h"
#include <cassert>
#include <complex>
#include <limits>

#include "qx/Common.hpp"
#include "qx/CompileTimeConfiguration.hpp"
#include "qx/MixedStateSimplifier.hpp"
#include "qx/utils/FloatComparison.hpp"
#include "qx/utils/StrongTypes.hpp"
#include "qx/Common.hpp"
#include "qx/CircuitInstruction.hpp"
#include "qx/SimulationResult.hpp"

namespace qx {
namespace core {

template <typename V>
struct StringMap {
    virtual void add(std::string const& key, V value) = 0;
};

template <typename V>
class AssociationVectorStringMap : public StringMap<V> {
public:
    using T = std::vector<std::pair<std::string, V>>;

    explicit AssociationVectorStringMap(T& v) : data(v) {
        assert(data.empty());
    }

    void add(std::string const& key, V value) override {
        // Could be sorted.
        auto it = std::ranges::find_last_if(data, [&key](auto const& kv) {
            return kv.first == key;
        });

        if (it == data.rend()) {
            data.push_back(std::make_pair(key, value));
        } else {
            it->second += value;
        }
    }

private:
    T& data;
};

template <std::uint64_t MaxNumberOfQubits>
class MixedState;

class MixedStateBase {
public:
    virtual void operator()(CircuitInstruction const& circuitInstruction) = 0;
    virtual void getPureState(StringMap<std::complex<double>>& m) const = 0;
    virtual void getMeasurementRegisterStatistics(StringMap<double>& m) const = 0;
};

template <std::uint64_t MaxNumberOfQubits = 64>
class MixedState final : public MixedStateBase {
public:
    using BasisVector = utils::BasisVector<MaxNumberOfQubits>;
    using OperandsVector = absl::InlinedVector<QubitIndex, config::MAX_INLINED_OPERANDS>;
    using ClassicalBitReferenceVector = absl::InlinedVector<utils::BitReference, config::MAX_INLINED_OPERANDS>;
    //
    using BTreeMapKeyDouble = absl::btree_map<Key<MaxNumberOfQubits>, std::complex<double>>;

    explicit MixedState(std::size_t n) : numberOfQubits(n) {
        data.insert({Key<MaxNumberOfQubits>{.ensembleIndex = {}, .measurementRegister = {}, .basisVector = {}}, 1.});
    }

    std::size_t getNumberOfQubits() const {
        return numberOfQubits;
    }

    BTreeMapKeyDouble const& getData() const {
        return data;
    }

    void getPureState(StringMap<std::complex<double>>& m) const override {
        assert(!data.empty());
        
        if (data.rbegin()->first.ensembleIndex >= 1) {
            // There are cases where, when you trace out the measurement register, this is still a pure state.
            // So this function returns some pure states, but not all.
            return;
        }

        for (auto const& [key, factor]: data) {
            m.add(getStateString(key.basisVector), factor);
        }
    }

    bool operator==(MixedState &other) const {
        if (getNumberOfQubits() != other.getNumberOfQubits()) {
            return false;
        }

        if (toDensityMatrix() != other.toDensityMatrix()) {
            return false;
        }

        auto left = getMeasurementRegisterStatistics();
        auto right = other.getMeasurementRegisterStatistics();

        return areMeasurementRegisterStatisticsEqual(left, right);
    }

    std::size_t getSize() const {
        return data.size();
    }

    void simplify() {
    #ifndef NDEBUG
        std::unique_ptr<MixedState<MaxNumberOfQubits>> mixedStateBefore;
        if (getNumberOfQubits() < config::DEBUG_MAX_QUBITS) {
            mixedStateBefore = std::make_unique<MixedState<MaxNumberOfQubits>>(*this);
        }
    #endif

        MixedStateSimplifier::sparsifyGivens<MaxNumberOfQubits>(data);

    #ifndef NDEBUG
        if (getNumberOfQubits() < config::DEBUG_MAX_QUBITS) {
            assert(*this == *mixedStateBefore);
        }
        assert(isConsistent());
    #endif
    }

    SquareMatrix toDensityMatrix() const { // FIXME: share this with the Python version.
        // Trace out the measurement register.

        EditableSquareMatrix result(1 << getNumberOfQubits());

        EnsembleIndex currentEnsembleIndex = data.begin()->first.ensembleIndex;
        BasisVector currentMeasurementRegister = data.begin()->first.measurementRegister;

        for (auto it = data.begin(); it != data.end(); ++it) {
            if (it->first.measurementRegister != currentMeasurementRegister ||
                it->first.ensembleIndex != currentEnsembleIndex) {

                currentEnsembleIndex = it->first.ensembleIndex;
                currentMeasurementRegister = it->first.measurementRegister;
            }

            result.add(it->first.basisVector.toUInt64(), it->first.basisVector.toUInt64(), std::norm(it->second));

            auto it2 = std::next(it);
            while (it2 != data.end() && it2->first.ensembleIndex == currentEnsembleIndex &&
                   it2->first.measurementRegister == currentMeasurementRegister) {

                result.add(it2->first.basisVector.toUInt64(), it->first.basisVector.toUInt64(), it2->second * std::conj(it->second));
                result.add(it->first.basisVector.toUInt64(), it2->first.basisVector.toUInt64(), it->second * std::conj(it2->second));
                ++it2;
            }
        }

        return result; // slicing
    }

    Matrix toMatrix() {
        MixedStateSimplifier::tidy<MaxNumberOfQubits>(data);

        assert(getNumberOfQubits() <= 32);

        Matrix result(data.rbegin()->first.ensembleIndex.value + 1, 1 << (2 * getNumberOfQubits()));
        
        for (auto const& [key, factor]: data) {
            result.set(
                key.ensembleIndex.value,
                (key.measurementRegister.toUInt64() << getNumberOfQubits()) + key.basisVector.toUInt64(),
                factor);
        }

        return result; // slicing
    }

    void getMeasurementRegisterStatistics(StringMap<double>& m) const override {
        for (auto it = data.rbegin(); it != data.rend(); ++it) {
            m.add(getStateString(it->first.measurementRegister), std::norm(it->second));
        }
    }

    absl::btree_map<std::string, double> getMeasurementRegisterStatistics() const { // FIXME: delete this duplicate
        absl::btree_map<std::string, double> result;

        for (auto it = data.rbegin(); it != data.rend(); ++it) {
            result[getStateString(it->first.measurementRegister)] += std::norm(it->second);
        }

        return result;
    }

    void operator()(CircuitInstruction const& circuitInstruction) override {
        applyCircuitInstruction(circuitInstruction);
        simplify();
    }

    void applyCircuitInstruction(CircuitInstruction const& circuitInstruction) {
        assert(isConsistent());

        auto const operands = circuitInstruction.getOperands();
        auto firstKey = data.begin()->first;
        auto lastKey = data.rbegin()->first;
        auto const maxEnsembleIndex = lastKey.ensembleIndex;
        auto insertionEnsembleIndex = lastKey.ensembleIndex + 1;
        auto currentEnsembleIndex = firstKey.ensembleIndex;
        auto currentMeasurementRegister = firstKey.measurementRegister;

        for (;;) {
            if (firstKey = data.begin()->first; firstKey.ensembleIndex <= maxEnsembleIndex) {
                break;
            }

            // Turn this block into a function?
            if (firstKey.ensembleIndex != currentEnsembleIndex ||
                firstKey.measurementRegister != currentMeasurementRegister) {

                currentMeasurementRegister = firstKey.measurementRegister;
                currentEnsembleIndex = firstKey.ensembleIndex;
                insertionEnsembleIndex.value += applyKrausOperators
                    ? circuitInstruction.getNumberOfKrausOperators()
                    : 1;
            }

            // Turn this block into a function?
            bool applyKrausOperators = firstKey.measurementRegister.test(circuitInstruction.getControlBits());
            if (!applyKrausOperators) {
                auto extracted = data.extract(data.begin());
                extracted.key().ensembleIndex = insertionEnsembleIndex;
                data.insert(data.end(), std::move(extracted));
                continue;
            }

            // Turn this block into a function?
            BasisVector reduced;
            for (std::uint64_t i = 0; i < operands.size(); ++i) {
                auto const& operand = operands[operands.size() - i - 1];
                if (std::holds_alternative<QubitIndex>(operand)) {
                    reduced.set(i, firstKey.basisVector.test(std::get<QubitIndex>(operand).value));
                } else {
                    assert(std::holds_alternative<MeasurementRegisterIndex>(operand));
                    reduced.set(i, firstKey.measurementRegister.test(std::get<MeasurementRegisterIndex>(operand).value));
                }
            }

            for (KrausOperatorIndex operatorIndex; operatorIndex < circuitInstruction.getNumberOfKrausOperators(); ++operatorIndex) {
                for (BasisVector ket = {}; !ket.test(operands.size()); ++ket) {
                    auto const& krausValue = circuitInstruction.getKrausValue(operatorIndex, ket.toUInt64(), reduced.toUInt64());
                    if (krausValue == 0.) { // No epsilon double comparison here, since Kraus matrices are constants.
                        continue;
                    }

                    // At this point we are in a level 4 of nesting (for-for-for-for)
                    // I would try to reduce the complexity of this function,
                    // maybe turning each of this loops into a function itself
                    // That should help also to understand what each loop is doing
                    for (std::size_t operandIndex = 0; operandIndex < operands.size(); ++operandIndex) {
                        auto const& operand = operands[operands.size() - operandIndex - 1];
                        
                        if (std::holds_alternative<QubitIndex>(operand)) {
                            assert(std::get<QubitIndex>(operand).value < getNumberOfQubits());
                            auto modifiedBasisVector = firstKey.basisVector;
                            modifiedBasisVector.set(std::get<QubitIndex>(operand).value, utils::getBit(ket.toUInt64(), operandIndex));
                        } else {
                            assert(std::holds_alternative<MeasurementRegisterIndex>(operand));
                            auto modifiedMeasurementRegister = data.begin()->first.measurementRegister;
                            modifiedMeasurementRegister.set(std::get<MeasurementRegisterIndex>(operand).value, utils::getBit(ket.toUInt64(), operandIndex));
                        }
                    }

                    auto insertionIt = data.try_emplace(data.end(),
                        {
                            .ensembleIndex = insertionEnsembleIndex + operatorIndex.value,
                            .measurementRegister = modifiedMeasurementRegister,
                            .basisVector = modifiedBasisVector
                        },
                        0.
                    );
                    std::complex<double> addedValue = data.begin()->second * krausValue;
                    insertionIt->second += addedValue;
                }
            }

            data.erase(data.begin());
        }

        assert(isConsistent());
    }

private:
    static bool areMeasurementRegisterStatisticsEqual(
        absl::btree_map<std::string, double> left, absl::btree_map<std::string, double> right) {

        // Possible alternative code
        /*
        return (left.size() == right.size()) &&
            std::ranges::all_of(left, [rightIt = right.begin()](auto const &[leftKey, leftValue]) {
                return (leftKey == rightIt->first) &&
                    utils::isNull(leftValue - rightIt++->second);
            });
        */
        //
        if (left.size() != right.size()) {
            return false;
        }

        for (auto leftIt = left.begin(), rightIt = right.begin(); leftIt != left.end(); ++leftit, ++rightIt) {
            if (leftIt->first != rightIt->first ||
                utils::isNotNull(leftIt->second - rightIt->second)) {
                return false;
            }
        }
        return true;
        //
    }

    std::string getStateString(BasisVector s) const {
        auto str = s.toString();

        return str.substr(str.size() - getNumberOfQubits(), str.size());
    }

    bool isConsistent() const;

    // I would define a using alias for this type
    // E.g.
    // using BTreeMapQubitDouble = absl::btree_map<Key<MaxNumberOfQubits>, std::complex<double>>
    absl::btree_map<Key<MaxNumberOfQubits>, std::complex<double>> data;
    std::size_t const numberOfQubits = 1;
};

} // namespace core
} // namespace qx