#pragma once

#include "absl/container/btree_map.h"
#include "absl/container/inlined_vector.h"
#include <cassert>
#include <complex>
#include <limits>
#include <ranges>
#include <algorithm>

#include "qx/CircuitInstruction.hpp"
#include "qx/Common.hpp"
#include "qx/CompileTimeConfiguration.hpp"
#include "qx/MixedStateSimplifier.hpp"
#include "qx/SimulationResult.hpp"
#include "qx/utils/FloatComparison.hpp"
#include "qx/utils/StrongTypes.hpp"

namespace qx {
namespace core {

template <typename V> class StringMap {
public:
    virtual void add(std::string const &key, V value) = 0;
};

template <typename V> class AssociationVectorStringMap : public StringMap<V> {
public:
    using T = std::vector<std::pair<std::string, V>>;

    explicit AssociationVectorStringMap(T &v) : data(v) { assert(data.empty()); }

    void add(std::string const &key, V value) override {
        // Could be sorted.
        auto it =
            std::find_if(data.rbegin(), data.rend(),
                         [&key](auto const &kv) { return kv.first == key; });

        if (it == data.rend()) {
            data.push_back(std::make_pair(key, value));
        } else {
            it->second += value;
        }
    }

private:
    T &data;
};

template <std::uint64_t MaxNumberOfQubits> class MixedState;

class MixedStateBase {
public:
    virtual void operator()(CircuitInstruction const &circuitInstruction) = 0;
    virtual void getPureState(StringMap<std::complex<double>> &m) const = 0;
    virtual void
    getMeasurementRegisterStatistics(StringMap<double> &m) const = 0;

    virtual ~MixedStateBase() = default;
};

template <std::uint64_t MaxNumberOfQubits = 64>
class MixedState final : public MixedStateBase {
public:
    using UnderlyingT = absl::btree_map<Key<MaxNumberOfQubits>, std::complex<double>>;
    using BasisVector = utils::BasisVector<MaxNumberOfQubits>;
    using DynamicOperands =
        std::span<CircuitInstruction::DynamicOperand const>;

    static BasisVector getReducedBasisVector(Key<MaxNumberOfQubits> const& key, DynamicOperands const& operands) {
        BasisVector reduced;
        for (std::uint64_t i = 0; i < operands.size(); ++i) {
            auto const &operand = operands[operands.size() - i - 1];
            if (std::holds_alternative<QubitIndex>(operand)) {
                reduced.set(i, key.basisVector.test(
                                    std::get<QubitIndex>(operand).value));
            } else {
                assert(std::holds_alternative<MeasurementRegisterIndex>(
                    operand));
                reduced.set(
                    i,
                    key.measurementRegister.test(
                        std::get<MeasurementRegisterIndex>(operand).value));
            }
        }

        return reduced;
    }

    static std::pair<BasisVector, BasisVector> expand(Key<MaxNumberOfQubits> const& key, DynamicOperands const& operands, BasisVector const& toReplace) {
        auto modifiedBasisVector = key.basisVector;
        auto modifiedMeasurementRegister =
            key.measurementRegister;

        for (std::uint64_t operandIndex = 0;
                operandIndex < operands.size(); ++operandIndex) {
            auto const &operand =
                operands[operands.size() - operandIndex - 1];

            if (std::holds_alternative<QubitIndex>(operand)) {
                modifiedBasisVector.set(
                    std::get<QubitIndex>(operand).value,
                    toReplace.test(operandIndex));
            } else {
                assert(std::holds_alternative<
                        MeasurementRegisterIndex>(operand));
                modifiedMeasurementRegister.set(
                    std::get<MeasurementRegisterIndex>(operand)
                        .value,
                    toReplace.test(operandIndex));
            }
        }

        return std::make_pair(modifiedMeasurementRegister, modifiedBasisVector);
    }

    explicit MixedState(std::uint64_t n) : numberOfQubits(n) {
        data.insert({Key<MaxNumberOfQubits>{.ensembleIndex = {},
                                            .measurementRegister = {},
                                            .basisVector = {}},
                     1.});
    }

    std::uint64_t getNumberOfQubits() const { return numberOfQubits; }

    UnderlyingT const &
    getData() const {
        return data;
    }

    void getPureState(StringMap<std::complex<double>> &m) const override {
        assert(!data.empty());

        if (data.rbegin()->first.ensembleIndex >= EnsembleIndex{1}) {
            // There are cases where, when you trace out the measurement
            // register, this is still a pure state. So this function returns
            // some pure states, but not all.
            return;
        }

        for (auto const &[key, factor] : data) {
            m.add(getStateString(key.basisVector), factor);
        }
    }

    bool operator==(MixedState const &other) const {
        if (getNumberOfQubits() != other.getNumberOfQubits()) {
            return false;
        }

        if (toDensityMatrix() != other.toDensityMatrix()) {
            return false;
        }

        auto leftMeas = getMeasurementRegisterStatistics();
        auto rightMeas = other.getMeasurementRegisterStatistics();

        return std::ranges::equal(leftMeas, rightMeas, [](auto l, auto r) {
            return l.first == r.first && utils::isNull(l.second - r.second);
        });
    }

    std::uint64_t getSize() const { return data.size(); }

    void simplify() {
#ifndef NDEBUG
        std::unique_ptr<MixedState<MaxNumberOfQubits>> mixedStateBefore;
        if (getNumberOfQubits() < config::DEBUG_MAX_QUBITS) {
            mixedStateBefore =
                std::make_unique<MixedState<MaxNumberOfQubits>>(*this);
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

    SquareMatrix
    toDensityMatrix() const { // FIXME: share this with the Python version.
        // Trace out the measurement register.

        EditableSquareMatrix result(1 << getNumberOfQubits());

        EnsembleIndex currentEnsembleIndex = data.begin()->first.ensembleIndex;
        BasisVector currentMeasurementRegister =
            data.begin()->first.measurementRegister;

        for (auto it = data.begin(); it != data.end(); ++it) {
            if (it->first.measurementRegister != currentMeasurementRegister ||
                it->first.ensembleIndex != currentEnsembleIndex) {
                currentEnsembleIndex = it->first.ensembleIndex;
                currentMeasurementRegister = it->first.measurementRegister;
            }

            result.add(it->first.basisVector.toUInt64(),
                       it->first.basisVector.toUInt64(), std::norm(it->second));

            auto it2 = std::next(it);
            while (it2 != data.end() &&
                   it2->first.ensembleIndex == currentEnsembleIndex &&
                   it2->first.measurementRegister ==
                       currentMeasurementRegister) {
                result.add(it2->first.basisVector.toUInt64(),
                           it->first.basisVector.toUInt64(),
                           it2->second * std::conj(it->second));
                result.add(it->first.basisVector.toUInt64(),
                           it2->first.basisVector.toUInt64(),
                           it->second * std::conj(it2->second));
                ++it2;
            }
        }

        return result; // slicing
    }

    Matrix toMatrix() {
        MixedStateSimplifier::tidy<MaxNumberOfQubits>(data);

        assert(getNumberOfQubits() <= 32);

        EditableMatrix result(data.rbegin()->first.ensembleIndex.value + 1,
                              1 << (2 * getNumberOfQubits()));

        for (auto const &[key, factor] : data) {
            result.set(
                key.ensembleIndex.value,
                (key.measurementRegister.toUInt64() << getNumberOfQubits()) +
                    key.basisVector.toUInt64(),
                factor);
        }

        return result; // slicing
    }

    void getMeasurementRegisterStatistics(StringMap<double> &m) const override {
        for (auto it = data.rbegin(); it != data.rend(); ++it) {
            auto v = std::norm(it->second);
            if (utils::isNotNull(v)) {
                m.add(getStateString(it->first.measurementRegister), v);
            }
        }
    }

    absl::btree_map<std::string, double>
    getMeasurementRegisterStatistics() const { // FIXME: delete this duplicate
        absl::btree_map<std::string, double> result;

        for (auto it = data.rbegin(); it != data.rend(); ++it) {
            auto v = std::norm(it->second);
            if (utils::isNotNull(v)) {
                result[getStateString(it->first.measurementRegister)] += v;
            }
        }

        return result;
    }

    void operator()(CircuitInstruction const &circuitInstruction) override {
        applyCircuitInstruction(circuitInstruction);
        simplify();
    }

    void applyCircuitInstruction(CircuitInstruction const &circuitInstruction) {
        assert(isConsistent());

        auto const operands = circuitInstruction.getOperands();

        auto const maxEnsembleIndex = data.rbegin()->first.ensembleIndex;
        auto insertionEnsembleIndex = data.rbegin()->first.ensembleIndex + 1;

        EnsembleIndex currentEnsembleIndex = data.begin()->first.ensembleIndex;
        auto currentMeasurementRegister =
            data.begin()->first.measurementRegister;
        while (data.begin()->first.ensembleIndex <= maxEnsembleIndex) {
            bool applyKrausOperators =
                data.begin()->first.measurementRegister.test(
                    circuitInstruction.getControlBits());

            if (data.begin()->first.ensembleIndex != currentEnsembleIndex ||
                data.begin()->first.measurementRegister !=
                    currentMeasurementRegister) {
                currentMeasurementRegister =
                    data.begin()->first.measurementRegister;
                currentEnsembleIndex = data.begin()->first.ensembleIndex;
                insertionEnsembleIndex.value +=
                    applyKrausOperators
                        ? circuitInstruction.getNumberOfKrausOperators()
                        : 1;
            }

            if (!applyKrausOperators) {
                auto extracted = data.extract(data.begin());
                extracted.key().ensembleIndex = insertionEnsembleIndex;
                data.insert(data.end(), std::move(extracted));
                continue;
            }

            auto reduced = getReducedBasisVector(data.begin()->first, operands);

            for (KrausOperatorIndex operatorIndex;
                 operatorIndex < KrausOperatorIndex{circuitInstruction.getNumberOfKrausOperators()};
                 ++operatorIndex) {
                for (BasisVector ket = {}; !ket.test(operands.size()); ++ket) {
                    auto const &krausValue = circuitInstruction.getKrausValue(
                        operatorIndex, ket.toUInt64(), reduced.toUInt64());

                    if (krausValue == 0.) {
                        // No epsilon double comparison here, since Kraus
                        // matrices are constants.
                        continue;
                    }

                    std::complex<double> addedValue =
                        data.begin()->second * krausValue;

                    auto [newMeasurementRegister, newBasisVector] = expand(data.begin()->first, operands, ket);

                    auto insertionIt = data.try_emplace(
                        data.end(),
                        {.ensembleIndex =
                             insertionEnsembleIndex + operatorIndex.value,
                         .measurementRegister = newMeasurementRegister,
                         .basisVector = newBasisVector},
                        0.);
                    insertionIt->second += addedValue;
                }
            }

            data.erase(data.begin());
        }

        assert(isConsistent());
    }

private:
    std::string getStateString(BasisVector s) const {
        auto str = s.toString();

        return str.substr(str.size() - getNumberOfQubits(), str.size());
    }

    bool isConsistent() const;

    UnderlyingT data;
    std::uint64_t const numberOfQubits = 1;
};

} // namespace core
} // namespace qx