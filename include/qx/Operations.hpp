#pragma once

#include <variant>
#include <string_view>
#include "absl/container/inlined_vector.h"
#include "absl/container/flat_hash_map.h"
#include "qx/Matrix.hpp"
#include "qx/utils/BasisVector.hpp"
#include "qx/utils/StrongTypes.hpp"

namespace qx {
class Operations {
public:
    enum class OperandType {
        Double,
        Int,
        Qubit,
        ClassicalBit
    };

    // The fmt library is very convenient for these things
    // It would turn signatureToString implementation into a single liner
    # include <fmt/format.h>
    auto format_as(OperandType operandType) {
        switch (operandType) {
            case OperandType::Qubit: return "Qubit";
            case OperandType::ClassicalBit: return "ClassicalBit";
            case OperandType::Double: return "Double";
            case OperandType::Int: return "Int";
        }
        return "";
    }

    using Signature = absl::InlinedVector<OperandType, config::MAX_INLINED_OPERAND_TYPES>;
    using KrausOperators = absl::InlinedVector<SquareMatrix, config::MAX_INLINED_KRAUS_OPERATORS>;
    using StaticOperand = std::variant<double, std::int64_t>;
    using StaticOperands = absl::InlinedVector<StaticOperand, 2>;
    using KrausOperatorsGenerationFunction = std::function<KrausOperators(StaticOperands const&)>;

    static std::string signatureToString(Signature const& signature) {
        // This is basically saying:
        // - I want to return a string (fmt::format)
        // - The format of that string is: open parens, something, close parens
        // - That something is a vector of OperandType separated by ", "
        // - You know how to format an OperandType because I'm providing a format_as(OperandType) function above
        return fmt::format("({})", fmt::join(signature, ", "));
        /*
        std::stringstream res;
        res << "(";
        bool first = true;
        for (auto t: signature) {
            if (!first) {
                res << ", ";
            } else {
                first = false;
            }

            switch (t) {
                case OperandType::Qubit:
                    res << "Qubit";
                    break;
                case OperandType::ClassicalBit:
                    res << "ClassicalBit";
                    break;
                case OperandType::Double:
                    res << "Double";
                    break;
                case OperandType::Int:
                    res << "Int";
                    break;
            }
        }

        res << ")";
        return res.str();
        */
    }

    static void checkValidKrausOperatorSet(std::string const& name, std::int64_t numberOfDynamicOperands, KrausOperators const& krausOperators) {
        if (numberOfDynamicOperands <= 0) {
            throw std::runtime_error(std::string("Invalid number of dynamic operands for operation ") + name);
        }
        
        if (krausOperators.empty()) {
            throw std::runtime_error(std::string("Kraus operators set for operation ") + name + " is empty");
        }

        // All Kraus operator matrices need to have the same size.
        for (auto const& m: krausOperators) {
            if (m.getSize() != krausOperators.begin()->getSize()) {
                throw std::runtime_error(std::string("Kraus operators for operation ") + name + " don't all have the same size");
            }
        }

        // Size of Kraus operators need to be a power of two.
        if (!std::has_single_bit(krausOperators.begin()->getSize())) {
            throw std::runtime_error(std::string("Kraus operators size for operation ") + name + " is not a power of two");
        }

        if (static_cast<std::int64_t>(std::countr_zero(krausOperators.begin()->getSize())) != numberOfDynamicOperands) {
            throw std::runtime_error(
                std::string("Kraus operators for operation ") + name +
                " have size " + std::to_string(krausOperators.begin()->getSize()) +
                ", which do not match the number of operands (" + std::to_string(numberOfDynamicOperands) + ")");
        }

        // Kraus operators need to satisfy the completeness relation.
        // That is: all the eigenvalues of the sum of E_k^t * E_k need to <= 1.
        EditableSquareMatrix accumulator(krausOperators.begin()->getSize());
        for (auto const& m: krausOperators) {
            accumulator += m.dagger() * m;
        }

        double spectralRadius = computeSpectralRadius(accumulator);

        if (spectralRadius > 1. + config::EPS) {
            throw std::runtime_error(std::string("Kraus operators for operation ") + name + " are not non-trace-increasing");
        }
    }

    KrausOperators get(std::string const& name, Signature signature, StaticOperands const& staticOperands) const {
#ifndef NDEBUG
        assert(signature.size() >= staticOperands.size());
        auto signatureIt = signature.begin();
        auto operandIt = staticOperands.begin();

        while (signatureIt != signature.end()) {
            if (*signatureIt == OperandType::Qubit || *signatureIt == OperandType::ClassicalBit) {
                ++signatureIt;
                continue;
            }

            assert(
                (operandIt != staticOperands.end()) &&
                ((*signatureIt == OperandType::Double && std::holds_alternative<double>(*operandIt)) ||
                 (*signatureIt == OperandType::Int && std::holds_alternative<std::int64_t>(*operandIt)))
            );
            
            ++operandIt;
            ++signatureIt;
        }
#endif

        auto it = krausOperatorsGeneration.find(name);

        if (it == krausOperatorsGeneration.end()) {
            throw std::runtime_error("No registered operation with name " + name);
        }

        if (it->second.first != signature) {
            throw std::runtime_error(
                "Requested operation " + name +
                " with signature " + Operations::signatureToString(signature) +
                ", but the registered operation with that name has signature " + Operations::signatureToString(it->second.first));
        }
        
        auto result = it->second.second(staticOperands);

        return result;
    }

    void add(std::string const& name, Signature signature, UnitaryMatrix const& m) {
        add(name, std::move(signature), KrausOperators{m});
    }

    void add(std::string const& name, Signature signature, KrausOperators const& ks) {
        checkValidKrausOperatorSet(name, signature.size(), ks);
        
        assert(signature.size() == static_cast<std::uint64_t>(std::countr_zero(ks.begin()->getSize())));
        assert(std::ranges::all_of(signature, [](auto x) {
            return x == OperandType::Qubit || x == OperandType::ClassicalBit;
        }));

        add(name, signature, [ks](Operations::StaticOperands const& staticOperands) -> KrausOperators {
            assert(staticOperands.empty());
            return ks;
        });
    }

    using SingleFloatOperand = KrausOperators (*) (double); // Here return type templated on # of qubit/bit operands would be beneficial...
    void add(std::string const& name, Signature signature, SingleFloatOperand f) {
#ifndef NDEBUG
        bool seenDouble = false;
        for (auto t: signature) {
            if (t == OperandType::Double) {
                assert(!seenDouble);
                seenDouble = true;
            } else {
                assert(t == OperandType::Qubit || t == OperandType::ClassicalBit);
            }
        }
        assert(seenDouble && "No double in signature");
#endif

        add(name, signature, [f, signature](Operations::StaticOperands const& staticOperands) -> KrausOperators {
            assert (staticOperands.size() == 1 && std::holds_alternative<double>(staticOperands[0]));

            auto res = f(std::get<double>(staticOperands[0]));

            return res;
        });
    }

    // Possible alternative code
    //
    // void add(std::string const& name, Signature signature, KrausOperatorsGenerationFunction&& f) {
    //     auto blah{ name, std::make_pair(signature, std::forward<KrausOperatorsGenerationFunction>(f)) };
    //     const auto& inserted = krausOperatorsGeneration.insert(std::move(blah));
    //     if (!inserted.second)
    void add(std::string const& name, Signature signature, KrausOperatorsGenerationFunction f) {
        auto inserted = krausOperatorsGeneration.insert({name, std::make_pair(signature, std::move(f))}).second;

        if (!inserted) {
            throw std::runtime_error("Registering operation with name " + name + " twice");
        }
    }

private:
    // I would define a using alias for this type
    // E.g.
    // using KrausOperatorGeneration = absl::flat_hash_map<std::string, std::pair<Signature, KrausOperatorsGenerationFunction>>
    KrausOperatorsGeneration krausOperatorsGeneration;
};
}