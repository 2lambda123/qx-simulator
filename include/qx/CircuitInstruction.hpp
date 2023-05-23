#include <span>
#include <bit>
#include <cmath>

#include "qx/Matrix.hpp"
#include "qx/Common.hpp"
#include "qx/CompileTimeConfiguration.hpp"
#include "qx/Operations.hpp"

namespace qx {

class CircuitInstruction {
public:
    using DynamicOperand = std::variant<QubitIndex, MeasurementRegisterIndex>;
    using DynamicOperandsVector = absl::InlinedVector<DynamicOperand, config::MAX_INLINED_OPERANDS>;
    using ControlBitsVector = absl::InlinedVector<MeasurementRegisterIndex, config::MAX_INLINED_CONTROL_BITS>;
    using KrausOperators = Operations::KrausOperators;

    CircuitInstruction(KrausOperators ks, DynamicOperandsVector o, ControlBitsVector b = {})
        : krausOperators(std::move(ks)), dynamicOperands(std::move(o)), controlBits(std::move(b)) {
#ifndef NDEBUG
        Operations::checkValidKrausOperatorSet("unknown", dynamicOperands.size(), krausOperators);
#endif
    }

    std::uint64_t getNumberOfKrausOperators() const {
        return krausOperators.size();
    }

    std::complex<double> getKrausValue(KrausOperatorIndex operatorIndex, std::uint64_t i, std::uint64_t j) const {
        assert(operatorIndex.value < getNumberOfKrausOperators());
        assert(i < krausOperators[operatorIndex.value].getNumberOfRows());
        assert(j < krausOperators[operatorIndex.value].getNumberOfCols());
        return krausOperators[operatorIndex.value].get(i, j);
    }

    std::span<DynamicOperand const> getOperands() const {
        return {dynamicOperands.begin(), dynamicOperands.end()};
    }

    std::span<MeasurementRegisterIndex const> getControlBits() const {
        return {controlBits.begin(), controlBits.end()};
    }

private:
    KrausOperators const krausOperators;
    DynamicOperandsVector dynamicOperands;
    ControlBitsVector controlBits;
};

}