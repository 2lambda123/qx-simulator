#include "qx/Circuit.h"

#include "qx/Random.h"
#include <algorithm>

namespace qx {
namespace {
struct InstructionExecutor {
public:
    InstructionExecutor(core::QuantumState &s) : quantumState(s){};

    void operator()(Circuit::Measure const &m) {
        quantumState.measure(m.qubitIndex, [this]() { return rand.next(); });
    }

    void operator()(Circuit::MeasureAll const &) {
        quantumState.measureAll([this]() { return rand.next(); });
    }

    void operator()(Circuit::PrepZ const &r) {
        quantumState.prep(r.qubitIndex, [this]() { return rand.next(); });
    }
    
    void operator()(Circuit::MeasurementRegisterOperation const &op) {
        op.operation(quantumState.getMeasurementRegister());
    }

    template <std::size_t N> void operator()(Circuit::Unitary<N> const &u) {
        quantumState.apply(u.matrix, u.operands);
    }

private:
    core::QuantumState &quantumState;
    random::UniformRandomNumberGenerator rand{0.0, 1.0};
};
} // namespace

void Circuit::execute(core::QuantumState &quantumState) const {
    std::size_t it = iterations;
    InstructionExecutor instructionExecutor(quantumState);
    while (it-- > 0) {
        for (auto const& controlledInstruction : controlledInstructions) {
            auto const& controlBits = controlledInstruction.controlBits;
            if (controlBits) {
                auto measurementRegister = quantumState.getMeasurementRegister();
                auto isBitNotSet = [&measurementRegister](auto const& cb) { return !measurementRegister.test(cb.value); };
                if (std::any_of(controlBits->begin(), controlBits->end(), isBitNotSet)) {
                        continue;
                }
            }

            auto const& instruction = controlledInstruction.instruction;

            // AppleClang doesn't support std::visit
            // std::visit(instructionExecutor, instruction);
            if (auto *measure = std::get_if<Circuit::Measure>(&instruction)) {
                instructionExecutor(*measure);
            } else if (auto *measureAll =
                           std::get_if<Circuit::MeasureAll>(&instruction)) {
                instructionExecutor(*measureAll);
            } else if (auto *prepZ =
                           std::get_if<Circuit::PrepZ>(&instruction)) {
                instructionExecutor(*prepZ);
            } else if (auto *classicalOp =
                           std::get_if<Circuit::MeasurementRegisterOperation>(
                               &instruction)) {
                instructionExecutor(*classicalOp);
            } else if (auto *instruction1 =
                           std::get_if<Circuit::Unitary<1>>(&instruction)) {
                instructionExecutor(*instruction1);
            } else if (auto *instruction2 =
                           std::get_if<Circuit::Unitary<2>>(&instruction)) {
                instructionExecutor(*instruction2);
            } else if (auto *instruction3 =
                           std::get_if<Circuit::Unitary<3>>(&instruction)) {
                instructionExecutor(*instruction3);
            } else {
                assert(false && "Unimplemented circuit instruction");
            }
        }
    }
}

} // namespace qx