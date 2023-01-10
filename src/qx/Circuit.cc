#include "qx/Circuit.h"

#include "qx/Random.h"

namespace qx {
namespace {
struct InstructionExecutor {
public:
    InstructionExecutor(core::QuantumState& s) : quantumState(s) {};

    void operator() (Circuit::Measure const& m) { quantumState.measure(m.qubitIndex, [this]() { return rand.next(); }); }
    void operator() (Circuit::MeasureAll const&) { quantumState.measureAll([this]() { return rand.next(); }); }
    void operator() (Circuit::PrepZ const& r) { quantumState.prep(r.qubitIndex); }

    template <std::size_t N>
    void operator() (Circuit::Unitary<N>& u) {
        if (u.controlBits) {
            auto measurementRegister = quantumState.getMeasurementRegister();
            for (auto b: *u.controlBits) {
                if (!measurementRegister.test(b.value)) {
                    return;
                }
            }
        }

        quantumState.apply(u.matrix, u.operands);
    }

private:
    core::QuantumState& quantumState;
    random::uniform_random_number_generator rand{0.0, 1.0};
};
}

void Circuit::execute(core::QuantumState& quantumState) const {
    std::size_t it = iterations;
    InstructionExecutor instructionExecutor(quantumState);
    while (it-- > 0) {
        for (auto instruction: instructions) {
            std::visit(instructionExecutor, instruction);
        }
    }
}

}