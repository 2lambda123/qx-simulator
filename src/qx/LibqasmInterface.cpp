#include "qx/LibqasmInterface.hpp"

#include "cqasm-v1-semantic-gen.hpp"
#include "qx/Circuit.hpp"
#include "qx/Core.hpp"
#include "qx/Operations.hpp"
#include "absl/container/inlined_vector.h"

namespace qx {

namespace cq = ::cqasm::v1::semantic;
namespace values = ::cqasm::v1::values;

namespace {
Operations::Signature getSignature(cq::Instruction const& instruction) {
    Operations::Signature res;
    for (auto const& op: instruction.operands) {
        if (op->as_const_real()) {
            res.push_back(Operations::OperandType::Double);
        } else if (op->as_const_int()) {
            res.push_back(Operations::OperandType::Int);
        } else if (op->as_qubit_refs()) {
            res.push_back(Operations::OperandType::Qubit);
        } else if (op->as_bit_refs()) {
            res.push_back(Operations::OperandType::ClassicalBit);
        } else {
            throw std::runtime_error("Unhandled operand type");
        }
    }
    return res;

}

absl::InlinedVector<CircuitInstruction::DynamicOperandsVector, config::MAX_INLINED_OPERANDS> getDynamicOperands(cq::Instruction const& instruction) {
    absl::InlinedVector<CircuitInstruction::DynamicOperandsVector, config::MAX_INLINED_OPERANDS> res;
    for (auto const& op: instruction.operands) {
        if (auto qubitRefs = op->as_qubit_refs()) {
            assert(res.size() == 0 || res.size() == qubitRefs->index.size());
            res.resize(std::max(res.size(), qubitRefs->index.size()));
            for (std::uint64_t i = 0; i < qubitRefs->index.size(); ++i) {
                res[i].push_back(QubitIndex{ static_cast<std::uint64_t>(qubitRefs->index[i]->value) });
            }
        } else if (auto bitRefs = op->as_bit_refs()) {
            assert(res.size() == 0 || res.size() == bitRefs->index.size());
            res.resize(std::max(res.size(), bitRefs->index.size()));
            for (std::uint64_t i = 0; i < bitRefs->index.size(); ++i) {
                res[i].push_back(MeasurementRegisterIndex{ static_cast<std::uint64_t>(bitRefs->index[i]->value) });
            }
        }
    }
    return res;
}

Operations::StaticOperands getStaticOperands(cq::Instruction const& instruction) {
    Operations::StaticOperands res;
    for (auto const& op: instruction.operands) {
        if (auto real = op->as_const_real()) {
            res.push_back(real->value);
        } else if (auto integer = op->as_const_int()) {
            res.push_back(integer->value);
        } else if (op->as_qubit_refs()) {
            // No-op
        } else if (op->as_bit_refs()) {
            // No-op
        } else {
            throw std::runtime_error("Unhandled operand type");
        }
    }
    return res;
}

std::string to_string(cq::NodeType nodeType) {
    switch (nodeType) {
    case cq::NodeType::AnnotationData:
        return "AnnotationData";
    case cq::NodeType::Block:
        return "Block";
    case cq::NodeType::BreakStatement:
        return "BreakStatement";
    case cq::NodeType::Bundle:
        return "Bundle";
    case cq::NodeType::BundleExt:
        return "BundleExt";
    case cq::NodeType::ContinueStatement:
        return "ContinueStatement";
    case cq::NodeType::ErrorModel:
        return "ErrorModel";
    case cq::NodeType::ForLoop:
        return "ForLoop";
    case cq::NodeType::ForeachLoop:
        return "ForeachLoop";
    case cq::NodeType::GotoInstruction:
        return "GotoInstruction";
    case cq::NodeType::IfElse:
        return "IfElse";
    case cq::NodeType::IfElseBranch:
        return "IfElseBranch";
    case cq::NodeType::Instruction:
        return "Instruction";
    case cq::NodeType::Mapping:
        return "Mapping";
    case cq::NodeType::Program:
        return "Program";
    case cq::NodeType::RepeatUntilLoop:
        return "RepeatUntilLoop";
    case cq::NodeType::SetInstruction:
        return "SetInstruction";
    case cq::NodeType::Subcircuit:
        return "Subcircuit";
    case cq::NodeType::Variable:
        return "Variable";
    case cq::NodeType::Version:
        return "Version";
    case cq::NodeType::WhileLoop:
        return "WhileLoop";
    }

    return "Unknown";
}

class GateConvertor : public cq::RecursiveVisitor {
public:
    GateConvertor(qx::Circuit &c, Operations const& ops) : circuit(c), operations(ops) {}

    void visit_instruction(cq::Instruction &instr) override { addQuantumOperation(instr); }

    void visit_bundle_ext(cq::BundleExt &node) override {
        node.items.visit(*this);
    }

    void visit_node(cq::Node &node) override {
        throw std::runtime_error("Statements of the following type are not "
                                 "supported by the simulator: " +
                                 to_string(node.type()));
    }

private:
    void addQuantumOperation(const cq::Instruction &instruction) {
        auto b = instruction.condition->as_const_bool();
        if (b && !(b->value)) {
            return;
        }

        auto name = instruction.name;
        std::transform(name.begin(), name.end(), name.begin(), [](auto c){ return std::tolower(c); });

        auto signature = getSignature(instruction);

        Operations::StaticOperands nonQubitOperands = getStaticOperands(instruction);

        auto krausOperators = operations.get(name, signature, nonQubitOperands);

        Operations::checkValidKrausOperatorSet(instruction.name, signature.size() - nonQubitOperands.size(), krausOperators);

        absl::InlinedVector<MeasurementRegisterIndex, config::MAX_INLINED_CONTROL_BITS> controlBits;

        if (auto bitref = instruction.condition->as_bit_refs()) {
            for (auto const &b : bitref->index) {
                controlBits.push_back(
                    MeasurementRegisterIndex{static_cast<std::size_t>(b->value)});
            }
        }

        auto dynamicOperands = getDynamicOperands(instruction);

        for (auto ops: dynamicOperands) {
            circuit.addInstruction(CircuitInstruction(krausOperators, ops, controlBits));
        }
    }

    qx::Circuit &circuit;
    Operations const& operations;
};
} // namespace

qx::Circuit loadCqasmCode(cq::Subcircuit const &subcircuit, Operations const& operations, std::uint64_t qubitCount) {
    qx::Circuit circuit(qubitCount, subcircuit.name, subcircuit.iterations);

    for (const auto &statement : subcircuit.body->statements) {
        GateConvertor gateConvertor(circuit, operations);

        statement->visit(gateConvertor);
    }

    return circuit;
}

} // namespace qx