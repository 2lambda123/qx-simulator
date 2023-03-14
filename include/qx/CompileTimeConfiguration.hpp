#pragma once

namespace qx {
namespace config {

// Absolute TOLerance for double comparison
static constexpr double ATOL = 1e-12;

// Number of decimals in output
static constexpr std::uint64_t const OUTPUT_DECIMALS = 8;

static constexpr std::uint64_t DEBUG_MAX_QUBITS = 5;

static constexpr std::uint64_t MAX_INLINED_OPERANDS = 5;

static constexpr std::uint64_t MAX_INLINED_CONTROL_BITS = 2;

static constexpr std::uint64_t MAX_INLINED_MATRIX_ELEMENTS = 16;

static constexpr std::uint64_t MAX_INLINED_KRAUS_OPERATORS = 2;

static constexpr std::uint64_t MAX_INLINED_OPERAND_TYPES = 6;

} // namespace config

} // namespace qx