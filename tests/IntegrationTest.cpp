#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"
#include "qx/Random.hpp"
#include "qx/Simulator.hpp"
#include "qx/DefaultOperations.hpp"

namespace qx {

class IntegrationTest {
public:
    SimulationResult runFromString(std::string s,
                                   std::uint64_t iterations = 1) {
        auto output = executeString(s, default_operations::defaultOperations, iterations);
        REQUIRE(std::holds_alternative<SimulationResult>(output));
        return std::get<SimulationResult>(output);
    }

    template <typename T>
    void checkEquality(std::vector<std::pair<std::string, T>> const& left, std::vector<std::pair<std::string, T>> const& right) {
        REQUIRE_EQ(left.size(), right.size());

        for (std::uint64_t i = 0; i < left.size(); ++i) {
            CHECK_EQ(left[i].first, right[i].first);
            CHECK(utils::isNull(left[i].second - right[i].second));
        }
    }
};

bool operator==(std::complex<double> const &left,
                std::complex<double> const &right) {
    return utils::isNull(left - right);
}

TEST_CASE_FIXTURE(IntegrationTest, "Syntax error") {
    auto cqasm = R"(
version 1.0

qubits 2,

x q[0]
)";
    auto result = executeString(cqasm);
    CHECK(std::holds_alternative<SimulationError>(result));
#if defined(_MSC_VER)
    CHECK_EQ(std::get<SimulationError>(result).message, "<unknown>:4:9: syntax error, unexpected ',', expecting $end");
#else
    CHECK_EQ(std::get<SimulationError>(result).message, "<unknown>:4:9: syntax error, unexpected ',', expecting end of file");
#endif
}

TEST_CASE_FIXTURE(IntegrationTest, "Bell pair") {
    auto cqasm = R"(
version 1.0

qubits 2

h q[0]
cnot q[0], q[1]
)";
    auto actual = runFromString(cqasm);

    CHECK_EQ(actual.iterations, 1);
    checkEquality(actual.measurementRegisterStatistics, {{"00", 1.}});
    checkEquality(actual.state, {{"00", 1 / std::sqrt(2)}, {"11", 1 / std::sqrt(2)}});
}

TEST_CASE_FIXTURE(IntegrationTest, "Controlled gates") {
    auto cqasm = R"(
version 1.0

qubits 2

cond(b[0]) h q[0] # not executed
x q[0:1]
measure q[0:1], b[0:1]
cond(b[0]) h q[0]
)";

    auto actual = runFromString(cqasm, 2);

    checkEquality(actual.measurementRegisterStatistics, {{"11", 1.}});
    checkEquality(actual.state, {{"10", 1 / std::sqrt(2)}, {"11", -1 / std::sqrt(2)}});
}

TEST_CASE_FIXTURE(IntegrationTest, "Range operands") {
    auto cqasm = R"(
version 1.0

qubits 6
x q[0:2]
cnot q[0:2], q[3:5]
)";

    auto actual = runFromString(cqasm, 2);

    checkEquality(actual.state, {{"111111", 1.}});
}

TEST_CASE_FIXTURE(IntegrationTest, "Repeated subcircuit") {
    auto cqasm = R"(
version 1.0

qubits 1

.subcircuit(3)
    x q[0]

.another_subcircuit(5)
    h q[0]
)";

    auto actual = runFromString(cqasm, 2);

    checkEquality(actual.state, {{"0", 1 / std::sqrt(2)}, {"1", -1 / std::sqrt(2)}});
}

TEST_CASE_FIXTURE(IntegrationTest, "Classical not gate") {
    auto cqasm = R"(
version 1.0

qubits 1

x q[0]
measure q[0], b[0]
not b[0]
)";
    auto actual = runFromString(cqasm);

    checkEquality(actual.measurementRegisterStatistics, {{"0", 1.}});
    checkEquality(actual.state, {{"1", 1}});
}

TEST_CASE_FIXTURE(IntegrationTest,
                  "Shor-9 error correction code with bit flip") {
    auto cqasm = R"(
version 1.0

qubits 17
prep_z q[0:8]

.encode
CNOT q[0],q[3]
CNOT q[0],q[6]
{ H q[0] | H q[3] | H q[6]}
{ CNOT q[0],q[1] | CNOT q[3],q[4] | CNOT q[6],q[7] }
{ CNOT q[0],q[2] | CNOT q[3],q[5] | CNOT q[6],q[8] }

.error
X q[2]

.stabilizers
prep_z q[9:16]
H q[9:16]

CZ q[8],q[16]
CZ q[7],q[16]
CZ q[7],q[15]
CZ q[6],q[15]
CZ q[5],q[14]
CZ q[4],q[14]
CZ q[4],q[13]
CZ q[3],q[13]
CZ q[2],q[12]
CZ q[1],q[12]
CZ q[1],q[11]
CZ q[0],q[11]

CNOT q[10],q[8]
CNOT q[10],q[7]
CNOT q[10],q[6]
CNOT q[10],q[5]
CNOT q[10],q[4]
CNOT q[10],q[3]

CNOT q[9],q[5]
CNOT q[9],q[4]
CNOT q[9],q[3]
CNOT q[9],q[2]
CNOT q[9],q[1]
CNOT q[9],q[0]

H q[9:16]

measure_z q[9:16], b[9:16]
)";

    std::size_t iterations = 7;
    auto actual = runFromString(cqasm, iterations);

    checkEquality(actual.measurementRegisterStatistics, {{"00001000000000000", 1.}});

    checkEquality(actual.state, {
                 {"00001000000000011", 1 / std::sqrt(8)},
                 {"00001000000000100", 1 / std::sqrt(8)},
                 {"00001000000111011", 1 / std::sqrt(8)},
                 {"00001000000111100", 1 / std::sqrt(8)},
                 {"00001000111000011", 1 / std::sqrt(8)},
                 {"00001000111000100", 1 / std::sqrt(8)},
                 {"00001000111111011", 1 / std::sqrt(8)},
                 {"00001000111111100", 1 / std::sqrt(8)}});
}

TEST_CASE_FIXTURE(IntegrationTest, "Probabilistic error models are not supported") {
    auto cqasm = R"(
version 1.0

qubits 1

error_model depolarizing_channel, 0.1
)";

    CHECK(std::holds_alternative<SimulationError>(executeString(cqasm)));
}

} // namespace qx