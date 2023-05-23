#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"
#include "qx/Core.hpp"
#include "qx/DefaultOperations.hpp"

namespace qx {
namespace core {

using namespace std::complex_literals;

class MixedStateTest {
protected:
    static void checkEqComplex(std::complex<double> left, std::complex<double> right) {
        CHECK_EQ(left.real(), ::doctest::Approx(right.real()));
        CHECK_EQ(left.imag(), ::doctest::Approx(right.imag()));
    }

    static void checkMatrix(MixedState<>& victim,
                            std::initializer_list<std::initializer_list<std::complex<double>>> expected) {
        auto actual = victim.toMatrix();
        auto expectedMatrix = Matrix(expected);

        CHECK_EQ(actual, expectedMatrix);
    }
};

TEST_CASE_FIXTURE(MixedStateTest, "Identity") {
    MixedState victim(1);

    SquareMatrix expected{{1, 0},
                          {0, 0}};

    CHECK_EQ(victim.toDensityMatrix(), expected);

    UnitaryMatrix identity{{1, 0},
                           {0, 1}};

    victim.applyCircuitInstruction(CircuitInstruction({identity}, {QubitIndex{0}}));

    CHECK_EQ(victim.toDensityMatrix(), expected);
    
    auto measurementRegisterStatistics = victim.getMeasurementRegisterStatistics();

    CHECK_EQ(measurementRegisterStatistics["0"], 1.);
}

TEST_CASE_FIXTURE(MixedStateTest, "Unitary gate and pure state") {
    MixedState victim(1);

    double angle = 1.23;
    auto rx = default_operations::RX(angle);

    victim(CircuitInstruction(rx, {QubitIndex{0}}));

    AssociationVectorStringMap<std::complex<double>>::T state;
    AssociationVectorStringMap wrappedState(state);
    victim.getPureState(wrappedState);

    REQUIRE_EQ(state.size(), 2);
    CHECK_EQ(state[0].first, "0");
    CHECK_EQ(state[0].second, std::cos(angle / 2));
    CHECK_EQ(state[1].first, "1");
    CHECK_EQ(state[1].second, - 1i * std::sin(angle / 2));
}

TEST_CASE_FIXTURE(MixedStateTest, "Measure in computational basis") {
    MixedState victim(1);

    victim.applyCircuitInstruction(CircuitInstruction(default_operations::MEAS_Z, {QubitIndex{0}, MeasurementRegisterIndex{0}}));

    SquareMatrix expected{{1, 0},
                          {0, 0}};

    CHECK_EQ(victim.toDensityMatrix(), expected);

    checkMatrix(victim, {{1., 0., 0., 0.}});

    auto measurementRegisterStatistics = victim.getMeasurementRegisterStatistics();

    CHECK_EQ(measurementRegisterStatistics["0"], 1.);
}

TEST_CASE_FIXTURE(MixedStateTest, "Measure in computational basis after Hadamard") {
    MixedState victim(1);

    victim.applyCircuitInstruction(CircuitInstruction({default_operations::H}, {QubitIndex{0}}));
    victim.applyCircuitInstruction(CircuitInstruction(default_operations::MEAS_Z, {QubitIndex{0}, MeasurementRegisterIndex{0}}));

    SquareMatrix expected{{1./2, 0.},
                          {0., 1./2}};

    CHECK_EQ(victim.toDensityMatrix(), expected);

    checkMatrix(victim, {{1/std::sqrt(2), 0., 0., 0.}, {0., 0., 0., 1/std::sqrt(2)}});

    AssociationVectorStringMap<std::complex<double>>::T state;
    AssociationVectorStringMap wrappedState(state);
    victim.getPureState(wrappedState);

    CHECK_EQ(state.size(), 0);

    auto measurementRegisterStatistics = victim.getMeasurementRegisterStatistics();

    CHECK_EQ(measurementRegisterStatistics["0"], ::doctest::Approx(.5));
    CHECK_EQ(measurementRegisterStatistics["1"], ::doctest::Approx(.5));
}

TEST_CASE_FIXTURE(MixedStateTest, "Probabilistic bit-shift") {
    MixedState victim(1);

    double p = 0.2;
    SquareMatrix k0 = std::sqrt(1 - p) * SquareMatrix::identity(2);

    SquareMatrix k1 = std::sqrt(p) * default_operations::X;

    victim.applyCircuitInstruction(CircuitInstruction({k0, k1}, {QubitIndex{0}}));
    
    SquareMatrix expected{{1-p, 0},
                          {0,   p}};

    CHECK_EQ(victim.toDensityMatrix(), expected);
}

TEST_CASE_FIXTURE(MixedStateTest, "|+> state") {
    MixedState victim(1);

    victim.applyCircuitInstruction(CircuitInstruction({default_operations::H}, {QubitIndex{0}}));
    
    SquareMatrix expected{{1/2., 1/2.},
                          {1/2., 1/2.}};

    CHECK_EQ(victim.toDensityMatrix(), expected);
}

TEST_CASE_FIXTURE(MixedStateTest, "|+> state with 2 qubits") {
    MixedState victim(2);

    victim.applyCircuitInstruction(CircuitInstruction({default_operations::H}, {QubitIndex{1}}));
    
    SquareMatrix expected{{1/2., 0, 1/2., 0},
                          {0, 0, 0, 0},
                          {1/2., 0, 1/2., 0},
                          {0, 0, 0, 0}};

    CHECK_EQ(victim.toDensityMatrix(), expected);
}

TEST_CASE_FIXTURE(MixedStateTest, "Measure Bell pair separately") {
    MixedState victim(3);

    victim.applyCircuitInstruction(CircuitInstruction({default_operations::H}, {QubitIndex{1}}));
    
    victim.applyCircuitInstruction(CircuitInstruction({default_operations::CNOT}, {QubitIndex{1}, QubitIndex{2}}));

    victim.applyCircuitInstruction(CircuitInstruction(default_operations::MEAS_Z, {QubitIndex{2}, MeasurementRegisterIndex{2}}));

    auto measurementRegisterStatistics = victim.getMeasurementRegisterStatistics();
    CHECK_EQ(measurementRegisterStatistics["000"], ::doctest::Approx(.5));
    CHECK_EQ(measurementRegisterStatistics["100"], ::doctest::Approx(.5));

    victim.applyCircuitInstruction(CircuitInstruction(default_operations::MEAS_Z, {QubitIndex{0}, MeasurementRegisterIndex{0}}));
    victim.applyCircuitInstruction(CircuitInstruction(default_operations::MEAS_Z, {QubitIndex{1}, MeasurementRegisterIndex{1}}));

    measurementRegisterStatistics = victim.getMeasurementRegisterStatistics();
    CHECK_EQ(measurementRegisterStatistics["000"], ::doctest::Approx(.5));
    CHECK_EQ(measurementRegisterStatistics["110"], ::doctest::Approx(.5));
}

TEST_CASE_FIXTURE(MixedStateTest, "Prep") {
    MixedState victim(1);

    victim.applyCircuitInstruction(CircuitInstruction({default_operations::H}, {QubitIndex{0}}));
    
    victim.applyCircuitInstruction(CircuitInstruction(default_operations::MEAS_Z, {QubitIndex{0}, MeasurementRegisterIndex{0}}));

    victim.applyCircuitInstruction(CircuitInstruction(default_operations::PREP_Z, {QubitIndex{0}}));

    auto measurementRegisterStatistics = victim.getMeasurementRegisterStatistics();
    CHECK_EQ(measurementRegisterStatistics["0"], ::doctest::Approx(.5));
    CHECK_EQ(measurementRegisterStatistics["1"], ::doctest::Approx(.5));

    checkMatrix(victim, {{1/std::sqrt(2), 0, 0, 0}, {0, 0, 1/std::sqrt(2), 0}});
}

TEST_CASE_FIXTURE(MixedStateTest, "Measure mixed state") {
    MixedState victim(3);

    double p = 0.5;
    SquareMatrix k0 = std::sqrt(1 - p) * SquareMatrix::identity(2);

    SquareMatrix k1 = std::sqrt(p) * default_operations::H;

    victim.applyCircuitInstruction(CircuitInstruction({k0, k1}, {QubitIndex{1}}));

    victim.applyCircuitInstruction(CircuitInstruction(default_operations::MEAS_Z, {QubitIndex{1}, MeasurementRegisterIndex{1}}));

    auto measurementRegisterStatistics = victim.getMeasurementRegisterStatistics();
    CHECK_EQ(measurementRegisterStatistics["000"], ::doctest::Approx(.75));
    CHECK_EQ(measurementRegisterStatistics["010"], ::doctest::Approx(.25));
}

TEST_CASE_FIXTURE(MixedStateTest, "Cond gate not applied") {
    MixedState victim(1);

    victim.applyCircuitInstruction(CircuitInstruction({default_operations::X}, {QubitIndex{0}}, {MeasurementRegisterIndex{0}}));

    SquareMatrix expected{{1., 0.},
                          {0., 0.}};

    CHECK_EQ(victim.toDensityMatrix(), expected);

    checkMatrix(victim, {{1., 0., 0., 0.}});

    auto measurementRegisterStatistics = victim.getMeasurementRegisterStatistics();

    CHECK_EQ(measurementRegisterStatistics["0"], ::doctest::Approx(1.));
}

TEST_CASE_FIXTURE(MixedStateTest, "Prep with measure + cond X") {
    MixedState victim(1);

    victim.applyCircuitInstruction(CircuitInstruction({default_operations::H}, {QubitIndex{0}}));

    victim.applyCircuitInstruction(CircuitInstruction(default_operations::MEAS_Z, {QubitIndex{0}, MeasurementRegisterIndex{0}}));

    victim.applyCircuitInstruction(CircuitInstruction({default_operations::X}, {QubitIndex{0}}, {MeasurementRegisterIndex{0}}));

    SquareMatrix expected{{1., 0.},
                          {0., 0.}};

    CHECK_EQ(victim.toDensityMatrix(), expected);

    checkMatrix(victim, {{1/std::sqrt(2), 0., 0., 0.}, {0., 0., 1/std::sqrt(2), 0.}});

    auto measurementRegisterStatistics = victim.getMeasurementRegisterStatistics();

    CHECK_EQ(measurementRegisterStatistics["0"], ::doctest::Approx(.5));
    CHECK_EQ(measurementRegisterStatistics["1"], ::doctest::Approx(.5));
}

TEST_CASE_FIXTURE(MixedStateTest, "Conditional X") {
    MixedState victim(2);

    victim.applyCircuitInstruction(CircuitInstruction({default_operations::X}, {MeasurementRegisterIndex{0}})); // Bit-flip of classical measurement bit...

    victim.applyCircuitInstruction(CircuitInstruction({default_operations::X}, {QubitIndex{1}}, {MeasurementRegisterIndex{0}}));

    CHECK_EQ(victim.getSize(), 1);
    checkMatrix(victim, {{0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.}});

    auto measurementRegisterStatistics = victim.getMeasurementRegisterStatistics();

    CHECK_EQ(measurementRegisterStatistics["01"], ::doctest::Approx(1.));
}

TEST_CASE_FIXTURE(MixedStateTest, "Simplify pure state" * doctest::skip()) {
    MixedState victim(2);

    victim.applyCircuitInstruction(CircuitInstruction({default_operations::H}, {QubitIndex{0}}));
    victim.applyCircuitInstruction(CircuitInstruction({default_operations::CNOT}, {QubitIndex{0}, QubitIndex{1}}));
    
    checkMatrix(victim, {{1/std::sqrt(2), 0., 0., 1/std::sqrt(2), 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}});

    victim.simplify();

    checkMatrix(victim, {{1/std::sqrt(2), 0., 0., 1/std::sqrt(2), 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}});

    auto measurementRegisterStatistics = victim.getMeasurementRegisterStatistics();
    CHECK_EQ(measurementRegisterStatistics["00"], ::doctest::Approx(1.));
}

TEST_CASE_FIXTURE(MixedStateTest, "Simplify mixed state" * doctest::skip()) {
    MixedState victim(1);

    double p = 0.5;
    SquareMatrix k0 = std::sqrt(1 - p) * SquareMatrix::identity(2);
    SquareMatrix k1 = std::sqrt(p) * default_operations::H;
    victim.applyCircuitInstruction(CircuitInstruction({k0, k1}, {QubitIndex{0}}));

    checkMatrix(victim, {{1/std::sqrt(2), 0, 0, 0}, {1/2., 1/2., 0, 0}});

    victim.applyCircuitInstruction(CircuitInstruction({k0, k1}, {QubitIndex{0}}));
    
    SquareMatrix expectedDensityMatrix{{3/4., 1/4.},
                                       {1/4., 1/4.}};

    CHECK_EQ(victim.toDensityMatrix(), expectedDensityMatrix);

    checkMatrix(victim, {{1/2., 0, 0, 0}, {0.5/std::sqrt(2), 0.5/std::sqrt(2), 0, 0}, {0.5/std::sqrt(2), 0.5/std::sqrt(2), 0, 0}, {1/2., 0, 0, 0}});

    victim.simplify();

    checkMatrix(victim, {{-0.86602540378443871, -0.28867513459481292, 0, 0}, {0, -0.40824829046386313, 0, 0}});
    CHECK_EQ(victim.toDensityMatrix(), expectedDensityMatrix);
}

TEST_CASE_FIXTURE(MixedStateTest, "Simplify very mixed state" * doctest::skip()) {
    MixedState victim(2);

    double p = 0.5;
    SquareMatrix k0 = std::sqrt(1 - p) * SquareMatrix::identity(2);
    SquareMatrix k1 = std::sqrt(p) * default_operations::H;
    victim.applyCircuitInstruction(CircuitInstruction({k0, k1}, {QubitIndex{0}}));
    victim.applyCircuitInstruction(CircuitInstruction({k0, k1}, {QubitIndex{1}}));
    victim.applyCircuitInstruction(CircuitInstruction({k0, k1}, {QubitIndex{1}}));
    victim.applyCircuitInstruction(CircuitInstruction({k0, k1}, {QubitIndex{0}}));
    victim.applyCircuitInstruction(CircuitInstruction({k0, k1}, {QubitIndex{0}}));
    victim.applyCircuitInstruction(CircuitInstruction({k0, k1}, {QubitIndex{1}}));
    victim.applyCircuitInstruction(CircuitInstruction({k0, k1}, {QubitIndex{0}}));
    victim.applyCircuitInstruction(CircuitInstruction({k0, k1}, {QubitIndex{1}}));
    victim.applyCircuitInstruction(CircuitInstruction({k0, k1}, {QubitIndex{0}}));
    
    CHECK_EQ(victim.getSize(), 1953);

    victim.simplify();

    CHECK_EQ(victim.getSize(), 9);
    checkMatrix(victim, {{-3/4., -1/4., -1/4., -1/12., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                         {0., -.5/std::sqrt(2), 0, -1./(6*std::sqrt(2)), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                         {0., 0., -.5/std::sqrt(2), -1./(6*std::sqrt(2)), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                         {0., 0., 0., -1/6., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}});
}

TEST_CASE_FIXTURE(MixedStateTest, "Simplify simplifies rho = .5 |+><+| + .5 |-><-| to rho = .5 |0><0| + .5 |1><1|" * doctest::skip()) {
    MixedState victim(3);

    SquareMatrix k0 = std::sqrt(.5) * SquareMatrix::identity(2);
    SquareMatrix k1 = std::sqrt(.5) * default_operations::X;
    
    victim.applyCircuitInstruction(CircuitInstruction({default_operations::CNOT}, {QubitIndex{0}, QubitIndex{1}}));
    victim.applyCircuitInstruction(CircuitInstruction({k0, k1}, {QubitIndex{2}}));
    victim.applyCircuitInstruction(CircuitInstruction({default_operations::H}, {QubitIndex{2}}));
    
    auto actualMatrix = victim.toMatrix();
    CHECK_EQ(actualMatrix.getNumberOfRows(), 2);
    checkEqComplex(actualMatrix.get(0, 0), .5);
    checkEqComplex(actualMatrix.get(0, 4), .5);
    checkEqComplex(actualMatrix.get(1, 0), .5);
    checkEqComplex(actualMatrix.get(1, 4), -.5);

    victim.simplify();

    auto newActualMatrix = victim.toMatrix();
    CHECK_EQ(newActualMatrix.getNumberOfRows(), 2);
    checkEqComplex(newActualMatrix.get(0, 0), -1/std::sqrt(2));
    checkEqComplex(newActualMatrix.get(1, 4), -1/std::sqrt(2));
}

} // namespace core
} // namespace qx