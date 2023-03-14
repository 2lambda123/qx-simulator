#pragma once

#include "qx/Core.hpp"
#include "qx/Operations.hpp"

namespace qx {
namespace default_operations {

using namespace std::complex_literals;

inline long double PI = 3.141592653589793238462643383279502884L;
inline double SQRT_2 = 1.414213562373095048801688724209698078L;

inline UnitaryMatrix IDENTITY = UnitaryMatrix::identity(2);

inline UnitaryMatrix X{{0, 1},
                       {1, 0}};

inline UnitaryMatrix Y{{0, -1i},
                       {1i, 0}};

inline UnitaryMatrix Z{{1, 0},
                       {0, -1}};

inline UnitaryMatrix S{{1, 0},
                       {0, 1i}};

inline UnitaryMatrix SDAG = S.dagger();

inline UnitaryMatrix T{{1, 0},
                       {0, 1 / SQRT_2 + 1i / SQRT_2}};

inline UnitaryMatrix TDAG = T.dagger();

inline Operations::KrausOperators RX(double theta) {
    return {UnitaryMatrix{{std::cos(theta / 2), -1i * std::sin(theta / 2)},
                          {-1i * std::sin(theta / 2), std::cos(theta / 2)}}};
}

inline auto X90 = RX(PI / 2);
inline auto MX90 = RX(-PI / 2);

inline Operations::KrausOperators RY(double theta) {
    return {UnitaryMatrix{{std::cos(theta / 2), -std::sin(theta / 2)},
                          {std::sin(theta / 2), std::cos(theta / 2)}}};
}

inline auto Y90 = RY(PI / 2);
inline auto MY90 = RY(-PI / 2);

inline Operations::KrausOperators RZ(double theta) {
    return {UnitaryMatrix{{std::cos(theta / 2) - 1i * std::sin(theta / 2), 0},
                          {0, std::cos(theta / 2) + 1i * std::sin(theta / 2)}}};
}

inline auto Z90 = RZ(PI / 2);
inline auto MZ90 = RZ(-PI / 2);

inline UnitaryMatrix H{{1 / SQRT_2, 1 / SQRT_2},
                       {1 / SQRT_2, -1 / SQRT_2}};

inline UnitaryMatrix CNOT{
    {1, 0, 0, 0},
    {0, 1, 0, 0},
    {0, 0, 0, 1},
    {0, 0, 1, 0}};

inline UnitaryMatrix SWAP{
    {1, 0, 0, 0},
    {0, 0, 1, 0},
    {0, 1, 0, 0},
    {0, 0, 0, 1}};

inline UnitaryMatrix CZ{
    {1, 0, 0, 0},
    {0, 1, 0, 0},
    {0, 0, 1, 0},
    {0, 0, 0, -1}};

inline Operations::KrausOperators CR(double theta) {
    return {UnitaryMatrix{{1, 0, 0, 0},
                          {0, 1, 0, 0},
                          {0, 0, 1, 0},
                          {0, 0, 0, std::cos(theta) + 1i * std::sin(theta)}}};
}

inline UnitaryMatrix TOFFOLI{
    {1, 0, 0, 0, 0, 0, 0, 0},
    {0, 1, 0, 0, 0, 0, 0, 0},
    {0, 0, 1, 0, 0, 0, 0, 0},
    {0, 0, 0, 1, 0, 0, 0, 0},
    {0, 0, 0, 0, 1, 0, 0, 0},
    {0, 0, 0, 0, 0, 1, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 1},
    {0, 0, 0, 0, 0, 0, 1, 0}};

inline std::initializer_list<SquareMatrix> MEAS_Z_NO_CLASSICAL_EFFECT{
    SquareMatrix{{1, 0},
                 {0, 0}},
    SquareMatrix{{0, 0},
                 {0, 1}}};

inline std::initializer_list<SquareMatrix> MEAS_Z{
    SquareMatrix{{1, 0, 0, 0},
                 {0, 0, 0, 0},
                 {0, 0, 0, 0},
                 {0, 0, 0, 0}},
    SquareMatrix{{0, 1, 0, 0},
                 {0, 0, 0, 0},
                 {0, 0, 0, 0},
                 {0, 0, 0, 0}},
    SquareMatrix{{0, 0, 0, 0},
                 {0, 0, 0, 0},
                 {0, 0, 0, 0},
                 {0, 0, 1, 0}},
    SquareMatrix{{0, 0, 0, 0},
                 {0, 0, 0, 0},
                 {0, 0, 0, 0},
                 {0, 0, 0, 1}}};

// FIXME
inline std::initializer_list<SquareMatrix> MEAS_X{
    SquareMatrix{{1 / 2., 1 / 2.},
                 {1 / 2., 1 / 2.}},
    SquareMatrix{{1 / 2., -1 / 2.},
                 {-1 / 2., 1 / 2.}}};
inline std::initializer_list<SquareMatrix> MEAS_Y{
    SquareMatrix{{1 / 2., 1i / 2.},
                 {1i / 2., -1 / 2.}},
    SquareMatrix{{1 / 2., -1i / 2.},
                 {-1i / 2., 1 / 2.}}};
inline std::initializer_list<SquareMatrix> PREP_Z{SquareMatrix{{1, 0},
                                                               {0, 0}},
                                                  SquareMatrix{{0, 1},
                                                               {0, 0}}};
inline std::initializer_list<SquareMatrix> PREP_Z_SWAP{
    SquareMatrix{{1, 1},
                 {0, 0}}};
inline std::initializer_list<SquareMatrix> PREP_Y{
    SquareMatrix{{1 / 2., 1i / 2.},
                 {1i / 2., -1 / 2.}},
    SquareMatrix{{1 / 2., -1i / 2.},
                 {1i / 2., -1 / 2.}}};

inline Operations::KrausOperators DEPOLARIZING_CHANNEL(double lambda) {
    auto k0 = std::sqrt(1. - 3 * lambda / 4) * SquareMatrix::identity(2);
    auto k1 = std::sqrt(lambda / 4) * X;
    auto k2 = std::sqrt(lambda / 4) * Y;
    auto k3 = std::sqrt(lambda / 4) * Z;

    return {k0, k1, k2, k3};
}

inline Operations::KrausOperators AMPLITUDE_DAMPING(double gamma) {
    auto e0 = SquareMatrix{{1, 0},
                           {0, std::sqrt(1 - gamma)}};
    auto e1 = SquareMatrix{{0, std::sqrt(gamma)},
                           {0, 0}};

    return {e0, e1};
}

inline Operations::KrausOperators PHASE_DAMPING(double lambda) {
    auto e0 = SquareMatrix{{1, 0},
                           {0, std::sqrt(1 - lambda)}};
    auto e1 = SquareMatrix{{0, 0},
                           {0, std::sqrt(lambda)}};

    return {e0, e1};
}

inline Operations createDefaultOperations() {
    auto Q = OperandType::Qubit;
    auto B = OperandType::ClassicalBit;
    auto D = OperandType::Double;

    Operations defaultOperations;

    defaultOperations.add("id", {Q}, IDENTITY);
    defaultOperations.add("x", {Q}, X);
    defaultOperations.add("not", {B}, X);
    defaultOperations.add("x90", {Q}, X90);
    defaultOperations.add("mx90", {Q}, MX90);
    defaultOperations.add("y", {Q}, Y);
    defaultOperations.add("y90", {Q}, Y90);
    defaultOperations.add("my90", {Q}, MY90);
    defaultOperations.add("z", {Q}, Z);
    defaultOperations.add("z90", {Q}, Z90);
    defaultOperations.add("mz90", {Q}, MZ90);
    defaultOperations.add("s", {Q}, S);
    defaultOperations.add("sdag", {Q}, SDAG);
    defaultOperations.add("t", {Q}, T);
    defaultOperations.add("tdag", {Q}, TDAG);
    defaultOperations.add("h", {Q}, H);
    defaultOperations.add("cnot", {Q, Q}, CNOT);
    defaultOperations.add("swap", {Q, Q}, SWAP);
    defaultOperations.add("cz", {Q, Q}, CZ);
    defaultOperations.add("toffoli", {Q, Q, Q}, TOFFOLI);
    defaultOperations.add("measure", {Q, B}, MEAS_Z);
    defaultOperations.add("measure_z", {Q, B}, MEAS_Z);
    // defaultOperations.add("measure_x", {Q, B}, MEAS_X); // FIXME
    // defaultOperations.add("measure_y", {Q, B}, MEAS_Y);
    defaultOperations.add("prep", {Q}, PREP_Z);
    defaultOperations.add("prep_z", {Q}, PREP_Z);
    defaultOperations.add("rx", {Q, D}, RX);
    defaultOperations.add("ry", {Q, D}, RY);
    defaultOperations.add("rz", {Q, D}, RZ);
    defaultOperations.add("cr", {Q, Q, D}, CR);
    defaultOperations.add("depolarizing_channel", {Q, D}, DEPOLARIZING_CHANNEL);
    defaultOperations.add("phase_damping", {Q, D}, PHASE_DAMPING);
    defaultOperations.add("amplitude_damping", {Q, D}, AMPLITUDE_DAMPING);

    return defaultOperations;
}

inline Operations defaultOperations = createDefaultOperations();

} // namespace default_operations
} // namespace qx