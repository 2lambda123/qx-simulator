#pragma once

#include "qx/Core.h"

namespace qx {
namespace gates {

template <std::size_t N>
using UnitaryMatrix = core::DenseUnitaryMatrix<N>;

using namespace std::complex_literals;

static constexpr long double PI = 3.141592653589793238462643383279502884L;
static constexpr double SQRT_2 = std::sqrt(2);

static constexpr UnitaryMatrix<2> IDENTITY = UnitaryMatrix<2>::identity();

static constexpr UnitaryMatrix<2> X({{{0, 1}, {1, 0}}});

static constexpr UnitaryMatrix<2> Y({{{0, -1i}, {1i, 0}}});

static constexpr UnitaryMatrix<2> Z({{{1, 0}, {0, -1}}});

static constexpr UnitaryMatrix<2> S({{{1, 0}, {0, 1i}}});
static_assert(S * S == Z);

static constexpr UnitaryMatrix<2> SDAG = S.dagger();

static constexpr UnitaryMatrix<2> T({{{1, 0}, {0, 1 / SQRT_2 + 1i / SQRT_2}}});
static_assert(T * T == S);

static constexpr UnitaryMatrix<2> TDAG = T.dagger();

UnitaryMatrix<2> constexpr RX(double theta) {
    return UnitaryMatrix<2>({{{std::cos(theta / 2), -1i * std::sin(theta / 2)}, {-1i * std::sin(theta / 2), std::cos(theta / 2)}}});
}

constexpr auto X90 = RX(PI / 2);
constexpr auto MX90 = RX(-PI / 2);

UnitaryMatrix<2> constexpr RY(double theta) {
    return UnitaryMatrix<2>({{{std::cos(theta / 2), -std::sin(theta / 2)}, {std::sin(theta / 2), std::cos(theta / 2)}}});
}

constexpr auto Y90 = RY(PI / 2);
constexpr auto MY90 = RY(-PI / 2);

UnitaryMatrix<2> constexpr RZ(double theta) {
    return UnitaryMatrix<2>({{{std::cos(theta / 2) - 1i * std::sin(theta / 2), 0}, {0, std::cos(theta / 2) + 1i * std::sin(theta / 2)}}});
}

static constexpr UnitaryMatrix<2> H({{{1 / SQRT_2, 1 / SQRT_2}, {1 / SQRT_2, -1 / SQRT_2}}});

static constexpr UnitaryMatrix<4> CNOT({{{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 0, 1}, {0, 0, 1, 0}}});

static constexpr UnitaryMatrix<4> SWAP({{{1, 0, 0, 0}, {0, 0, 1, 0}, {0, 1, 0, 0}, {0, 0, 0, 1}}});

static constexpr UnitaryMatrix<4> CZ({{{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, -1}}});

static constexpr UnitaryMatrix<8> TOFFOLI(
    {{{1, 0, 0, 0, 0, 0, 0, 0},
      {0, 1, 0, 0, 0, 0, 0, 0},
      {0, 0, 1, 0, 0, 0, 0, 0},
      {0, 0, 0, 1, 0, 0, 0, 0},
      {0, 0, 0, 0, 1, 0, 0, 0},
      {0, 0, 0, 0, 0, 1, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 1},
      {0, 0, 0, 0, 0, 0, 1, 0}}});

// TODO: make this even more user-friently with ctrl operator and multiplication by double.

}
}