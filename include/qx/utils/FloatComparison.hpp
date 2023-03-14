#pragma once

#include <complex>

#include "qx/CompileTimeConfiguration.hpp"

namespace qx {
namespace utils {

inline constexpr bool isNotNull(double d) {
#if defined(_MSC_VER)
    return d > config::EPS || -d > config::EPS;
#else
    return std::abs(d) > config::EPS;
#endif
}

inline constexpr bool isNotNull(std::complex<double> c) {
    return isNotNull(c.real()) || isNotNull(c.imag());
}

template <typename T>
inline constexpr bool isNull(T t) {
    return !isNotNull(t);
}

} // namespace utils
} // namespace qx