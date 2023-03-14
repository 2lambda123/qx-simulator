#pragma once

#include <array>
#include <cassert>
#include <climits>
#include <limits>
#include <complex>
#include <string>
#include <span>

#include "qx/CompileTimeConfiguration.hpp"

namespace qx {
namespace utils {

template <typename T>
struct TaggedInteger {
    std::uint64_t value = 0;

    template <typename H> friend H AbslHashValue(H h, TaggedInteger const &taggedInteger) {
        return H::combine(std::move(h), taggedInteger.value);
    }
    
    inline constexpr bool operator==(TaggedInteger const &other) const {
        return value == other.value;
    }

    inline constexpr bool operator<(TaggedInteger const &other) const {
        return value < other.value;
    }

    inline constexpr bool operator>(TaggedInteger const &other) const {
        return value > other.value;
    }

    inline constexpr bool operator<(std::uint64_t integer) const {
        return value < integer;
    }

    inline constexpr bool operator>=(TaggedInteger const &other) const {
        return value >= other.value;
    }

    inline constexpr bool operator<=(TaggedInteger const &other) const {
        return value <= other.value;
    }

    inline constexpr bool operator>=(std::uint64_t integer) const {
        return value >= integer;
    }

    inline constexpr void operator++() {
        ++value;
    }
    
    template <typename S>
    friend TaggedInteger<S> operator*(TaggedInteger<S> tagged, std::uint64_t integer);

    template <typename S>
    friend TaggedInteger<S> operator+(TaggedInteger<S> tagged, std::uint64_t integer);
};

template <typename T>
inline TaggedInteger<T> operator*(TaggedInteger<T> tagged, std::uint64_t integer) {
    return {tagged.value * integer};
}

template <typename T>
inline TaggedInteger<T> operator+(TaggedInteger<T> tagged, std::uint64_t integer) {
    return {tagged.value + integer};
}

} // namespace utils
} // namespace qx