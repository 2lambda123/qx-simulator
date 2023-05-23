#pragma once

#include <array>
#include <cassert>
#include <climits>
#include <limits>
#include <string>
#include <span>

#include "qx/utils/StrongTypes.hpp"

namespace qx {
namespace utils {

inline bool getBit(std::uint64_t x, std::uint64_t index) {
    return (x >> index) & 1;
}

inline void setBit(std::uint64_t &x, std::uint64_t index, bool value) {
    assert(index < 64);
    x = (x & ~(1UL << index)) | (static_cast<std::uint64_t>(value) << index);
}

// std::bitset is slightly slower than this.

class BitReference {
public:
    BitReference(std::uint64_t& v, std::uint64_t i) : value(v), index(i) {}

    operator bool() const { return getBit(value, index); }

    void operator=(bool bit) { setBit(value, index, bit); }

private:
    std::uint64_t& value;
    std::uint64_t index;
};

template <std::uint64_t NumberOfBits> class BasisVector {
private:
    static constexpr std::uint64_t BITS_PER_UNIT =
        CHAR_BIT * sizeof(std::uint64_t);
    static_assert(BITS_PER_UNIT == 64);
    static constexpr std::uint64_t STORAGE_SIZE =
        NumberOfBits / BITS_PER_UNIT +
        (NumberOfBits % BITS_PER_UNIT >= 1 ? 1 : 0);

public:
    static constexpr std::uint64_t getNumberOfBits() { return NumberOfBits; }

    BasisVector() = default;

    BasisVector(std::string_view s) {
        assert(s.size() <= NumberOfBits);
        std::uint64_t i = 0;
        for (auto it = s.rbegin(); it != s.rend(); ++it, ++i) {
            set(i, *it == '1' ? true : false);
        }
    }

    BasisVector(std::uint64_t s) {
        static_assert(NumberOfBits <= BITS_PER_UNIT);
        assert(NumberOfBits == BITS_PER_UNIT || (s >> NumberOfBits) == 0);
        data[0] = s;
    }

    inline void reset() { data = {}; }

    inline bool test(std::uint64_t index) const {
        assert(index < NumberOfBits && "BasisVector::test bit index out of range");
        return getBit(data[index / BITS_PER_UNIT], index % BITS_PER_UNIT);
    }

    inline bool test(std::span<MeasurementRegisterIndex const> indices) const {
        return std::ranges::all_of(indices, [](auto i){ return test(i.value); });
    }

    inline void set(std::uint64_t index, bool value = true) {
        assert(index < NumberOfBits && "BasisVector::set bit index out of range");
        setBit(data[index / BITS_PER_UNIT], index % BITS_PER_UNIT, value);
    }

    inline bool operator==(BasisVector<NumberOfBits> const &other) const {
        return data == other.data;
    }

    inline bool operator<(BasisVector<NumberOfBits> const &other) const {
        return LT(other);
    }

    inline bool operator<=(BasisVector<NumberOfBits> const &other) const {
        return *this == other || *this < other;
    }

    inline void operator^=(BasisVector<NumberOfBits> const &other) {
        for (std::uint64_t i = 0; i < data.size(); ++i) {
            data[i] ^= other.data[i];
        }
    }

    inline bool operator[](QubitIndex index) const {
        return test(index.value);
    }

    inline BitReference operator[](QubitIndex index) {
        return BitReference(data[index.value / BITS_PER_UNIT], index.value % BITS_PER_UNIT);
    }

    inline void operator++() {
        for (auto& d: data) {
            if (d < std::numeric_limits<std::uint64_t>::max()) {
                ++d;
                return;
            }
            d = 0;
        }
    }

    template <typename H> friend H AbslHashValue(H h, BasisVector const &basisVector) {
        return H::combine(std::move(h), basisVector.data);
    }

    std::uint64_t toUInt64() const {
#ifndef NDEBUG
        for (std::uint64_t i = 1; i < data.size(); ++i) {
            assert(data[i] == 0);
        }
#endif

        return data[0];
    }

    std::string toString() const {
        // Possible alternative code
        std::string result(NumberOfBits, '0');
        for (std::uint64_t i = 0; i < NumberOfBits; ++i) {
            if (test(NumberOfBits - i - 1)) {
                result[i] =  '1';
            }
        }
        // Or even
        /*
        std::ranges::for_each(std::ranges::iota(0, NumberOfBits), [&result](auto i) {
            if (test(NumberOfBits - i - 1)) {
                result[i] =  '1';
            }
        });
        */
        return result;
    }

private:
    template <std::uint64_t StorageSize = STORAGE_SIZE>
    inline bool LT(BasisVector<NumberOfBits> const &other) const {
        constexpr auto Index = StorageSize - 1;
        if constexpr (Index > 0) {
            if (data[Index] == other.data[Index]) {
                return LT<StorageSize - 1>(other);
            }
        }

        return data[Index] < other.data[Index];
    }

    std::array<std::uint64_t, STORAGE_SIZE> data{};
};

// The following might not be true on all systems.
#if !defined(_MSC_VER)
static_assert(BasisVector<32>::STORAGE_SIZE == 1);
static_assert(BasisVector<64>::STORAGE_SIZE == 1);
static_assert(BasisVector<65>::STORAGE_SIZE == 2);
#endif

} // namespace utils
} // namespace qx