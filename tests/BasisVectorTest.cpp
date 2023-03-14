#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "qx/utils/BasisVector.hpp"
#include "doctest/doctest.h"
// #include "absl/hash/hash_testing.h"

namespace qx {
namespace utils {

class BasisVectorTest {
public:
};

TEST_CASE_FIXTURE(BasisVectorTest, "Set/test") {
    BasisVector<15> victim{};
    CHECK(!victim.test(1));
    CHECK(!victim.test(13));

    CHECK_EQ(victim.toUInt64(), 0);
    victim.set(0);
    CHECK(!victim.test(1));
    CHECK_EQ(victim.toUInt64(), 1);
    victim.set(1);
    CHECK_EQ(victim.toUInt64(), 3);
    CHECK(victim.test(1));
    CHECK(victim.test(0));
    CHECK(!victim.test(10));

    victim.set(3);
    CHECK(victim.test(3));
    CHECK_EQ(victim.toUInt64(), 11);
}

TEST_CASE_FIXTURE(BasisVectorTest, "Set/test 64 bits") {
    BasisVector<64> victim{};

    victim.set(31);

    CHECK(victim.test(31));
    for (auto i = 0; i < 64; ++i) {
        if (i != 31) {
            CHECK(!victim.test(i));
        }
    }
}

TEST_CASE_FIXTURE(BasisVectorTest, "Set/test with a lot of bits") {
    BasisVector<150> victim{};
    CHECK(!victim.test(120));
    CHECK(!victim.test(130));

    victim.set(0);
    CHECK(victim.test(0));
    CHECK(!victim.test(1));
    victim.set(132);
    CHECK(victim.test(0));
    CHECK(victim.test(132));
    CHECK(!victim.test(149));

    CHECK_EQ(victim.toString(),
             "00000000000000000100000000000000000000000000000000000000000000000"
             "00000000000000000000000000000000000000000000000000000000000000000"
             "00000000000000000001");
}

TEST_CASE_FIXTURE(BasisVectorTest, "From string") {
    BasisVector<5> victim{"1010"};
    CHECK(!victim.test(0));
    CHECK(victim.test(1));
    CHECK(!victim.test(2));
    CHECK(victim.test(3));
    CHECK(!victim.test(4));
    CHECK_EQ(victim.toUInt64(), 10);
}

TEST_CASE_FIXTURE(BasisVectorTest, "toString") {
    BasisVector<15> victim{};
    victim.set(0);
    CHECK_EQ(victim.toString(), "000000000000001");

    victim.set(13);
    CHECK_EQ(victim.toString(), "010000000000001");
}

// Requires GMock
// TEST_CASE_FIXTURE(BasisVectorTest, "Hash") {
//     BasisVector<15> victim1{};
//     CHECK(absl::VerifyTypeImplementsAbslHashCorrectly({victim1}));

//     BasisVector<15456> victim2{};
//     victim2.set(4542);
//     victim2.set(8945);
//     CHECK(absl::VerifyTypeImplementsAbslHashCorrectly({victim2}));
// }

TEST_CASE_FIXTURE(BasisVectorTest, "operator^") {
    SUBCASE("Small bitset") {
        BasisVector<15> victim{"000010000010001"};
        BasisVector<15> mask{"000011001000001"};

        victim ^= mask;
        CHECK_EQ(victim.toString(), "000001001010000");
    }

    SUBCASE("Large bitset") {
        BasisVector<123456> victim{};
        victim.set(457);

        BasisVector<123456> mask{};

        mask.set(457);
        mask.set(654);

        victim ^= mask;

        CHECK(!victim.test(457));
        CHECK(victim.test(654));
        CHECK(!victim.test(875));
    }
}

} // namespace utils
} // namespace qx