#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "qx/MixedStateSimplifier.hpp"
#include "doctest/doctest.h"

namespace qx {

using namespace std::complex_literals;

class MixedStateSimplifierTest {
public:
    void checkTriangularization(AbstractMatrix const &input) {
        auto householderAccumulator =
            EditableSquareMatrix::identity(input.getNumberOfRows());

        EditableMatrix m(input);

        qrTriangularize(m, &householderAccumulator);

        CHECK_EQ(input, householderAccumulator.dagger() * m);
    }

    bool areEquals(absl::btree_map<Key<>, std::complex<double>> const &left,
                   absl::btree_map<Key<>, std::complex<double>> const &right) {
        auto leftIt = left.begin();
        auto rightIt = right.begin();

        while (leftIt != left.end() || rightIt != right.end()) {
            while (leftIt != left.end() && utils::isNull(leftIt->second)) {
                ++leftIt;
            }

            while (rightIt != right.end() && utils::isNull(rightIt->second)) {
                ++rightIt;
            }

            if ((leftIt != left.end() && rightIt == right.end()) ||
                (leftIt == left.end() && rightIt != right.end())) {
                CAPTURE("Sizes not equal");
                return false;
            }

            if (leftIt != left.end() && rightIt != right.end() &&
                (leftIt->first != rightIt->first ||
                 utils::isNotNull(leftIt->second - rightIt->second))) {
                CAPTURE(leftIt->first);
                CAPTURE(rightIt->first);
                return false;
            }

            ++leftIt;
            ++rightIt;
        }

        return true;
    }
};

TEST_CASE_FIXTURE(MixedStateSimplifierTest, "Identity Householder matrix") {
    EditableMatrix m{{1, 1}, {0, 1}};

    HouseholderMatrix householder{Matrix{{1}, {0}}};

    UnitaryMatrix expectedHouseholder{{-1, 0}, {0, 1}};

    CHECK_EQ(householder, expectedHouseholder);

    m.multiplyLeft(householder);

    CHECK_EQ(m, Matrix{{-1, -1}, {0, 1}});
}

TEST_CASE_FIXTURE(MixedStateSimplifierTest, "Small Householder matrix") {
    EditableMatrix m{{.5, .5, 0.}, {.5, .5, 0.}};

    HouseholderMatrix householder{Matrix{{.5}, {.5}}};

    UnitaryMatrix expectedHouseholder{{-1 / std::sqrt(2), -1 / std::sqrt(2)},
                                      {-1 / std::sqrt(2), 1 / std::sqrt(2)}};

    CHECK_EQ(householder, expectedHouseholder);

    m.multiplyLeft(householder);

    Matrix expectedProduct{{-1 / std::sqrt(2), -1 / std::sqrt(2), 0.},
                           {0., 0., 0.}};

    CHECK_EQ(m, expectedProduct);
}

TEST_CASE_FIXTURE(MixedStateSimplifierTest, "Triangularize") {
    EditableMatrix m{{1, 2}, {3, 4}};

    qrTriangularize(m);

    Matrix expectedTriangularized{{-3.1622776601683791, -4.4271887242357311},
                                  {0., 0.63245553203367644}};

    CHECK_EQ(m, expectedTriangularized);
}

TEST_CASE_FIXTURE(MixedStateSimplifierTest, "Triangularize large matrices") {
    checkTriangularization(
        Matrix{{1.1, 2, 3, 4, 5}, {1, 2, 3, 3, 5}, {1, 2, 3, 5, 5}});
    checkTriangularization(Matrix{{1.1, 2, 3, 4, 5, 1i},
                                  {1, 2, 3, 3, 5, 1i},
                                  {1, 2, 3, 5, 5, 1i},
                                  {12, 584, 84, 84, 51, 9.415i},
                                  {1i, 1i, 1i, 1i, 1i, 1i}});
    checkTriangularization(Matrix{
        {1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0},
        {1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0},
    });
    checkTriangularization(Matrix{{1i}, {1i}});
    checkTriangularization(Matrix{{1i}, {1i}, {1i}, {1i}, {1i}, {1i}, {1i},
                                  {1i}, {1i}, {1i}, {1i}, {1i}, {1i}, {1i},
                                  {1i}, {1i}, {1i}, {1i}, {1i}, {1i}, {1i}});
    checkTriangularization(Matrix{
        {1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i},
        {1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i},
        {1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i},
        {1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i},
        {1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i},
        {1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i},
        {1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i},
        {1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i},
        {1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i},
        {1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i},
        {1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i},
        {1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i},
        {1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i},
        {1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i},
        {1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i},
        {1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i},
        {1i, 1i, 53. + 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i},
        {1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i},
        {1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i},
        {1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i},
        {1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i, 1i},
    });
}

TEST_CASE_FIXTURE(MixedStateSimplifierTest, "SparseStateAdapter") {
    using T = absl::btree_map<Key<>, std::complex<double>>;

    SUBCASE("Triangularization does not result in fewer elements") {
        T data{{{.ensembleIndex = EnsembleIndex{0},
                 .measurementRegister = {"00"},
                 .basisVector = {"00"}},
                1.},
               {{.ensembleIndex = EnsembleIndex{0},
                 .measurementRegister = {"00"},
                 .basisVector = {"10"}},
                .5},
               {{.ensembleIndex = EnsembleIndex{1},
                 .measurementRegister = {"00"},
                 .basisVector = {"10"}},
                .5},
               {{.ensembleIndex = EnsembleIndex{1},
                 .measurementRegister = {"10"},
                 .basisVector = {"00"}},
                .2},
               {{.ensembleIndex = EnsembleIndex{1},
                 .measurementRegister = {"10"},
                 .basisVector = {"01"}},
                .1}};

        auto state = SparseStateAdapter<64>::makeIfNeeded(data);

        CHECK(!state.has_value());
    }

    SUBCASE("Triangularization sparsifies matrix") {
        T data{{{.ensembleIndex = EnsembleIndex{0},
                 .measurementRegister = {"10"},
                 .basisVector = {"00"}},
                1.},
               {{.ensembleIndex = EnsembleIndex{0},
                 .measurementRegister = {"10"},
                 .basisVector = {"10"}},
                .5},
               {{.ensembleIndex = EnsembleIndex{1},
                 .measurementRegister = {"10"},
                 .basisVector = {"00"}},
                .5},
               {{.ensembleIndex = EnsembleIndex{1},
                 .measurementRegister = {"10"},
                 .basisVector = {"10"}},
                .2}};

        auto state = SparseStateAdapter<64>::makeIfNeeded(data);

        REQUIRE(state.has_value());
        CHECK_EQ(*state, Matrix{{1., .5}, {.5, .2}});

        state->multiplyLeft(SquareMatrix::identity(2));

        CHECK_EQ(*state, Matrix{{1., .5}, {.5, .2}});

        state->multiplyLeft(Matrix{{1., 1.}, {1., 1.}});

        CHECK_EQ(*state, Matrix{{1.5, .7}, {1.5, .7}});

        std::complex<double> x = 0.;
        state->forEach([&](auto i, auto j, auto v) { x += v; });
        CHECK(utils::isNull(x - 4.4));
    }
}

TEST_CASE_FIXTURE(MixedStateSimplifierTest, "Simplify") {
    using T = absl::btree_map<Key<>, std::complex<double>>;
    T data{{{.ensembleIndex = EnsembleIndex{0},
             .measurementRegister = {"00"},
             .basisVector = {"00"}},
            .5},
           {{.ensembleIndex = EnsembleIndex{0},
             .measurementRegister = {"00"},
             .basisVector = {"01"}},
            .5},
           {{.ensembleIndex = EnsembleIndex{1},
             .measurementRegister = {"00"},
             .basisVector = {"00"}},
            .5},
           {{.ensembleIndex = EnsembleIndex{1},
             .measurementRegister = {"00"},
             .basisVector = {"01"}},
            .5},
           {{.ensembleIndex = EnsembleIndex{1},
             .measurementRegister = {"01"},
             .basisVector = {"01"}},
            .5}};

    auto oldData = data;

    MixedStateSimplifier::simplify(data);

    CHECK(areEquals(data, oldData));
}

TEST_CASE_FIXTURE(MixedStateSimplifierTest, "Sparsify 1" * doctest::skip()) {
    using T = absl::btree_map<Key<>, std::complex<double>>;
    T data{{{.ensembleIndex = EnsembleIndex{0},
             .measurementRegister = {"00"},
             .basisVector = {"00"}},
            1 / std::sqrt(2)},
           {{.ensembleIndex = EnsembleIndex{1},
             .measurementRegister = {"00"},
             .basisVector = {"00"}},
            1 / std::sqrt(2)}};

    MixedStateSimplifier::sparsify(data);

    CHECK(areEquals(data, T{{{.ensembleIndex = EnsembleIndex{0},
                              .measurementRegister = {"00"},
                              .basisVector = {"00"}},
                             -1.}}));
}

TEST_CASE_FIXTURE(MixedStateSimplifierTest, "Sparsify 2" * doctest::skip()) {
    using T = absl::btree_map<Key<>, std::complex<double>>;
    T data{{{.ensembleIndex = EnsembleIndex{0},
             .measurementRegister = {"00"},
             .basisVector = {"11"}},
            1 / std::sqrt(3)},
           {{.ensembleIndex = EnsembleIndex{1},
             .measurementRegister = {"00"},
             .basisVector = {"11"}},
            1 / std::sqrt(3)},
           {{.ensembleIndex = EnsembleIndex{3},
             .measurementRegister = {"00"},
             .basisVector = {"11"}},
            1 / std::sqrt(3)}};

    MixedStateSimplifier::sparsify(data);

    CHECK(areEquals(data, T{{{.ensembleIndex = EnsembleIndex{0},
                              .measurementRegister = {"00"},
                              .basisVector = {"11"}},
                             -1.}}));
}

TEST_CASE_FIXTURE(MixedStateSimplifierTest, "Sparsify 3" * doctest::skip()) {
    using T = absl::btree_map<Key<>, std::complex<double>>;
    T data{{{.ensembleIndex = EnsembleIndex{0},
             .measurementRegister = {"00"},
             .basisVector = {"00"}},
            .5},
           {{.ensembleIndex = EnsembleIndex{0},
             .measurementRegister = {"00"},
             .basisVector = {"01"}},
            .5},
           {{.ensembleIndex = EnsembleIndex{1},
             .measurementRegister = {"00"},
             .basisVector = {"00"}},
            .5},
           {{.ensembleIndex = EnsembleIndex{1},
             .measurementRegister = {"00"},
             .basisVector = {"01"}},
            .5}};

    MixedStateSimplifier::sparsify(data);

    CHECK(areEquals(data, T{{{.ensembleIndex = EnsembleIndex{0},
                              .measurementRegister = {"00"},
                              .basisVector = {"00"}},
                             -1 / std::sqrt(2)},
                            {{.ensembleIndex = EnsembleIndex{0},
                              .measurementRegister = {"00"},
                              .basisVector = {"01"}},
                             -1 / std::sqrt(2)}}));
}

TEST_CASE_FIXTURE(MixedStateSimplifierTest,
                  "Sparsify doesn't modify" * doctest::skip()) {
    using T = absl::btree_map<Key<>, std::complex<double>>;
    T data{{{.ensembleIndex = EnsembleIndex{0},
             .measurementRegister = {"00"},
             .basisVector = {"00"}},
            .5},
           {{.ensembleIndex = EnsembleIndex{0},
             .measurementRegister = {"00"},
             .basisVector = {"01"}},
            .5},
           {{.ensembleIndex = EnsembleIndex{1},
             .measurementRegister = {"00"},
             .basisVector = {"00"}},
            .5}};

    auto oldData = data;

    MixedStateSimplifier::sparsify(data);

    CHECK(areEquals(data, oldData));
}

TEST_CASE_FIXTURE(MixedStateSimplifierTest, "Sparsify Givens 1") {
    using T = absl::btree_map<Key<>, std::complex<double>>;
    T data{{{.ensembleIndex = EnsembleIndex{0},
             .measurementRegister = {"00"},
             .basisVector = {"00"}},
            1 / std::sqrt(2)},
           {{.ensembleIndex = EnsembleIndex{1},
             .measurementRegister = {"00"},
             .basisVector = {"00"}},
            1 / std::sqrt(2)}};

    MixedStateSimplifier::sparsifyGivens(data);

    CHECK(areEquals(data, T{{{.ensembleIndex = EnsembleIndex{0},
                              .measurementRegister = {"00"},
                              .basisVector = {"00"}},
                             1.}}));
}

TEST_CASE_FIXTURE(MixedStateSimplifierTest, "Sparsify Givens 2") {
    using T = absl::btree_map<Key<>, std::complex<double>>;
    T data{{{.ensembleIndex = EnsembleIndex{0},
             .measurementRegister = {"00"},
             .basisVector = {"11"}},
            1 / std::sqrt(3)},
           {{.ensembleIndex = EnsembleIndex{1},
             .measurementRegister = {"00"},
             .basisVector = {"11"}},
            1 / std::sqrt(3)},
           {{.ensembleIndex = EnsembleIndex{3},
             .measurementRegister = {"00"},
             .basisVector = {"11"}},
            1 / std::sqrt(3)}};

    MixedStateSimplifier::sparsifyGivens(data);

    CHECK(areEquals(data, T{{{.ensembleIndex = EnsembleIndex{0},
                              .measurementRegister = {"00"},
                              .basisVector = {"11"}},
                             1.}}));
}

TEST_CASE_FIXTURE(MixedStateSimplifierTest, "Sparsify Givens 3") {
    using T = absl::btree_map<Key<>, std::complex<double>>;
    T data{{{.ensembleIndex = EnsembleIndex{0},
             .measurementRegister = {"00"},
             .basisVector = {"00"}},
            .5},
           {{.ensembleIndex = EnsembleIndex{0},
             .measurementRegister = {"00"},
             .basisVector = {"01"}},
            .5},
           {{.ensembleIndex = EnsembleIndex{1},
             .measurementRegister = {"00"},
             .basisVector = {"00"}},
            .5},
           {{.ensembleIndex = EnsembleIndex{1},
             .measurementRegister = {"00"},
             .basisVector = {"01"}},
            .5}};

    MixedStateSimplifier::sparsifyGivens(data);

    CHECK(areEquals(data, T{{{.ensembleIndex = EnsembleIndex{0},
                              .measurementRegister = {"00"},
                              .basisVector = {"00"}},
                             1 / std::sqrt(2)},
                            {{.ensembleIndex = EnsembleIndex{0},
                              .measurementRegister = {"00"},
                              .basisVector = {"01"}},
                             1 / std::sqrt(2)}}));
}

} // namespace qx