#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"
#include "qx/Core.hpp"
#include "qx/Gates.hpp"

namespace qx {
namespace core {

using namespace std::complex_literals;

class StatisticalEnsembleTest {
protected:
    static constexpr DenseMatrix<2> const M0{{{{1, 0},
                                               {0, 0}}}};

    static constexpr DenseMatrix<2> const M1{{{{0, 0},
                                            {0, 1}}}};

    static void checkEqComplex(std::complex<double> left, std::complex<double> right) {
        CHECK_EQ(left.real(), ::doctest::Approx(right.real()));
        CHECK_EQ(left.imag(), ::doctest::Approx(right.imag()));
    }

    template <std::size_t N>
    static void checkMatrix(StatisticalEnsemble const& victim, std::vector<std::array<std::complex<double>, N>> expected) {
        auto actual = victim.testToMatrix<N>();

        REQUIRE_EQ(actual.size(), expected.size());

        for (std::size_t i = 0; i < actual.size(); ++i) {
            for (std::size_t j = 0; j < N; ++j) {
                CAPTURE(i);
                CAPTURE(j);
                checkEqComplex(actual[i][j], expected[i][j]);
            }
        }
    }
};

TEST_CASE_FIXTURE(StatisticalEnsembleTest, "Identity") {
    StatisticalEnsemble victim(1);

    DenseMatrix<2> identity({{{1, 0},
                              {0, 1}}});

    victim.applyKrausOperators<1>({identity}, {0});
    CHECK_EQ(victim, StatisticalEnsemble(1));
}

TEST_CASE_FIXTURE(StatisticalEnsembleTest, "Kraus operators not satisfying the completeness relation") {
    StatisticalEnsemble victim(1);

    DenseMatrix<2> k0({{{1, 0},
                        {0, 0}}});

    DenseMatrix<2> k1({{{0, 1},
                        {0, 1}}});

    CHECK(!StatisticalEnsemble::areValidKrausOperators<1>({k0, k1}));
}

TEST_CASE_FIXTURE(StatisticalEnsembleTest, "Measure in computational basis") {
    StatisticalEnsemble victim(1);

    auto probabilities = victim.applyKrausOperatorsGetOutcomeProbabilities<1>({M0, M1}, {0});

    DenseMatrix<2> expected({{{1, 0},
                              {0, 0}}});

    CHECK_EQ(victim.testToDensityMatrix<2>(), expected);

    CHECK_EQ(probabilities.size(), 2);
    CHECK_EQ(probabilities[0], ::doctest::Approx(1.));
    CHECK_EQ(probabilities[1], ::doctest::Approx(0.));
}

TEST_CASE_FIXTURE(StatisticalEnsembleTest, "Conditional bit-shift") {
    StatisticalEnsemble victim(1);

    double p = 0.2;
    DenseMatrix<2> k0 = std::sqrt(1 - p) * DenseMatrix<2>::identity();

    DenseMatrix<2> k1 = std::sqrt(p) * gates::X;

    victim.applyKrausOperators<1>({k0, k1}, {0});
    
    DenseMatrix<2> expected({{{1-p, 0},
                              {0,   p}}});

    CHECK_EQ(victim.testToDensityMatrix<2>(), expected);
}

TEST_CASE_FIXTURE(StatisticalEnsembleTest, "|+> state") {
    StatisticalEnsemble victim(1);

    victim.applyKrausOperators<1>({gates::H}, {0});
    
    DenseMatrix<2> expected({{{1/2., 1/2.},
                              {1/2., 1/2.}}});

    CHECK_EQ(victim.testToDensityMatrix<2>(), expected);
}

TEST_CASE_FIXTURE(StatisticalEnsembleTest, "|+> state with 2 qubits") {
    StatisticalEnsemble victim(2);

    victim.applyKrausOperators<1>({gates::H}, {1});
    
    DenseMatrix<4> expected({{{1/2., 0, 1/2., 0},
                              {0, 0, 0, 0},
                              {1/2., 0, 1/2., 0},
                              {0, 0, 0, 0}}});

    CHECK_EQ(victim.testToDensityMatrix<4>(), expected);
}

TEST_CASE_FIXTURE(StatisticalEnsembleTest, "Measure bell pair") {
    StatisticalEnsemble victim(3);

    victim.applyKrausOperators<1>({gates::H}, {1});
    victim.applyKrausOperators<2>({gates::CNOT}, {1, 2});

    auto probabilities = victim.applyKrausOperatorsGetOutcomeProbabilities<1>({M0, M1}, {2});

    CHECK_EQ(probabilities.size(), 2);
    CHECK_EQ(probabilities[0], ::doctest::Approx(0.5));
    CHECK_EQ(probabilities[1], ::doctest::Approx(0.5));
}

TEST_CASE_FIXTURE(StatisticalEnsembleTest, "Measure mixed state") {
    StatisticalEnsemble victim(3);

    double p = 0.5;
    DenseMatrix<2> k0 = std::sqrt(1 - p) * DenseMatrix<2>::identity();

    DenseMatrix<2> k1 = std::sqrt(p) * gates::H;

    victim.applyKrausOperators<1>({k0, k1}, {1});

    auto probabilities = victim.applyKrausOperatorsGetOutcomeProbabilities<1>({M0, M1}, {1});

    CHECK_EQ(probabilities.size(), 2);
    CHECK_EQ(probabilities[0], ::doctest::Approx(0.75));
    CHECK_EQ(probabilities[1], ::doctest::Approx(0.25));
}

TEST_CASE_FIXTURE(StatisticalEnsembleTest, "Shrink pure state") {
    StatisticalEnsemble victim(2);

    victim.applyKrausOperators<1>({gates::H}, {0});
    victim.applyKrausOperators<2>({gates::CNOT}, {0, 1});
    
    checkMatrix<4>(victim, {{1/std::sqrt(2), 0., 0., 1/std::sqrt(2)}});

    victim.shrink();

    checkMatrix<4>(victim, {{1/std::sqrt(2), 0., 0., 1/std::sqrt(2)}});
}

TEST_CASE_FIXTURE(StatisticalEnsembleTest, "Shrink mixed state") {
    StatisticalEnsemble victim(1);

    double p = 0.5;
    DenseMatrix<2> k0 = std::sqrt(1 - p) * DenseMatrix<2>::identity();
    DenseMatrix<2> k1 = std::sqrt(p) * gates::H;
    victim.applyKrausOperators<1>({k0, k1}, {0});

    checkMatrix<2>(victim, {{1/std::sqrt(2), 0}, {1/2., 1/2.}});

    victim.applyKrausOperators<1>({k0, k1}, {0});
    
    DenseMatrix<2> expectedDensityMatrix({{{3/4., 1/4.},
                                           {1/4., 1/4.}}});

    CHECK_EQ(victim.testToDensityMatrix<2>(), expectedDensityMatrix);

    checkMatrix<2>(victim, {{1/2., 0}, {0.5/std::sqrt(2), 0.5/std::sqrt(2)}, {0.5/std::sqrt(2), 0.5/std::sqrt(2)}, {1/2., 0}});

    victim.shrink();

    checkMatrix<2>(victim, {{0.86602540378443871, 0.28867513459481292}, {0, 0.40824829046386313}});
    CHECK_EQ(victim.testToDensityMatrix<2>(), expectedDensityMatrix);
}

TEST_CASE_FIXTURE(StatisticalEnsembleTest, "Shrink very mixed state") {
    StatisticalEnsemble victim(2);

    double p = 0.5;
    DenseMatrix<2> k0 = std::sqrt(1 - p) * DenseMatrix<2>::identity();
    DenseMatrix<2> k1 = std::sqrt(p) * gates::H;
    victim.applyKrausOperators<1>({k0, k1}, {0});
    victim.applyKrausOperators<1>({k0, k1}, {1});
    victim.applyKrausOperators<1>({k0, k1}, {1});
    victim.applyKrausOperators<1>({k0, k1}, {0});
    victim.applyKrausOperators<1>({k0, k1}, {0});
    victim.applyKrausOperators<1>({k0, k1}, {1});
    victim.applyKrausOperators<1>({k0, k1}, {0});
    victim.applyKrausOperators<1>({k0, k1}, {1});
    victim.applyKrausOperators<1>({k0, k1}, {0});
    
    CHECK_EQ(victim.getEnsembleSize(), 512);

    auto densityMatrixBeforeShrink = victim.testToDensityMatrix<4>();

    victim.shrink();

    CHECK_EQ(victim.getEnsembleSize(), 4);
    checkMatrix<4>(victim, {{3/4., 1/4., 1/4., 1/12.},
                            {0., .5/std::sqrt(2), 0, 1./(6*std::sqrt(2))}, {0., 0., .5/std::sqrt(2), 1./(6*std::sqrt(2))}, {0., 0., 0., 1/6.}});

    CHECK_EQ(victim.testToDensityMatrix<4>(), densityMatrixBeforeShrink);
}

TEST_CASE_FIXTURE(StatisticalEnsembleTest, "Shrink simplifies rho = .5 |+><+| + .5 |-><-| to rho = .5 |0><0| + .5 |1><1|") {
    StatisticalEnsemble victim(10);

    DenseMatrix<2> k0 = std::sqrt(.5) * DenseMatrix<2>::identity();
    DenseMatrix<2> k1 = std::sqrt(.5) * gates::X;
    victim.applyKrausOperators<1>({k0, k1}, {5});
    victim.applyKrausOperators<1>({gates::H}, {5});
    
    auto actualMatrix = victim.testToMatrix<1024>();
    CHECK_EQ(actualMatrix.size(), 2);
    checkEqComplex(actualMatrix[0][0], .5);
    checkEqComplex(actualMatrix[0][32], .5);
    checkEqComplex(actualMatrix[1][0], .5);
    checkEqComplex(actualMatrix[1][32], -.5);

    victim.shrink();

    actualMatrix = victim.testToMatrix<1024>();
    CHECK_EQ(actualMatrix.size(), 2);
    checkEqComplex(actualMatrix[0][0], 1/std::sqrt(2));
    checkEqComplex(actualMatrix[1][32], 1/std::sqrt(2));
}


} // namespace core
} // namespace qx