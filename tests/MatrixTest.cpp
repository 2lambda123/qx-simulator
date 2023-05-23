#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"
#include "qx/Core.hpp"

namespace qx {
namespace core {

using namespace std::complex_literals;

class MatrixTest {};

TEST_CASE_FIXTURE(MatrixTest, "Identity") {
    SquareMatrix expectedIdentity{{1, 0, 0, 0, 0},
                                  {0, 1, 0, 0, 0},
                                  {0, 0, 1, 0, 0},
                                  {0, 0, 0, 1, 0},
                                  {0, 0, 0, 0, 1}};

    CHECK_EQ(SquareMatrix::identity(5), expectedIdentity);
}

TEST_CASE_FIXTURE(MatrixTest, "Simple generic product") {
    EditableMatrix a{{1, 2}};
    
    Matrix b{{1, 2},
             {1, 2}};
    
    a *= b;

    Matrix aTimesB{{3, 6}};

    CHECK_EQ(a, aTimesB);
}

TEST_CASE_FIXTURE(MatrixTest, "Generic product") {
    EditableMatrix a{{1, 1, 1},
                     {1, 1, 1},
                     {1, 1, 1},
                     {1, 1, 1}};
    
    Matrix b{{1, 1, 1},
             {1, 1, 1},
             {1, 1, 1}};
    
    a *= b;

    Matrix aTimesB{{3, 3, 3},
                   {3, 3, 3},
                   {3, 3, 3},
                   {3, 3, 3}};

    CHECK_EQ(a, aTimesB);
}

TEST_CASE_FIXTURE(MatrixTest, "Multiply left") {
    EditableMatrix a{{1}, {2}};
    
    Matrix b{{1, 2},
             {1, 2}};
    
    a.multiplyLeft(b);

    Matrix aTimesB{{5}, {5}};

    CHECK_EQ(a, aTimesB);
}

TEST_CASE_FIXTURE(MatrixTest, "Product") {
    Matrix a{{1, 2},
             {3, 4},
             {5.5 + 1i, 6}};
    
    Matrix b{{1},
             {3}};
    
    Matrix aTimesB{{7},
                   {15},
                   {23.5 + 1i}};

    CHECK_EQ(a * b, aTimesB);
}

TEST_CASE_FIXTURE(MatrixTest, "Sum") {
    EditableMatrix a{{1, 2}};
    
    Matrix b{{2, 3.5}};
    
    Matrix aPlusB{{3, 5.5}};

    a += b;
    CHECK_EQ(a, aPlusB);
}

TEST_CASE_FIXTURE(MatrixTest, "Submatrix") {
    Matrix victim{{1, 2, 3},
                  {1, 2, 3},
                  {1, 4, 5},
                  {1, 2, 3}};

    SubMatrixAdapter::Coordinates coordinates{
        .startRow = 1,
        .startCol = 1,
        .endRow = 3,
        .endCol = 3
    };
    SubMatrixAdapter subMatrix(victim, coordinates);

    Matrix expected{{2, 3},
                    {4, 5}};

    CHECK_EQ(subMatrix, expected);
}

TEST_CASE_FIXTURE(MatrixTest, "Editable submatrix") {
    EditableMatrix victim{{1, 2},
                  {1, 2}};

    EditableSubMatrixAdapter::Coordinates coordinates{
        .startRow = 1,
        .startCol = 1,
        .endRow = 2,
        .endCol = 2
    };
    EditableSubMatrixAdapter subMatrix(victim, coordinates);

    subMatrix += Matrix{{2}};

    Matrix expected{{1, 2},
                    {1, 4}};

    CHECK_EQ(victim, expected);
}

TEST_CASE_FIXTURE(MatrixTest, "Reject non-unitary matrices") {
    std::initializer_list<std::initializer_list<std::complex<double>>> test1 = {{1, 0}, {0, 1.1}};
    CHECK_THROWS_WITH([test1] { UnitaryMatrix test(test1); }(),
        "Matrix is not unitary");

    std::initializer_list<std::initializer_list<std::complex<double>>> test2 = {{0, 0}, {0, 1}};
    CHECK_THROWS_WITH(
        [test2] { UnitaryMatrix test(test2); }(),
        "Matrix is not unitary");

    std::initializer_list<std::initializer_list<std::complex<double>>> test3 = {{1i, 0, 0}, {0, 1i, 1i}, {0, 1, 0}};
    CHECK_THROWS_WITH(
        [test3] { UnitaryMatrix test(test3); }(),
        "Matrix is not unitary");
}

TEST_CASE_FIXTURE(MatrixTest, "Dagger") {
    Matrix m{{1 / std::sqrt(2), 1 / std::sqrt(2)},
            {1i / std::sqrt(2), -1i / std::sqrt(2)}};

    Matrix mDag{{1 / std::sqrt(2), -1i / std::sqrt(2)},
            {1 / std::sqrt(2), 1i / std::sqrt(2)}};

    CHECK_EQ(m.dagger(), mDag);
}

TEST_CASE_FIXTURE(MatrixTest, "Spectral radius") {
    Matrix m0{{1 / std::sqrt(2), 0.},
              {0., -1i / std::sqrt(2)}};

    CHECK_EQ(computeSpectralRadius(m0), ::doctest::Approx(1 / std::sqrt(2)));

    Matrix m1{{1, 0.},
              {0., -2}};

    CHECK_EQ(computeSpectralRadius(m1), ::doctest::Approx(2));

    Matrix m2{{1, -10},
              {0., -2i}};

    CHECK_EQ(computeSpectralRadius(m2), ::doctest::Approx(2));

    Matrix m3{{1, -10, 1, 2},
              {0., -2, 4, 2i},
              {0, 0, 4i, 23},
              {0, 0, 0, 1i}};

    CHECK_EQ(computeSpectralRadius(m3), ::doctest::Approx(4));

    Matrix m4{{0, 1, -2},
              {0, 1, 0},
              {1, -1, 3}};

    CHECK_EQ(computeSpectralRadius(m4), ::doctest::Approx(2));
}

} // namespace core
} // namespace qx