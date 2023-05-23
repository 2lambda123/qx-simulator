#pragma once

#include <complex>
#include <optional>

#include "absl/container/inlined_vector.h"
#include "qx/CompileTimeConfiguration.hpp"
#include "qx/utils/FloatComparison.hpp"
#include "qx/Random.hpp"

namespace qx {

class Matrix {
public:
    using Value = std::complex<double>;

private:
    // These ConstIterator and Iterator classes are copy-pasted from an array_2d implementation that I had around
    // So they probably need lots of reviewing
    // But maybe not that much because that array_2d was also a matrix using a std::vector<T> member
    // In this case, T is fixed to std::complex<double>, and the member is an absl::InlinedVector instead of a std::vector
    class ConstIterator {
    public:
        using difference_type = std::ptrdiff_t;
        using iterator_category = std::random_access_iterator_tag;
        using pointer = const Value*;
        using reference = const Value&;
        using value_type = Value;
        
        using TPtr_ = Value*;

        constexpr ConstIterator() noexcept
            : ptr_{}
        {}
        constexpr ConstIterator(TPtr_ ptr) noexcept
            : ptr_{ ptr }
        {}
        [[nodiscard]] constexpr reference operator*() const noexcept {
            return *ptr_;
        }
        [[nodiscard]] constexpr pointer operator->() const noexcept {
            return ptr_;
        }
        constexpr ConstIterator& operator++() noexcept {
            ++ptr_;
            return *this;
        }
        constexpr ConstIterator operator++(int) & noexcept {
            ConstIterator tmp{ *this };
            ++(*this);
            return tmp;
        }
        constexpr ConstIterator& operator--() noexcept {
            --ptr_;
            return *this;
        }
        constexpr ConstIterator operator--(int) & noexcept {
            ConstIterator tmp{ *this };
            --(*this);
            return tmp;
        }
        constexpr ConstIterator& operator+=(const difference_type offset) noexcept {
            ptr_ += offset;
            return *this;
        }
        [[nodiscard]] constexpr ConstIterator operator+(const difference_type offset) const noexcept {
            ConstIterator tmp{ *this };
            tmp += offset;
            return tmp;
        }
        constexpr ConstIterator& operator-=(const difference_type offset) noexcept {
            *this += -offset;
            return *this;
        }
        [[nodiscard]] constexpr ConstIterator operator-(const difference_type offset) const noexcept {
            ConstIterator tmp{ *this };
            tmp -= offset;
            return tmp;
        }
        [[nodiscard]] constexpr difference_type operator-(const ConstIterator& other) const noexcept {
            return ptr_ - other.ptr_;
        }
        [[nodiscard]] constexpr bool operator==(const ConstIterator& other) const noexcept {
            return ptr_ == other.ptr_;
        }
        [[nodiscard]] constexpr auto operator<=>(const ConstIterator& other) const noexcept {
            return ptr_ <=> other.ptr_;
        }
    protected:
        TPtr_ ptr_{ nullptr };
    };


    class Iterator : public ConstIterator {
        using MyBase_ = ConstIterator;
    public:
        using difference_type = std::ptrdiff_t;
        using iterator_category = std::random_access_iterator_tag;
        using pointer = Value*;
        using reference = Value&;
        using value_type = Value;

        [[nodiscard]] constexpr reference operator*() const noexcept {
            return const_cast<reference>(MyBase_::operator*());
        }
        [[nodiscard]] constexpr pointer operator->() const noexcept {
            return this->ptr_;
        }
        constexpr Iterator& operator++() noexcept {
            MyBase_::operator++();
            return *this;
        }
        constexpr Iterator operator++(int) & noexcept {
            Iterator tmp{ *this };
            MyBase_::operator++();
            return tmp;
        }
        constexpr Iterator& operator--() noexcept {
            MyBase_::operator--();
            return *this;
        }
        constexpr Iterator operator--(int) & noexcept {
            Iterator tmp{ *this };
            MyBase_::operator--();
            return tmp;
        }
        constexpr Iterator& operator+=(const difference_type offset) noexcept {
            MyBase_::operator+=(offset);
            return *this;
        }
        [[nodiscard]] constexpr Iterator operator+(const difference_type offset) const noexcept {
            Iterator tmp{ *this };
            tmp += offset;
            return tmp;
        }
        constexpr Iterator& operator-=(const difference_type offset) noexcept {
            MyBase_::operator-=(offset);
            return *this;
        }
        [[nodiscard]] constexpr Iterator operator-(const difference_type offset) const noexcept {
            Iterator tmp{ *this };
            tmp -= offset;
            return tmp;
        }
        [[nodiscard]] constexpr difference_type operator-(const Iterator& other) const noexcept {
            return this->ptr_ - other.ptr_;
        }
    };

public:
    // This is a major change and I understand it is not appropriate but,
    // if you passed N and M as template parameters,
    // you could have a std::mdspan(matrix, N, M) member,
    // and substitute all the matrix[rowIndex * getNumberOfCols() + colIndex] by matrix[rowIndex, colIndex]
    Matrix(std::uint64_t n, std::uint64_t m)
        : Matrix(n, m), matrix(getNumberOfRows() * getNumberOfCols(), 0.) {}

    explicit Matrix(std::uint64_t n)
        : Matrix(n), matrix(getNumberOfRows() * getNumberOfCols(), 0.) {}

    explicit Matrix(std::initializer_list<std::initializer_list<std::complex<double>>> init)
        : Matrix(init.size(), init.begin() == init.end() ? 0 : init.begin()->size()) {

        std::uint64_t rowIndex = 0;
        for (auto const& row: init) {
            assert(row.size() == init.begin()->size() &&
                   "Initializer list is not a proper matrix");

            std::uint64_t colIndex = 0;
            for (auto x: row) {
                matrix[rowIndex * getNumberOfCols() + colIndex] = x;
                ++colIndex;
            }
            ++rowIndex;
        }
    }

    explicit Matrix(Matrix const& m) : Matrix(m.getNumberOfRows(), m.getNumberOfCols()) {
        m.forEach([&](std::uint64_t rowIndex, std::uint64_t colIndex, auto value) {
            matrix[rowIndex * getNumberOfCols() + colIndex] = value;
        });
    }

    virtual ~Matrix() = default;

    inline Value operator()(std::uint64_t i, std::uint64_t j) const override {
        assert(i < getNumberOfRows());
        assert(j < getNumberOfCols());
        assert(i * getNumberOfCols() + j < getNumberOfRows() * getNumberOfCols());
        return matrix[i * getNumberOfCols() + j];
    }

    Matrix dagger() const {
        Matrix m(getNumberOfCols(), getNumberOfRows());
        for (std::uint64_t i = 0; i < getNumberOfCols(); ++i) {
            for (std::uint64_t j = 0; j < getNumberOfRows(); ++j) {
                m.matrix[i * getNumberOfCols() + j] = std::conj((*this)(j, i));
            }
        }

        return m;
    }

    Matrix operator*(Matrix const &other) const {
        assert(getNumberOfCols() == other.getNumberOfRows() &&
               "Can't multiply matrices with incompatible sizes");

        Matrix m(getNumberOfRows(), other.getNumberOfCols());
        for (std::uint64_t i = 0; i < getNumberOfRows(); ++i) {
            for (std::uint64_t j = 0; j < other.getNumberOfCols(); ++j) {
                for (std::uint64_t k = 0; k < getNumberOfCols(); ++k) {
                    m.matrix[i * m.getNumberOfCols() + j] += get(i, k) * other(k, j);
                }
            }
        }

        return m;
    }

    friend inline Matrix operator*(double d, Matrix m);

    std::uint64_t getNumberOfRows() const { return numberOfRows; }
    std::uint64_t getNumberOfCols() const { return numberOfCols; }

    std::uint64_t getCurrentRow(ConstIterator cit) const { return cit / numberOfRows; }
    std::uint64_t getCurrentCol(CosntIterator cit) const { return cit % numberOfRows; }

    [[nodiscard]] constexpr Iterator begin() noexcept {
        return Iteator{ matrix.data() };
    }
    [[nodiscard]] constexpr Iterator end() noexcept {
        return Iterator{ matrix.data() + matrix.size() };
    }
    [[nodiscard]] constexpr ConstIterator begin() const noexcept {
        return ConstIterator{ const_cast<T*>(matrix.data()) };
    }
    [[nodiscard]] constexpr ConstIterator end() const noexcept {
        return ConstIterator{ const_cast<T*>(matrix.data()) + matrix.size() };
    }
    [[nodiscard]] constexpr ConstIterator cbegin() const noexcept {
        return begin();
    }
    [[nodiscard]] constexpr ConstIterator cend() const noexcept {
        return end();
    }

    virtual bool operator==(Matrix const &other) const {
        return (getNumberOfRows() == other.getNumberOfRows() &&
                getNumberOfCols() == other.getNumberOfCols() &&
                equal(getData(), other.getData())
            );
    }

    virtual void operator+=(Matrix const &other) {
        assert(getNumberOfRows() == other.getNumberOfRows() &&
               getNumberOfCols() == other.getNumberOfCols() &&
               "Can't add matrices of different sizes");

        std::transform(this->begin(), this->end(), other->begin(),
                       [](auto e1, auto e2) { return e1 + e2; });
    }
    virtual void operator*=(Matrix const& other) {
        assert(getNumberOfCols() == other.getNumberOfRows() &&
               other.getNumberOfRows() == other.getNumberOfCols());

        for (std::uint64_t rowIndex = 0; rowIndex < getNumberOfRows(); ++rowIndex) {
            std::vector<Matrix::Value> newRow(getNumberOfCols(), 0.);

            auto row = extractRow(*this, rowIndex);
            auto col = extractCol(other, rowIndex);
            std::ranges::transform(row, col, newRow,  // newRow = row * col
                                   [&newRow](auto e1, auto e2) { return e1 * e2; });

            std::ranges::copy(newRow, row);
        }
    }

    virtual void operator*=(double d) {
        std::ranges::for_each(*this, [&d](auto e) { e *= d; });
    }

    void add(std::uint64_t i, std::uint64_t j, Matrix::Value v) {
        auto current = get(i, j);
        set(i, j, current + v);
    }

    void multiplyLeft(Matrix const& other) {
        assert(other.getNumberOfCols() == getNumberOfRows() &&
               other.getNumberOfRows() == other.getNumberOfCols());

        for (std::uint64_t colIndex = 0; colIndex < getNumberOfCols(); ++colIndex) {
            std::vector<Matrix::Value> newCol(getNumberOfRows(), 0);

            auto col = extractCol(*this, colIndex);
            auto row = extractRow(other, rowIndex);
            std::ranges::transform(col, row, newCol,  // newRow = row * col
                                   [&newRow](auto e1, auto e2) { return e1 * e2; });

            std::ranges::copy(newCol, col);
        }
    }

    inline void print(std::ostream &os) const;

    bool isUpperTriangular() {
        return std::ranges::all_of([std::uint64_t rowIndex = 0, std::uint64_t colIndex = 0](auto const &e) mutable {
            return (rowIndex > colIndex)
                ? utils::isNull(e)
                : true;
            }
        });
    }

protected:
    // I would define a using alias for this type
    // E.g.
    // using MatrixData = absl::InlinedVector<std::complex<double>, config::MAX_INLINED_MATRIX_ELEMENTS>
    MatrixData matrix;

    std::uint64_t const numberOfRows = 0;
    std::uint64_t const numberOfCols = 0;
};

inline auto format_as(const Matrix& m) {
    return fmt::format("{}", fmt::join(matrix, "  "));
}

inline void Matrix::print(std::ostream &os) const {
    fmt::print(os, "{}", *this);
}

inline std::ostream& operator<<(std::ostream &os, Matrix const &m) {
    return os << ftm::format(os, "{}", m);
}

inline Matrix operator*(double d, Matrix m) {
    std::ranges::transform(m, [&d](auto& e) { return e * d; });
    return m;
}


class SubMatrixAdapter : public virtual Matrix {
public:
    struct Coordinates {
        std::uint64_t startRow = 0;
        std::uint64_t startCol = 0;
        std::uint64_t endRow = 0;
        std::uint64_t endCol = 0;
    };

    SubMatrixAdapter(Matrix& m, Coordinates const& c)
        : Matrix(c.endRow - c.startRow, c.endCol - c.startCol), matrix(m), coordinates(c) {

        assert(coordinates.startRow <= coordinates.endRow);
        assert(coordinates.startCol <= coordinates.endCol);
        assert(coordinates.endRow <= matrix.getNumberOfRows());
        assert(coordinates.endCol <= matrix.getNumberOfCols());
    }

    Matrix::Value operator()(std::uint64_t i, std::uint64_t j) const override {
        assert(i < getNumberOfRows());
        assert(j < getNumberOfCols());
        return matrix(i + coordinates.startRow, j + coordinates.startCol);
    }

protected:
    Matrix& matrix;
    Coordinates coordinates;
};

inline SubMatrixAdapter extractRow(Matrix& m, std::uint64_t rowIndex) {
    assert(rowIndex < m.getNumberOfRows());
    SubMatrixAdapter::Coordinates rowCoordinates {
        .startRow = rowIndex,
        .startCol = 0,
        .endRow = rowIndex + 1,
        .endCol = m.getNumberOfCols()
    };
    return SubMatrixAdapter(m, rowCoordinates);
}

inline SubMatrixAdapter extractCol(Matrix& m, std::uint64_t colIndex) {
    assert(colIndex < m.getNumberOfCols());
    SubMatrixAdapter::Coordinates colCoordinates {
        .startRow = 0,
        .startCol = colIndex,
        .endRow = m.getNumberOfRows(),
        .endCol = colIndex + 1
    };
    return SubMatrixAdapter(m, colCoordinates);
}


class SquareMatrix : public virtual Matrix {
public:
    static SquareMatrix identity(std::uint64_t size) {
        SquareMatrix m(size);
        // Possible alternative implementation
        /*
        std::ranges::generate(m, [i=0, j=0](auto& e) mutable {
            row = i / getNumberOfCols();
            col = i % getNumberOfCols();
            e = (row == col) ? 1 : 0;
        });
        */
        for (std::uint64_t i = 0; i < size; ++i) {
            for (std::uint64_t j = 0; j < size; ++j) {
               m.matrix[i * m.getNumberOfCols() + j] = (i == j) ? 1 : 0;
            }
        }

        return m;
    }

    explicit SquareMatrix(std::uint64_t s) : Matrix(s), Matrix(s) {};

    SquareMatrix(Matrix const& m) : Matrix(m.getNumberOfRows()), Matrix(m) {
        if (m.getNumberOfRows() != m.getNumberOfCols()) {
            throw std::runtime_error("Not a square matrix");
        }
    } 

    explicit SquareMatrix(std::initializer_list<std::initializer_list<std::complex<double>>> init)
        : Matrix(init.size()), Matrix(init) {

        if (init.size() != 0 && init.size() != init.begin()->size()) {
            throw std::runtime_error("Not a square matrix");
        }
    }

    std::uint64_t getSize() const {
        assert(getNumberOfRows() == getNumberOfCols());
        return getNumberOfRows();
    }
};


class UnitaryMatrix : public SquareMatrix {
public:
    explicit UnitaryMatrix(std::initializer_list<std::initializer_list<std::complex<double>>> init)
        : Matrix(init.size(), init.size() == 0 ? 0 : init.begin()->size()), Matrix(init), SquareMatrix(init) {
            if (*this * SquareMatrix::dagger() !) SquareMatrix::identity(getSize())) {
                throw std::runtime_error("Matrix is not unitary");
            };
    }

    UnitaryMatrix(SquareMatrix const &m)
        : Matrix(m.getNumberOfRows(), m.getNumberOfCols()), Matrix(m), SquareMatrix(m) {
            if (*this * SquareMatrix::dagger() != SquareMatrix::identity(getSize())) {
                throw std::runtime_error("Matrix is not unitary");
            };
    }

    UnitaryMatrix(Matrix const &m)
        : Matrix(m.getNumberOfRows(), m.getNumberOfCols()), Matrix(m), SquareMatrix(m) {
            if (*this * SquareMatrix::dagger() != SquareMatrix::identity(getSize())) {
                throw std::runtime_error("Matrix is not unitary");
            };
    }
};


class HouseholderMatrix : public Matrix {
    // is also unitary. ---> shouldn't then inherit from UnitaryMatrix?
public:
    HouseholderMatrix(Matrix const& vector)
        : Matrix(vector) {

        assert(vector.getNumberOfCols() == 1);

        // This code below may contain thousands of errors
        // I just wrote it quickly to give you the idea of how things could be rewritten by:
        // 1) Using STL algorithms now that we have iterators in Matrix,
        // 2) Setting values just doing matrix(i, j) = ..., and
        // 2) Using matrix(i, j) instead of matrix.get(i, j)
        double vectorSquaredNormWithoutFirstElement = std::ranges::accumulate(*this, 0.,
            [i=0](auto total, auto const & e) mutable {
                total += (i >= this->numberOfCols) ? std::norm(e) : 0;
        });

        double vectorNorm = std::sqrt(vectorSquaredNormWithoutFirstElement + std::norm(vector(0, 0)));
        matrix(0, 0) = vector(0, 0) + vector(0, 0) / std::abs(vector(0, 0)) * vectorNorm);
        for (std::uint64_t rowIndex = 1; rowIndex < u.getNumberOfRows(); ++rowIndex) {
            matrix(rowIndex, 0) = vector(rowIndex, 0);
        }

        double uSquaredNorm = vectorSquaredNormWithoutFirstElement + std::norm(u.get(0, 0));
        uInvSquaredNorm = utils::isNull(uSquaredNorm)
            ? std::nullopt
            : std::make_optional(1 / uSquaredNorm);
    }

    Matrix::Value get(std::uint64_t i, std::uint64_t j) const override {
        assert(i < getNumberOfRows());
        assert(j < getNumberOfCols());
        std::complex<double> value = ((i == j) ? 1. : 0.);

        if (uInvSquaredNorm) {
            value -= 2 * uInvSquaredNorm.value() * u.get(i, 0) * std::conj(u.get(j, 0));
        }

        return value;
    }

private:
    std::optional<double> uInvSquaredNorm = 0.;
};


// Is this intended to be the opposite to a submatrix, e.g. a supermatrix or something like that?
class MatrixExtender : public Matrix {
public:
    MatrixExtender(Matrix const& m, std::uint64_t size)
        : Matrix(size, size)
        , matrix(m)
        , minRow(getNumberOfRows() - matrix.getNumberOfRows())
        , minCol(getNumberOfCols() - matrix.getNumberOfCols()) {

        assert(matrix.getNumberOfCols() == matrix.getNumberOfRows());
        assert(matrix.getNumberOfCols() <= size);
    }

    Value get(std::uint64_t i, std::uint64_t j) const override {
        if (i >= minRow && j >= minCol) {
            return matrix.get(i - minRow, j - minCol);
        }

        return i == j ? 1. : 0.;
    }

private:
    Matrix const& matrix;
    std::uint64_t const minRow = 0;
    std::uint64_t const minCol = 0;
};


inline void qrTriangularize(EditableAbstractMatrix& matrix, EditableAbstractMatrix* auxiliary = nullptr) {
    SubMatrixAdapter::Coordinates subMatrixCoordinates{
        .startRow = 0,
        .startCol = 0,
        .endRow = matrix.getNumberOfRows(),
        .endCol = matrix.getNumberOfCols()
    };
    while (subMatrixCoordinates.startRow < std::min(matrix.getNumberOfRows(), matrix.getNumberOfCols())) {
        assert(subMatrixCoordinates.startRow == subMatrixCoordinates.startCol);

        EditableSubMatrixAdapter subMatrix(matrix, subMatrixCoordinates);

        auto firstCol = extractCol(subMatrix, 0);
        HouseholderMatrix householderMatrix(firstCol);

        subMatrix.multiplyLeft(householderMatrix);

        if (auxiliary) {
            auxiliary->multiplyLeft(MatrixExtender(householderMatrix, auxiliary->getNumberOfRows()));
        }

        ++subMatrixCoordinates.startRow;
        ++subMatrixCoordinates.startCol;
    }

    assert(matrix.isUpperTriangular());
}


inline double computeSpectralRadius(Matrix const& matrix) {
    // Power iterations.
    
    assert(matrix.getNumberOfRows() == matrix.getNumberOfCols());
    static std::uint64_t const MAX_ITERATIONS = 1000;

    EditableMatrix v(matrix.getNumberOfCols(), 1);

    for (std::uint64_t i = 0; i < v.getNumberOfRows(); ++i) {
        v.set(i, 0, std::complex<double>{random::randomZeroOneDouble(), random::randomZeroOneDouble()});
    }

    std::uint64_t iteration = 0;
    double previousVNorm = 0.;
    while (true) {
        v.multiplyLeft(matrix);

        double vNorm = 0.;
        v.forEach([&vNorm](auto, auto, std::complex<double> v) {
            vNorm += std::norm(v);
        });
        vNorm = std::sqrt(vNorm);
        assert(utils::isNotNull(vNorm));

        if (iteration >= MAX_ITERATIONS) {
            return vNorm;
        } else if (utils::isNull(previousVNorm - vNorm)) {
            return vNorm;
        } else {
            ++iteration;
            previousVNorm = vNorm;
        }

        double vNormInv = 1 / vNorm;

        v *= vNormInv;

    }
}

}  // namespace qx