#pragma once

#include <complex>
#include <optional>

#include "absl/container/inlined_vector.h"
#include "qx/CompileTimeConfiguration.hpp"
#include "qx/utils/FloatComparison.hpp"
#include "qx/Random.hpp"

namespace qx {

// General comments
// * forEach, map, and so on, may not be needed if Matrix would implement begin/end/cbegin/cend...
//   - This is a bit complex because it requires the classes iterator/const_iterator/maybe reverse iterator to be implemented
//   - But it would open the whole of STL algorithms (or maybe some of them, e.g. the ones that use forward iterators) to be available to be used with Matrix
//   - And it would render the set of Editable classes (and the multiple inheritance) as unnecessary (non-const matrices would be editable)
//   - Also, it may render AbstractMatrix as unnecessary as well (just merge it with Matrix)
// * Why is HouseholdMatrix not inheriting from UnitaryMatrix?
// * I've seen there are some C++ linear algebra libraries implementing dense matrices (e.g. armadillo, Boost uBLAS)
//   - Would it be overkill to use them here? Or maybe not really suited for what we want to do?

class AbstractMatrix {
public:
    using Value = std::complex<double>;

    AbstractMatrix(std::uint64_t n, std::uint64_t m) : numberOfRows(n), numberOfCols(m) {}

    explicit AbstractMatrix(std::uint64_t n) : numberOfRows(n), numberOfCols(n) {}

    virtual ~AbstractMatrix() = default;

    std::uint64_t getNumberOfRows() const {
        return numberOfRows;
    }

    std::uint64_t getNumberOfCols() const {
        return numberOfCols;
    }

    virtual Value get(std::uint64_t i, std::uint64_t j) const = 0;

    // I think the proper way to do this would be to have begin/end/++ methods in the matrix
    // Then you could feed a matrix to any std algorithm
    /*
    virtual void forEach(std::function<void(std::uint64_t, std::uint64_t, AbstractMatrix::Value)> f) const { // FIXME: find a way to avoid std::function
        for (std::uint64_t rowIndex = 0; rowIndex < getNumberOfRows(); ++rowIndex) {
            for (std::uint64_t colIndex = 0; colIndex < getNumberOfCols(); ++colIndex) {
                f(rowIndex, colIndex, get(rowIndex, colIndex));
            }
        }
    }
    */

    [[nodiscard]] virtual constexpr iterator begin() noexcept = 0;
    [[nodiscard]] virtual constexpr iterator end() noexcept = 0;
    [[nodiscard]] virtual constexpr const_iterator begin() const noexcept = 0;
    [[nodiscard]] virtual constexpr const_iterator end() const noexcept = 0;
    [[nodiscard]] virtual constexpr const_iterator cbegin() const noexcept = 0;
    [[nodiscard]] virtual constexpr const_iterator cend() const noexcept = 0;

    virtual bool operator==(AbstractMatrix const &other) const {
        return (getNumberOfRows() == other.getNumberOfRows() &&
                getNumberOfCols() == other.getNumberOfCols() &&
                equal(getData(), other.getData())
            );
    }

    inline void print(std::ostream &os) const {
        for (std::uint64_t i = 0; i < getNumberOfRows(); ++i) {
            bool first = true;
            for (std::uint64_t j = 0; j < getNumberOfCols(); ++j) {
                if (!first) {
                    os << "  ";
                } else {
                    first = false;
                }

                os << get(i, j);
            }
        }
    }
    
    bool isUpperTriangular() {
        bool result = true;
        forEach([&result](std::uint64_t rowIndex, std::uint64_t colIndex, auto value) {
            if (rowIndex > colIndex && utils::isNotNull(value)) {
                result = false;
            }
        });

        return result;
    }

private:
    std::uint64_t const numberOfRows = 0;
    std::uint64_t const numberOfCols = 0;
};


class SubMatrixAdapter : public virtual AbstractMatrix {
public:
    struct Coordinates {
        std::uint64_t startRow = 0;
        std::uint64_t startCol = 0;
        std::uint64_t endRow = 0;
        std::uint64_t endCol = 0;
    };

    SubMatrixAdapter(AbstractMatrix& m, Coordinates const& c)
        : AbstractMatrix(c.endRow - c.startRow, c.endCol - c.startCol), matrix(m), coordinates(c) {

        assert(coordinates.startRow <= coordinates.endRow);
        assert(coordinates.startCol <= coordinates.endCol);
        assert(coordinates.endRow <= matrix.getNumberOfRows());
        assert(coordinates.endCol <= matrix.getNumberOfCols());
    }

    AbstractMatrix::Value get(std::uint64_t i, std::uint64_t j) const override {
        assert(i < getNumberOfRows());
        assert(j < getNumberOfCols());
        return matrix.get(i + coordinates.startRow, j + coordinates.startCol);
    }
    
/*
    void forEach(std::function<void(std::uint64_t, std::uint64_t, AbstractMatrix::Value)> f) const override {
        matrix.forEach([&](std::uint64_t rowIndex, std::uint64_t colIndex, AbstractMatrix::Value value) {
            if (rowIndex >= coordinates.startRow && rowIndex < coordinates.endRow // FIXME: not efficient. Maybe forEach should get begin and end indices optionally? Iteration should always be ordered.
                    && colIndex >= coordinates.startCol && colIndex < coordinates.endCol) {
                        f(rowIndex - coordinates.startRow, colIndex - coordinates.startCol, value);
            }
        });
    }
*/

protected:
    AbstractMatrix& matrix;
    Coordinates coordinates;
};

inline SubMatrixAdapter extractRow(AbstractMatrix& m, std::uint64_t rowIndex) {
    assert(rowIndex < m.getNumberOfRows());
    SubMatrixAdapter::Coordinates rowCoordinates {
        .startRow = rowIndex,
        .startCol = 0,
        .endRow = rowIndex + 1,
        .endCol = m.getNumberOfCols()
    };
    return SubMatrixAdapter(m, rowCoordinates);
}

inline SubMatrixAdapter extractCol(AbstractMatrix& m, std::uint64_t colIndex) {
    assert(colIndex < m.getNumberOfCols());
    SubMatrixAdapter::Coordinates colCoordinates {
        .startRow = 0,
        .startCol = colIndex,
        .endRow = m.getNumberOfRows(),
        .endCol = colIndex + 1
    };
    return SubMatrixAdapter(m, colCoordinates);
}


// If 1) you can declare a matrix as const or non-const,
// and 2) you provide const and non-const iterators for a matrix,
// then you could get rid of all the editable matrices in the hierarchy, couldn't you?
class EditableAbstractMatrix : public virtual AbstractMatrix {
public:
    virtual void set(std::uint64_t, std::uint64_t, AbstractMatrix::Value) = 0;

    void add(std::uint64_t i, std::uint64_t j, AbstractMatrix::Value v) {
        auto current = get(i, j);
        set(i, j, current + v);
    }

    /*
    virtual void map(std::function<std::optional<AbstractMatrix::Value>(std::uint64_t, std::uint64_t, AbstractMatrix::Value)> f) { // FIXME: find a way to avoid std::function
        forEach([&](std::uint64_t rowIndex, std::uint64_t colIndex, AbstractMatrix::Value v) {
            auto result = f(rowIndex, colIndex, v);
            if (result) {
                set(rowIndex, colIndex, *result);
            }
        });
    }
    */

    virtual void operator*=(AbstractMatrix const& other) {
        assert(getNumberOfCols() == other.getNumberOfRows() &&
               other.getNumberOfRows() == other.getNumberOfCols());

        for (std::uint64_t rowIndex = 0; rowIndex < getNumberOfRows(); ++rowIndex) {
            // I guess the compiler reuses newRow's memory in every iteration
            // apart from resetting all of its elements to 0.?
            std::vector<AbstractMatrix::Value> newRow(getNumberOfCols(), 0.);

            auto row = extractRow(*this, rowIndex);
            auto col = extractCol(other, rowIndex);
            std::ranges::transform(row, col, newRow,  // newRow = row * col
                                   [&newRow](auto e1, auto e2) { return e1 * e2; });
            /*
            row.forEach([&](std::uint64_t, std::uint64_t colIndex, auto value) {
                for (std::uint64_t otherColIndex = 0; otherColIndex < other.getNumberOfCols(); ++otherColIndex) {
                    newRow[otherColIndex] += value * other.get(colIndex, otherColIndex);
                }
            });
            */

            std::ranges::copy(newRow, row);
            /*
            for (std::uint64_t colIndex = 0; colIndex < getNumberOfCols(); ++colIndex) {
                set(rowIndex, colIndex, newRow[colIndex]);
            }
            */
        }
    }
    
    virtual void operator*=(double d) {
        std::ranges::for_each(*this, [&d](auto e) { e *= d; });
        /*
        map([d](auto, auto, std::complex<double> v) {
            return v * d;
        });
        */
    }
    
    void multiplyLeft(AbstractMatrix const& other) {
        assert(other.getNumberOfCols() == getNumberOfRows() &&
               other.getNumberOfRows() == other.getNumberOfCols());

        std::vector<AbstractMatrix::Value> newCol(getNumberOfRows(), 0);
        for (std::uint64_t colIndex = 0; colIndex < getNumberOfCols(); ++colIndex) {
            for (auto& x: newCol) {
                x = 0.;
            }

            auto col = extractCol(*this, colIndex);
            col.forEach([&](std::uint64_t rowIndex, std::uint64_t, auto value) {
                for (std::uint64_t otherRowIndex = 0; otherRowIndex < other.getNumberOfRows(); ++otherRowIndex) {
                    newCol[otherRowIndex] += value * other.get(otherRowIndex, rowIndex);
                }
            });

            for (std::uint64_t rowIndex = 0; rowIndex < getNumberOfRows(); ++rowIndex) {
                set(rowIndex, colIndex, newCol[rowIndex]);
            }
        }
    }
    
    virtual void operator+=(AbstractMatrix const &other) {
        assert(getNumberOfRows() == other.getNumberOfRows() &&
               getNumberOfCols() == other.getNumberOfCols() &&
               "Can't add matrices of different sizes");

        std::transform(this->begin(), this->end(), other->begin(),
                       [](auto e1, auto e2) { return e1 + e2; });
        for (std::uint64_t i = 0; i < getNumberOfRows(); ++i) { // FIXME: forEach?
            for (std::uint64_t j = 0; j < getNumberOfCols(); ++j) {
                add(i, j, other.get(i, j));
            }
        }
    }
};


class EditableSubMatrixAdapter : public SubMatrixAdapter, public EditableAbstractMatrix {
public:
    using SubMatrixAdapter::Coordinates;
    EditableSubMatrixAdapter(EditableAbstractMatrix& m, Coordinates const& c)
        : AbstractMatrix(c.endRow - c.startRow, c.endCol - c.startCol), SubMatrixAdapter(m, c) {}

    void set(std::uint64_t i, std::uint64_t j, AbstractMatrix::Value v) override {
        assert(i < getNumberOfRows());
        assert(j < getNumberOfCols());
        return dynamic_cast<EditableAbstractMatrix&>(matrix).set(i + coordinates.startRow, j + coordinates.startCol, v);
    }

    void map(std::function<std::optional<Value>(std::uint64_t, std::uint64_t, Value)> f) override {
        dynamic_cast<EditableAbstractMatrix&>(matrix).map([&](std::uint64_t rowIndex, std::uint64_t colIndex, Value value) {
            if (rowIndex >= coordinates.startRow && rowIndex < coordinates.endRow
                    && colIndex >= coordinates.startCol && colIndex < coordinates.endCol) {
                        return f(rowIndex - coordinates.startRow, colIndex - coordinates.startCol, value);
            }
            return std::optional<Value>();
        });
    }
};

inline std::ostream& operator<<(std::ostream &os, AbstractMatrix const &m) {
    m.print(os);
    return os;
}


class Matrix : public virtual AbstractMatrix {
public:
    Matrix(std::uint64_t n, std::uint64_t m)
        : AbstractMatrix(n, m), matrix(getNumberOfRows() * getNumberOfCols(), 0.) {}

    Matrix(std::uint64_t n)
        : AbstractMatrix(n), matrix(getNumberOfRows() * getNumberOfCols(), 0.) {}

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

    explicit Matrix(AbstractMatrix const& m) : Matrix(m.getNumberOfRows(), m.getNumberOfCols()) {
        m.forEach([&](std::uint64_t rowIndex, std::uint64_t colIndex, auto value) {
            matrix[rowIndex * getNumberOfCols() + colIndex] = value;
        });
    }

    // You could overload operator()
    // https://isocpp.org/wiki/faq/operator-overloading#matrix-subscript-op
    // And access elements as m(i, j) instead of m.get(i, j)
    inline Value get(std::uint64_t i, std::uint64_t j) const override {
        assert(i < getNumberOfRows());
        assert(j < getNumberOfCols());
        assert(i * getNumberOfCols() + j < getNumberOfRows() * getNumberOfCols());
        return matrix[i * getNumberOfCols() + j];
    }

    void forEach(std::function<void(std::uint64_t, std::uint64_t, Value)> f) const override {
        std::uint64_t i = 0;
        for (auto const& v: matrix) {
            f(i / getNumberOfCols(), i % getNumberOfCols(), v);
            ++i;
        }
    }

    Matrix dagger() const {
        Matrix m(getNumberOfCols(), getNumberOfRows());
        for (std::uint64_t i = 0; i < getNumberOfCols(); ++i) {
            for (std::uint64_t j = 0; j < getNumberOfRows(); ++j) {
                m.matrix[i * getNumberOfCols() + j] = std::conj(get(j, i));
            }
        }

        return m;
    }

    Matrix operator*(Matrix const &other) const {
        assert(getNumberOfCols() == other.getNumberOfRows() && "Can't multiply matrices with incompatible sizes");
        Matrix m(getNumberOfRows(), other.getNumberOfCols());
        for (std::uint64_t i = 0; i < getNumberOfRows(); ++i) {
            for (std::uint64_t j = 0; j < other.getNumberOfCols(); ++j) {
                for (std::uint64_t k = 0; k < getNumberOfCols(); ++k) {
                    m.matrix[i * m.getNumberOfCols() + j] += get(i, k) * other.get(k, j);
                }
            }
        }

        return m;
    }

    friend inline Matrix operator*(double d, Matrix m);

protected:
    absl::InlinedVector<std::complex<double>, config::MAX_INLINED_MATRIX_ELEMENTS> matrix;
};

inline Matrix operator*(double d, Matrix m) {
    for (auto& v: m.matrix) {
        v *= d;
    }

    return m;
}


class EditableMatrix : public virtual Matrix, public EditableAbstractMatrix {
public:
    EditableMatrix(std::uint64_t n, std::uint64_t m) : AbstractMatrix(n, m), Matrix(n, m) {}

    EditableMatrix(std::uint64_t n) : AbstractMatrix(n), Matrix(n) {}

    EditableMatrix(AbstractMatrix const& m) : AbstractMatrix(m.getNumberOfRows(), m.getNumberOfCols()), Matrix(m) {}

    explicit EditableMatrix(std::initializer_list<std::initializer_list<std::complex<double>>> init)
        : AbstractMatrix(init.size(), init.begin() == init.end() ? 0 : init.begin()->size()), Matrix(init) {};

    inline void set(std::uint64_t i, std::uint64_t j, Value v) override {
        assert(i < getNumberOfRows());
        assert(j < getNumberOfCols());
        assert(i * getNumberOfCols() + j < getNumberOfRows() * getNumberOfCols());
        matrix[i * getNumberOfCols() + j] = v;
    }
};


class SquareMatrix : public virtual Matrix {
public:
    static SquareMatrix identity(std::uint64_t size) {
        SquareMatrix m(size);
        for (std::uint64_t i = 0; i < size; ++i) {
            for (std::uint64_t j = 0; j < size; ++j) {
               m.matrix[i * m.getNumberOfCols() + j] = (i == j) ? 1 : 0;
            }
        }

        return m;
    }

    explicit SquareMatrix(std::uint64_t s) : AbstractMatrix(s), Matrix(s) {};

    SquareMatrix(Matrix const& m) : AbstractMatrix(m.getNumberOfRows()), Matrix(m) {
        if (m.getNumberOfRows() != m.getNumberOfCols()) {
            throw std::runtime_error("Not a square matrix");
        }
    } 

    explicit SquareMatrix(std::initializer_list<std::initializer_list<std::complex<double>>> init)
        : AbstractMatrix(init.size()), Matrix(init) {

        if (init.size() != 0 && init.size() != init.begin()->size()) {
            throw std::runtime_error("Not a square matrix");
        }
    }

    std::uint64_t getSize() const {
        assert(getNumberOfRows() == getNumberOfCols());
        return getNumberOfRows();
    }
};


class EditableSquareMatrix : public SquareMatrix, public EditableMatrix {
public:
    static EditableSquareMatrix identity(std::uint64_t size) {
        EditableSquareMatrix m(size);
        for (std::uint64_t i = 0; i < size; ++i) {
            for (std::uint64_t j = 0; j < size; ++j) {
               m.set(i, j, (i == j) ? 1 : 0);
            }
        }

        return m;
    }

    explicit EditableSquareMatrix(std::uint64_t s)
        : AbstractMatrix(s), Matrix(s), SquareMatrix(s), EditableMatrix(s) {};

    explicit EditableSquareMatrix(std::initializer_list<std::initializer_list<std::complex<double>>> init)
        : AbstractMatrix(init.size()), Matrix(init), SquareMatrix(init), EditableMatrix(init) {}
};



class UnitaryMatrix : public SquareMatrix {
public:
    explicit UnitaryMatrix(std::initializer_list<std::initializer_list<std::complex<double>>> init)
        : AbstractMatrix(init.size(), init.size() == 0 ? 0 : init.begin()->size()), Matrix(init), SquareMatrix(init) {
            if (*this * SquareMatrix::dagger() !) SquareMatrix::identity(getSize())) {
                throw std::runtime_error("Matrix is not unitary");
            };
    }

    UnitaryMatrix(SquareMatrix const &m)
        : AbstractMatrix(m.getNumberOfRows(), m.getNumberOfCols()), Matrix(m), SquareMatrix(m) {
            if (*this * SquareMatrix::dagger() != SquareMatrix::identity(getSize())) {
                throw std::runtime_error("Matrix is not unitary");
            };
    }

    UnitaryMatrix(Matrix const &m)
        : AbstractMatrix(m.getNumberOfRows(), m.getNumberOfCols()), Matrix(m), SquareMatrix(m) {
            if (*this * SquareMatrix::dagger() != SquareMatrix::identity(getSize())) {
                throw std::runtime_error("Matrix is not unitary");
            };
    }
};


class HouseholderMatrix : public AbstractMatrix {
    // is also unitary. ---> shouldn't then inherit from UnitaryMatrix?
public:
    HouseholderMatrix(AbstractMatrix const& vector)
        : AbstractMatrix(vector.getNumberOfRows(), vector.getNumberOfRows())
        , u(vector.getNumberOfRows(), 1) {

        assert(vector.getNumberOfCols() == 1);

        double vectorSquaredNormWithoutFirstElement = 0.;
        vector.forEach([&vectorSquaredNormWithoutFirstElement](std::uint64_t rowIndex, std::uint64_t colIndex, AbstractMatrix::Value value) {
            assert(colIndex == 0);
            if (rowIndex == 0) {
                return;
            }
            vectorSquaredNormWithoutFirstElement += std::norm(value);
        });

        double vectorNorm = std::sqrt(vectorSquaredNormWithoutFirstElement + std::norm(vector.get(0, 0)));
        u.set(0, 0, vector.get(0, 0) + vector.get(0, 0) / std::abs(vector.get(0, 0)) * vectorNorm);
        for (std::uint64_t rowIndex = 1; rowIndex < u.getNumberOfRows(); ++rowIndex) {
            u.set(rowIndex, 0, vector.get(rowIndex, 0));
        }

        double uSquaredNorm = vectorSquaredNormWithoutFirstElement + std::norm(u.get(0, 0));
        uInvSquaredNorm = utils::isNull(uSquaredNorm) ? std::nullopt : std::make_optional(1 / uSquaredNorm);
    }

    AbstractMatrix::Value get(std::uint64_t i, std::uint64_t j) const override {
        assert(i < getNumberOfRows());
        assert(j < getNumberOfCols());
        std::complex<double> value = ((i == j) ? 1. : 0.);

        if (uInvSquaredNorm) {
            value -= 2 * uInvSquaredNorm.value() * u.get(i, 0) * std::conj(u.get(j, 0));
        }

        return value;
    }

private:
    EditableMatrix u;
    std::optional<double> uInvSquaredNorm = 0.;
};


// Is this intended to be the opposite to a submatrix, e.g. a supermatrix or something like that?
class MatrixExtender : public AbstractMatrix {
public:
    MatrixExtender(AbstractMatrix const& m, std::uint64_t size)
        : AbstractMatrix(size, size)
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
    AbstractMatrix const& matrix;
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


inline double computeSpectralRadius(AbstractMatrix const& matrix) {
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

}