#pragma once

#include "qx/Compat.hpp"
#include "absl/container/btree_set.h"
#include "absl/container/btree_map.h"
#include "qx/Matrix.hpp" 
#include "qx/Common.hpp" 
#include <optional>
#include <cmath>

namespace qx {

template <std::uint64_t MaxNumberOfQubits = 64>
struct Column {
    using BasisVector = utils::BasisVector<MaxNumberOfQubits>;

    BasisVector measurementRegister;
    BasisVector basisVector;

    bool operator<(Column const& other) const {
        if (measurementRegister == other.measurementRegister) {
            return basisVector < other.basisVector;
        }
        return measurementRegister < other.measurementRegister;
    }
};

template <std::uint64_t MaxNumberOfQubits = 64>
class SparseStateAdapter : public EditableAbstractMatrix {
public:
    using UsedColumn = Column<MaxNumberOfQubits>;

    // I would define a using alias for this type
    // E.g.
    // using BTreeMapQubitDouble = absl::btree_map<Key<MaxNumberOfQubits>, std::complex<double>>
    //
    // Maybe turning this into a free function? E.g. make_sparse_state
    // Alike make_unique or make_shared
    static std::optional<SparseStateAdapter> makeIfNeeded(BTreeMapQubitDouble& data) {
        absl::btree_set<UsedColumn> usedColumns;
        for (auto const& [key, complexAmplitude]: data) {
            usedColumns.insert({key.measurementRegister, key.basisVector});
        }

        std::uint64_t const& dimension1 = usedColumns.size();
        std::uint64_t dimension2 = std::min(data.rbegin()->first.ensembleIndex.value + 1, dimension1);
        std::uint64_t maxSize = dimension2 * (2 * dimension1 - dimension2 + 1) / 2;

        if (data.size() <= maxSize) {
            return std::nullopt;
        }

        return SparseStateAdapter(std::move(usedColumns), data);
    }

    AbstractMatrix::Value get(std::uint64_t i, std::uint64_t j) const override {
        assert(i <= getNumberOfRows() && j <= getNumberOfCols());

        auto correspondingState = std::next(columns.begin(), j);

        auto it = data.find(Key<MaxNumberOfQubits>{ .ensembleIndex = i, .measurementRegister = correspondingState->measurementRegister, .basisVector = correspondingState->basisVector });
        if (it == data.end()) {
            return 0.;
        }

        return it->second;
    }

    void set(std::uint64_t i, std::uint64_t j, Value v) override {
        assert(i <= getNumberOfRows() && j <= getNumberOfCols());
        auto correspondingState = std::next(columns.begin(), j);

        // For clarity, I would  write this as a separate variable and a move
        // E.g.
        auto blah{
            .ensembleIndex = i,
            .measurementRegister = correspondingState->measurementRegister,
            .basisVector = correspondingState->basisVector
        };
        data[std::move(blah)] = v;
    }

    void forEach(std::function<void(std::uint64_t, std::uint64_t, AbstractMatrix::Value)> f) const override {
        for (auto const& [k, v]: data) {
            // For clarity, I would  write this as a separate variable and a move
            // E.g.
            auto blah = UsedColumn{ .measurementRegister = k.measurementRegister, .basisVector = k.basisVector };
            auto it = columns.find(std::move(blah));
            f(k.ensembleIndex.value, it - columns.begin(), v);
        }
    }

private:
    SparseStateAdapter(absl::btree_set<UsedColumn>&& c, BTreeMapQubitDouble& d)
        : AbstractMatrix(d.rbegin()->first.ensembleIndex.value + 1, c.size()), data(d), columns(c) {}

    absl::btree_map<Key<MaxNumberOfQubits>, std::complex<double>>& data;
    absl::btree_set<UsedColumn> const columns;
};

class MixedStateSimplifier {
public:
    template <std::uint64_t MaxNumberOfQubits = 64>
    static void simplify(BTreeMapQubitDouble& data) {
        auto stateMatrix = SparseStateAdapter<MaxNumberOfQubits>::makeIfNeeded(data);

        if (stateMatrix) {
            qrTriangularize(*stateMatrix);
        }
    }

    template <std::uint64_t MaxNumberOfQubits = 64>
    static void tidy(BTreeMapQubitDouble& data) {
        auto it = data.begin();
        std::optional<EnsembleIndex> currentEnsembleIndex{};
        EnsembleIndex targetEnsembleIndex{0};
        while (it != data.end()) {
            if (utils::isNull(it->second)) {
                it = data.extract_and_get_next(it).next;
                continue;
            }

            if (!currentEnsembleIndex) {
                currentEnsembleIndex = it->first.ensembleIndex;
            }

            if (it->first.ensembleIndex != currentEnsembleIndex) {
                currentEnsembleIndex = it->first.ensembleIndex;
                ++targetEnsembleIndex;
            }

            if (it->first.ensembleIndex != targetEnsembleIndex) {
                auto [node, next] = data.extract_and_get_next(it);
                node.key().ensembleIndex = targetEnsembleIndex;
                it = data.insert(next, std::move(node));
            }

            ++it;
        }
    }

    template <std::uint64_t MaxNumberOfQubits = 64>
    static void sparsifyGivens(BTreeMapQubitDouble& data) {
        while(sparsifyGivensImpl(data));
    }

    template <std::uint64_t MaxNumberOfQubits = 64>
    static bool sparsifyGivensImpl(BTreeMapQubitDouble& data) {
        tidy(data);

        if (data.rbegin()->first.ensembleIndex == EnsembleIndex{0}) {
            return false;
        }

        // I think this block should go into a different function
        using DataIt = BTreeMapQubitDouble::iterator;
        absl::btree_multimap<std::uint64_t, DataIt> numberOfNonZerosToDataRow;
        auto start = data.begin();
        while(start != data.end()) {
            EnsembleIndex currentEnsembleIndex = start->first.ensembleIndex;
            auto it = start;
            std::uint64_t numberOfNonZeros = 0;
            do {
                ++numberOfNonZeros;
                std::advance(it, 1);
            } while (it != data.end() && it->first.ensembleIndex == currentEnsembleIndex);
            
            numberOfNonZerosToDataRow.insert({numberOfNonZeros, start});

            start = it;
        }

        bool hasProgressed = false;

        // I think block this should go into a different function
        // And I would also try to split that other function into a few of other functions
        auto mainIt = numberOfNonZerosToDataRow.begin();
        while (mainIt != numberOfNonZerosToDataRow.end()) {
            auto otherCol = std::next(mainIt);
            while (otherCol != numberOfNonZerosToDataRow.end() && otherCol->first == mainIt->first) {
                bool haveSameShape = true;
                auto col1 = mainIt->second;
                auto col2 = otherCol->second;
                for (std::uint64_t i = 0; i < mainIt->first; ++i) {
                    assert(col1 != data.end());
                    assert(col2 != data.end());
                    
                    if (col1->first.measurementRegister != col2->first.measurementRegister ||
                        col1->first.basisVector != col2->first.basisVector) {
                        haveSameShape = false;
                        break;
                    }
                    std::advance(col1, 1);
                    std::advance(col2, 1);
                }

                if (!haveSameShape) {
                    std::advance(otherCol, 1);
                    continue;
                }

                hasProgressed = true;

                std::complex<double> const a = mainIt->second->second;
                std::complex<double> const b = otherCol->second->second;
                assert(std::hypot(std::abs(a), std::abs(b)) > config::EPS);
                double const invr = 1 / std::hypot(std::abs(a), std::abs(b));

                std::complex<double> const c = b * invr;
                std::complex<double> const s = a * invr;
                std::complex<double> const cConj = std::conj(c);
                std::complex<double> const sConj = std::conj(s);

                col1 = mainIt->second;
                col2 = otherCol->second;
                for (std::uint64_t i = 0; i < mainIt->first; ++i) {
                    assert(col1 != data.end());
                    assert(col2 != data.end());

                    auto const oldCol1 = col1->second;
                    col1->second = c * oldCol1 - s * col2->second;
                    col2->second = sConj * oldCol1 + cConj * col2->second;

                    assert(i != 0 || utils::isNull(col1->second));
                    std::advance(col1, 1);
                    std::advance(col2, 1);
                }

                break;
            }
            std::advance(mainIt, 1);
        }

        return hasProgressed;
    }


    template <std::uint64_t MaxNumberOfQubits = 64>
    static void sparsify(BTreeMapQubitDouble& data) {
        tidy(data);

        if (data.size() <= 2000) {
            return;
        }

        absl::btree_set<Column<MaxNumberOfQubits>> columns;
        for (auto const& [key, complexAmplitude]: data) {
            columns.insert({key.measurementRegister, key.basisVector});
        }

        absl::InlinedVector<std::complex<double>, 10> x(data.rbegin()->first.ensembleIndex.value + 1, 0.);

        auto minEnsembleIndex = data.begin()->first.ensembleIndex;
        assert(minEnsembleIndex == EnsembleIndex{0});

        // I think block this should go into a different function
        // And I would also try to split that other function into a few of other functions
        while(!columns.empty()) {
            auto endRow = data.rbegin()->first.ensembleIndex;

            if (minEnsembleIndex >= endRow) {
                return;
            }

            std::uint64_t const& dimension1 = columns.size();
            std::uint64_t dimension2 = std::min(endRow.value - minEnsembleIndex.value + 1, dimension1);
            std::uint64_t maxSize = dimension2 * (2 * dimension1 - dimension2 + 1) / 2;

            auto column = *columns.begin();

            // For clarity, I would  write this as a separate variable and a move
            // E.g.
            auto blah{
                .ensembleIndex = minEnsembleIndex,
                .measurementRegister = column.measurementRegister,
                .basisVector = column.basisVector
            };
            auto lb = data.lower_bound(std::move(blah));
            assert(data.end() - lb >= 0);
            if (static_cast<std::uint64_t>(std::abs(data.end() - lb)) <= maxSize) {
                return;
            }

            columns.erase(columns.begin());
            double xSquaredNormWithoutFirstElement = 0.;
            for (auto row = minEnsembleIndex; row <= endRow; ++row) {
                // For clarity, I would  write this as a separate variable and a move
                // E.g.
                auto blah{
                    .ensembleIndex = row,
                    .measurementRegister = column.measurementRegister,
                    .basisVector = column.basisVector
                };
                auto it = data.find(std::move(blah));
                if (it != data.end()) {
                    x[row.value] = it->second;
                    if (row != minEnsembleIndex) {
                        xSquaredNormWithoutFirstElement += std::norm(it->second);
                        data.erase(it);
                    }
                } else {
                    x[row.value] = 0.;
                }
            }

            if (utils::isNull(xSquaredNormWithoutFirstElement)) {
                ++minEnsembleIndex;
                continue;
            }

            // Would it be useful to have this as a separate function (useful meaning, this code is likely to be used somewhere else)?
            // Even if not, just for this function simplicity and readability, it would be a good idea, I think
            std::complex<double> alpha = [&]() -> std::complex<double> {
                if (x[minEnsembleIndex.value] == 0.) {
                    return std::sqrt(xSquaredNormWithoutFirstElement);
                }

                return x[minEnsembleIndex.value] /
                       std::abs(x[minEnsembleIndex.value]) *
                       std::sqrt(xSquaredNormWithoutFirstElement + std::norm(x[minEnsembleIndex.value]));
            }();

            // For clarity, I would  write this as a separate variable and a move
            // E.g.
            auto blah{
                .ensembleIndex = minEnsembleIndex,
                .measurementRegister = column.measurementRegister,
                .basisVector = column.basisVector
            };
            data[std::move(blah)] = -alpha;

            auto u = [&](EnsembleIndex i) {
                assert(i >= minEnsembleIndex);

                if (i == minEnsembleIndex) {
                    return x[minEnsembleIndex.value] + alpha;
                }

                return x[i.value];
            };

            auto uSquaredNorm = x[minEnsembleIndex.value] == 0.
                ? 2 * xSquaredNormWithoutFirstElement
                : (xSquaredNormWithoutFirstElement + std::norm(x[minEnsembleIndex.value] + alpha));

            assert(utils::isNotNull(uSquaredNorm));
            auto uInvSquaredNormTimesTwo = 2 / uSquaredNorm;

            // Since this is a function, I would rename it as an action, e.g. getHouseHolderMatrix
            auto householderMatrix = [&](EnsembleIndex i, EnsembleIndex j) -> std::complex<double> {
                assert(i >= minEnsembleIndex);
                assert(j >= minEnsembleIndex);
                return (i == j ? 1. : 0.) - uInvSquaredNormTimesTwo * u(i) * std::conj(u(j));                
            };

            // I would define a using alias for this type
            // E.g.
            // using BTreeMap = absl::btree_map<Key<MaxNumberOfQubits>, std::complex<double>>
            BTreeMap temp;
            
            // For clarity, I would  write this as a separate variable and a move
            // E.g.
            auto blah{
                .ensembleIndex = minEnsembleIndex,
                .measurementRegister = column.measurementRegister,
                .basisVector = column.basisVector
            };
            auto multiplyIt = std::next(data.find(std::move(blah)));

            while (multiplyIt != data.end()) {
                if (multiplyIt->first.measurementRegister == column.measurementRegister &&
                    multiplyIt->first.basisVector == column.basisVector) {

                    ++multiplyIt;
                    continue;
                }

                for (EnsembleIndex row = minEnsembleIndex; row <= endRow; ++row) {
                    auto added = householderMatrix(row, multiplyIt->first.ensembleIndex) * multiplyIt->second;
                    if (utils::isNotNull(added)) {
                        // For clarity, I would write this as a separate variable and a move
                        // E.g.
                        auto blah{
                            .ensembleIndex = row,
                            .measurementRegister = multiplyIt->first.measurementRegister,
                            .basisVector = multiplyIt->first.basisVector
                        };
                        temp[std::move(blah)] += added;
                    }
                }
                multiplyIt = data.erase(multiplyIt);
            }

            data.merge(std::move(temp));
            tidy(data);

            ++minEnsembleIndex;
        }
    }
};

}