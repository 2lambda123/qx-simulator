#pragma once

#include "absl/container/btree_map.h"
#include "absl/container/btree_set.h"
#include "qx/Common.hpp"
#include "qx/Compat.hpp"
#include "qx/Matrix.hpp"
#include <cmath>
#include <optional>

namespace qx {

template <std::uint64_t MaxNumberOfQubits = 64> struct Column {
    using BasisVector = utils::BasisVector<MaxNumberOfQubits>;

    BasisVector measurementRegister;
    BasisVector basisVector;

    bool operator<(Column const &other) const {
        if (measurementRegister == other.measurementRegister) {
            return basisVector < other.basisVector;
        }
        return measurementRegister < other.measurementRegister;
    }
};

template <std::uint64_t MaxNumberOfQubits = 64>
class SparseStateAdapter : public EditableAbstractMatrix {
public:
    using MixedStateData = absl::btree_map<Key<MaxNumberOfQubits>, std::complex<double>>;
    using UsedColumn = Column<MaxNumberOfQubits>;

    static std::optional<SparseStateAdapter> makeIfNeeded(
        MixedStateData &data) {
        absl::btree_set<UsedColumn> usedColumns;
        for (auto const &[key, complexAmplitude] : data) {
            usedColumns.insert({key.measurementRegister, key.basisVector});
        }

        std::uint64_t const &dimension1 = usedColumns.size();
        std::uint64_t dimension2 =
            std::min(data.rbegin()->first.ensembleIndex.value + 1, dimension1);
        std::uint64_t maxSize =
            dimension2 * (2 * dimension1 - dimension2 + 1) / 2;

        if (data.size() <= maxSize) {
            return std::nullopt;
        }

        return SparseStateAdapter(std::move(usedColumns), data);
    }

    AbstractMatrix::Value get(std::uint64_t i, std::uint64_t j) const override {
        assert(i <= getNumberOfRows() && j <= getNumberOfCols());

        auto correspondingState = std::next(columns.begin(), j);

        auto it = data.find(Key<MaxNumberOfQubits>{
            .ensembleIndex = {i},
            .measurementRegister = correspondingState->measurementRegister,
            .basisVector = correspondingState->basisVector});
        if (it == data.end()) {
            return 0.;
        }

        return it->second;
    }

    void set(std::uint64_t i, std::uint64_t j, Value v) override {
        assert(i <= getNumberOfRows() && j <= getNumberOfCols());
        auto correspondingState = std::next(columns.begin(), j);

        data[{.ensembleIndex = {i},
              .measurementRegister = correspondingState->measurementRegister,
              .basisVector = correspondingState->basisVector}] = v;
    }

    void forEach(
        std::function<void(std::uint64_t, std::uint64_t, AbstractMatrix::Value)>
            f) const override {
        for (auto const &[k, v] : data) {
            auto it = columns.find(
                UsedColumn{.measurementRegister = k.measurementRegister,
                           .basisVector = k.basisVector});
            f(k.ensembleIndex.value, it - columns.begin(), v);
        }
    }

private:
    SparseStateAdapter(
        absl::btree_set<UsedColumn> &&c,
        MixedStateData &d)
        : AbstractMatrix(d.rbegin()->first.ensembleIndex.value + 1, c.size()),
          data(d), columns(c) {}

    MixedStateData &data;
    absl::btree_set<UsedColumn> const columns;
};

class MixedStateSimplifier {
public:
    template <std::uint64_t MaxNumberOfQubits = 64>
    using MixedStateData = typename SparseStateAdapter<MaxNumberOfQubits>::MixedStateData;

    template <std::uint64_t MaxNumberOfQubits = 64>
    static void simplify(
        MixedStateData<MaxNumberOfQubits> &data) {
        auto stateMatrix =
            SparseStateAdapter<MaxNumberOfQubits>::makeIfNeeded(data);

        if (stateMatrix) {
            qrTriangularize(*stateMatrix);
        }
    }

    template <std::uint64_t MaxNumberOfQubits = 64>
    static void
    tidy(MixedStateData<MaxNumberOfQubits> &data) {
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
    static void sparsifyGivens(
        MixedStateData<MaxNumberOfQubits> &data) {
        while (sparsifyGivensImpl<MaxNumberOfQubits>(data));
    }

    template <std::uint64_t MaxNumberOfQubits = 64>
    static bool sparsifyGivensImpl(
        MixedStateData<MaxNumberOfQubits> &data) {
        tidy<MaxNumberOfQubits>(data);

        if (data.rbegin()->first.ensembleIndex == EnsembleIndex{0}) {
            return false;
        }

        using DataIt = typename absl::btree_map<Key<MaxNumberOfQubits>,
                                       std::complex<double>>::iterator;
        absl::btree_multimap<std::uint64_t, DataIt> numberOfNonZerosToDataRow;
        auto start = data.begin();
        while (start != data.end()) {
            EnsembleIndex currentEnsembleIndex = start->first.ensembleIndex;
            auto it = start;
            std::uint64_t numberOfNonZeros = 0;
            do {
                ++numberOfNonZeros;
                std::advance(it, 1);
            } while (it != data.end() &&
                     it->first.ensembleIndex == currentEnsembleIndex);

            numberOfNonZerosToDataRow.insert({numberOfNonZeros, start});

            start = it;
        }

        bool hasProgressed = false;

        auto mainIt = numberOfNonZerosToDataRow.begin();
        while (mainIt != numberOfNonZerosToDataRow.end()) {
            auto otherCol = std::next(mainIt);
            while (otherCol != numberOfNonZerosToDataRow.end() &&
                   otherCol->first == mainIt->first) {
                bool haveSameShape = true;
                auto col1 = mainIt->second;
                auto col2 = otherCol->second;
                for (std::uint64_t i = 0; i < mainIt->first; ++i) {
                    assert(col1 != data.end());
                    assert(col2 != data.end());

                    if (col1->first.measurementRegister !=
                            col2->first.measurementRegister ||
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
                assert(std::hypot(std::abs(a), std::abs(b)) > config::ATOL);
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
    static void sparsify(
        MixedStateData<MaxNumberOfQubits> &data) {
        tidy(data);

        if (data.size() <= 2000) {
            return;
        }

        absl::btree_set<Column<MaxNumberOfQubits>> columns;
        for (auto const &[key, complexAmplitude] : data) {
            columns.insert({key.measurementRegister, key.basisVector});
        }

        absl::InlinedVector<std::complex<double>, 10> x(
            data.rbegin()->first.ensembleIndex.value + 1, 0.);

        auto minEnsembleIndex = data.begin()->first.ensembleIndex;
        assert(minEnsembleIndex == EnsembleIndex{0});

        while (!columns.empty()) {
            auto endRow = data.rbegin()->first.ensembleIndex;

            if (minEnsembleIndex >= endRow) {
                return;
            }

            std::uint64_t const &dimension1 = columns.size();
            std::uint64_t dimension2 =
                std::min(endRow.value - minEnsembleIndex.value + 1, dimension1);
            std::uint64_t maxSize =
                dimension2 * (2 * dimension1 - dimension2 + 1) / 2;

            auto column = *columns.begin();

            auto lb = data.lower_bound(
                {.ensembleIndex = minEnsembleIndex,
                 .measurementRegister = column.measurementRegister,
                 .basisVector = column.basisVector});
            assert(data.end() - lb >= 0);
            if (static_cast<std::uint64_t>(std::abs(data.end() - lb)) <=
                maxSize) {
                return;
            }

            columns.erase(columns.begin());
            double xSquaredNormWithoutFirstElement = 0.;
            for (auto row = minEnsembleIndex; row <= endRow; ++row) {
                auto it = data.find(
                    {.ensembleIndex = row,
                     .measurementRegister = column.measurementRegister,
                     .basisVector = column.basisVector});
                if (it != data.end()) {
                    x[row.value] = it->second;
                    if (row != minEnsembleIndex) {
                        xSquaredNormWithoutFirstElement +=
                            std::norm(it->second);
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

            std::complex<double> alpha = [&]() -> std::complex<double> {
                if (x[minEnsembleIndex.value] == 0.) {
                    return std::sqrt(xSquaredNormWithoutFirstElement);
                }

                return x[minEnsembleIndex.value] /
                       std::abs(x[minEnsembleIndex.value]) *
                       std::sqrt(xSquaredNormWithoutFirstElement +
                                 std::norm(x[minEnsembleIndex.value]));
            }();

            data[{.ensembleIndex = minEnsembleIndex,
                  .measurementRegister = column.measurementRegister,
                  .basisVector = column.basisVector}] = -alpha;

            auto u = [&](EnsembleIndex i) {
                assert(i >= minEnsembleIndex);

                if (i == minEnsembleIndex) {
                    return x[minEnsembleIndex.value] + alpha;
                }

                return x[i.value];
            };

            auto uSquaredNorm =
                x[minEnsembleIndex.value] == 0.
                    ? 2 * xSquaredNormWithoutFirstElement
                    : (xSquaredNormWithoutFirstElement +
                       std::norm(x[minEnsembleIndex.value] + alpha));

            assert(utils::isNotNull(uSquaredNorm));
            auto uInvSquaredNormTimesTwo = 2 / uSquaredNorm;

            auto householderMatrix =
                [&](EnsembleIndex i, EnsembleIndex j) -> std::complex<double> {
                assert(i >= minEnsembleIndex);
                assert(j >= minEnsembleIndex);
                return (i == j ? 1. : 0.) -
                       uInvSquaredNormTimesTwo * u(i) * std::conj(u(j));
            };

            MixedStateData<MaxNumberOfQubits> temp;

            auto multiplyIt = std::next(
                data.find({.ensembleIndex = minEnsembleIndex,
                           .measurementRegister = column.measurementRegister,
                           .basisVector = column.basisVector}));

            while (multiplyIt != data.end()) {
                if (multiplyIt->first.measurementRegister ==
                        column.measurementRegister &&
                    multiplyIt->first.basisVector == column.basisVector) {
                    ++multiplyIt;
                    continue;
                }

                for (EnsembleIndex row = minEnsembleIndex; row <= endRow;
                     ++row) {
                    auto added = householderMatrix(
                                     row, multiplyIt->first.ensembleIndex) *
                                 multiplyIt->second;
                    if (utils::isNotNull(added)) {
                        temp[{.ensembleIndex = row,
                              .measurementRegister =
                                  multiplyIt->first.measurementRegister,
                              .basisVector = multiplyIt->first.basisVector}] +=
                            added;
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

} // namespace qx