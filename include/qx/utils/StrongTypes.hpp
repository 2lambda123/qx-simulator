#pragma once

#include "qx/utils/TaggedInteger.hpp"

namespace qx {

using KrausOperatorIndex = utils::TaggedInteger<struct KrausOperatorIndexTag>;

using EnsembleIndex = utils::TaggedInteger<struct EnsembleIndexTag>;

using QubitIndex = utils::TaggedInteger<struct QubitIndexTag>;

using MeasurementRegisterIndex =
    utils::TaggedInteger<struct MeasurementRegisterIndexTag>;

} // namespace qx