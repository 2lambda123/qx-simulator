#pragma once

#include <optional>
#include <Python.h>
#include <complex.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include "qx/Core.hpp"
#include "qx/Operations.hpp"

namespace qxelarator {

std::optional<qx::Operations::Signature> getSignatureFromPython(PyObject* pObject) {
    if (!pObject) {
        PyErr_SetString(PyExc_ValueError, "Internal error");
        return std::nullopt;
    }

    if (!PyUnicode_Check(pObject)) {
        PyErr_SetString(PyExc_ValueError, "Operation signature must be a string");
        return std::nullopt;
    }

    auto signatureLength = PyUnicode_GET_LENGTH(pObject);
    if (signatureLength == 0) {
        PyErr_SetString(PyExc_ValueError, "Operation signature cannot be empty");
        return std::nullopt;
    }

    qx::Operations::Signature cppSignature;

    auto kind = PyUnicode_KIND(pObject);
    auto* data = PyUnicode_DATA(pObject);
    for (auto index = 0; index < signatureLength; ++index) {
        auto c = Py_UNICODE_TOLOWER(PyUnicode_READ(kind, data, index));

        using T = qx::Operations::OperandType;

        switch (c) {
            case 'q':
                cppSignature.push_back(T::Qubit);
                break;
            case 'b':
                cppSignature.push_back(T::ClassicalBit);
                break;
            case 'i':
                cppSignature.push_back(T::Int);
                break;
            case 'd':
                cppSignature.push_back(T::Double);
                break;
            default:
                PyErr_SetString(PyExc_ValueError, "Unknown char in signature, must be 'q', 'b', 'i' or 'd'");
                return std::nullopt;
        }
    }

    return cppSignature;
}

std::string getStaticOperandsFormat(qx::Operations::Signature signature) {
    std::stringstream format;
    using T = qx::Operations::OperandType;
    for (auto t: signature) {
        if (t == T::Double) {
            format << "d";
        } else if (t == T::Int) {
            format << "l";
        }
    }

    return format.str();
}

std::optional<std::string> getPythonString(PyObject *o) {
    if (!o) {
        return std::nullopt;
    }
    
    auto* pNameAsciiBytes = PyUnicode_AsASCIIString(o);
    if (!pNameAsciiBytes) {
        PyErr_SetString(PyExc_ValueError, "Operation name is not an ASCII string");
        return std::nullopt;
    }

    char* pNameAscii = PyBytes_AsString(pNameAsciiBytes); // FIXME: check
    Py_DECREF(pNameAsciiBytes);

    if (!pNameAscii) {
        PyErr_SetString(PyExc_ValueError, "Internal error with string");
        return std::nullopt;
    }
    std::string name(pNameAscii);

    return name;
}

qx::SquareMatrix numpyArrayToQXSquareMatrix(PyObject *o) {
    // Assumes o is a numpy.ndarray!

    auto* shape = PyObject_GetAttrString(o, "shape");

    auto nRows = PyLong_AsLong(PyTuple_GetItem(shape, 0));
    auto nCols = PyLong_AsLong(PyTuple_GetItem(shape, 1));

    if (nRows != nCols) {
        throw std::runtime_error("Not a square matrix");
    }

    qx::EditableSquareMatrix squareMatrix(nRows);
    for (std::uint64_t rowIndex = 0; rowIndex < nRows; ++rowIndex) {
        for (std::uint64_t colIndex = 0; colIndex < nCols; ++colIndex) {
            auto* matrixIndex = PyTuple_New(2);
            PyTuple_SetItem(matrixIndex, 0, PyLong_FromLong(rowIndex));
            PyTuple_SetItem(matrixIndex, 1, PyLong_FromLong(colIndex));
            auto* value = PyObject_GetItem(o, matrixIndex);

            auto valueComplex = PyComplex_AsCComplex(value); // FIXME: check

            squareMatrix.set(rowIndex, colIndex, std::complex{valueComplex.real, valueComplex.imag});
        }
    }

    return squareMatrix;
}

PyObject* getPythonSimulationError(char const* errorMessage) {
    auto* pQxelarator = PyImport_ImportModule("qxelarator");
    if (!pQxelarator) {
        return NULL;
    }
    
    auto* pSimulationError = PyObject_GetAttrString(pQxelarator, "SimulationError");
    Py_DECREF(pQxelarator);

    if (!pSimulationError) {
        return NULL;
    }

    auto* errorString = PyUnicode_FromString(errorMessage);

    if (!errorString) {
        return NULL;
    }

    auto* args = PyTuple_Pack(1, errorString);
    Py_DECREF(errorString);

    if (!args) {
        return NULL;
    }

    auto* pSimulationErrorInstance = PyObject_CallObject(pSimulationError, args);
    Py_DECREF(args);

    if (!pSimulationErrorInstance) {
        return NULL;
    }
    Py_DECREF(pSimulationError);

    return pSimulationErrorInstance;
}

// Same comment as below, could be a numpy structured array?
class PythonDictStringComplexMap : public qx::core::StringMap<std::complex<double>> {
public:
    PythonDictStringComplexMap() {
        underlying = PyDict_New();
        if (!underlying) {
            throw std::runtime_error("Could not create Python dictionary");
        }
    }

    void add(std::string const& key, std::complex<double> value) override {
        auto* existing = PyDict_GetItemString(underlying, key.c_str());

        if (existing) {
            if (!PyComplex_Check(existing)) {
                throw std::runtime_error("Dictionary value is not a PyComplex");
            }
            
            auto* c = reinterpret_cast<PyComplexObject*>(existing);
            c->cval = _Py_c_sum(c->cval, Py_complex{.real = value.real(), .imag = value.imag()});
        } else {
            auto* pPyComplex = PyComplex_FromDoubles(value.real(), value.imag());
            if (!pPyComplex) {
                throw std::runtime_error("Cannot create PyComplex");
            }

            PyDict_SetItemString(underlying, key.c_str(), pPyComplex);
            Py_DECREF(pPyComplex);
        }
    }

    PyObject* getResult() {
        if (PyDict_Size(underlying) <= 0) {
            Py_DECREF(underlying);
            return NULL;
        }

        // Gives ownership to caller.
        return underlying;
    }

private:
    PyObject* underlying;
};

PyObject* getPythonState(qx::core::MixedStateBase const& mixedState) {
    auto m = PythonDictStringComplexMap();
    mixedState.getPureState(m);
    return m.getResult();
}

PyObject* getPythonDensityMatrix(qx::core::MixedStateBase const* mixedState) {
    auto const* mixedState64 = dynamic_cast<qx::core::MixedState<64> const*>(mixedState);
    if (!mixedState64) {
        auto* explanation = PyUnicode_FromString("Number of qubits is too large to output the density matrix; support for sparse density matrix is still to be added");
        if (!explanation) {
            throw std::runtime_error("Cannot create Python string to explain that the density matrix is too large");
        }

        return explanation;
    }

    npy_intp shape[2];
    shape[0] = 1 << mixedState64->getNumberOfQubits();
    shape[1] = 1 << mixedState64->getNumberOfQubits();

    PyObject* pNumpyArray = PyArray_ZEROS(2, shape, NPY_COMPLEX128, 0);
    if (!pNumpyArray) {
        throw std::runtime_error("Cannot create the Python numpy array to store the density matrix, maybe it is too large?");
    }

    auto const& data = mixedState64->getData();

    qx::EnsembleIndex currentEnsembleIndex = data.begin()->first.ensembleIndex;
    qx::utils::BasisVector currentMeasurementRegister = data.begin()->first.measurementRegister;

    for (auto it = data.begin(); it != data.end(); ++it) {
        if (it->first.measurementRegister != currentMeasurementRegister || it->first.ensembleIndex != currentEnsembleIndex) {
            currentEnsembleIndex = it->first.ensembleIndex;
            currentMeasurementRegister = it->first.measurementRegister;
        }
        
        auto* numpyArrayElement = reinterpret_cast<std::complex<double>*>(PyArray_GETPTR2(reinterpret_cast<PyArrayObject*>(pNumpyArray), it->first.basisVector.toUInt64(), it->first.basisVector.toUInt64()));
        *numpyArrayElement += std::norm(it->second);

        auto it2 = std::next(it);
        while (it2 != data.end() && it2->first.ensembleIndex == currentEnsembleIndex && it2->first.measurementRegister == currentMeasurementRegister) {
            numpyArrayElement = reinterpret_cast<std::complex<double>*>(PyArray_GETPTR2(reinterpret_cast<PyArrayObject*>(pNumpyArray), it2->first.basisVector.toUInt64(), it->first.basisVector.toUInt64()));
            *numpyArrayElement += it2->second * std::conj(it->second);
            numpyArrayElement = reinterpret_cast<std::complex<double>*>(PyArray_GETPTR2(reinterpret_cast<PyArrayObject*>(pNumpyArray), it->first.basisVector.toUInt64(), it2->first.basisVector.toUInt64()));
            *numpyArrayElement += it->second * std::conj(it2->second);
            ++it2;
        }
    }

    return pNumpyArray;
}

// This could be done in a numpy structured array instead of a Python dict but it's complicated.
// Also string keys have fixed sizes so they don't really need to be strings, could be inline.
class PythonDictStringDoubleMap : public qx::core::StringMap<double> {
public:
    PythonDictStringDoubleMap() {
        underlying = PyDict_New();
        if (!underlying) {
            throw std::runtime_error("Could not create PyDict");
        }
    }

    void add(std::string const& key, double value) override {
        auto* existing = PyDict_GetItemString(underlying, key.c_str());

        if (existing) {
            if (!PyFloat_Check(existing)) {
                throw std::runtime_error("Dictionary value is not a PyFloat");
            }

            reinterpret_cast<PyFloatObject*>(existing)->ob_fval += value;
        } else {
            auto* pPyFloat = PyFloat_FromDouble(value);
            if (!pPyFloat) {
                throw std::runtime_error("Cannot create PyFloat");
            }

            PyDict_SetItemString(underlying, key.c_str(), pPyFloat);
            Py_DECREF(pPyFloat);
        }
    }

    PyObject* getUnderlying() {
        // Gives ownership to caller.
        return underlying;
    }

private:
    PyObject* underlying;
};

PyObject* getPythonMeasurementRegisterStatistics(qx::core::MixedStateBase const& mixedState) {
    auto m = PythonDictStringDoubleMap();
    mixedState.getMeasurementRegisterStatistics(m);
    return m.getUnderlying();
}

PyObject* getPythonSimulationResult(std::unique_ptr<qx::core::MixedStateBase> quantumState) {
    if (!quantumState) {
        throw std::runtime_error("Simulator didn't return a valid quantum state");
    }

    auto* pQxelarator = PyImport_ImportModule("qxelarator");
    if (!pQxelarator) {
        throw std::runtime_error("Cannot import qxelarator Python module");
    }

    auto* pSimulationResult = PyObject_GetAttrString(pQxelarator, "SimulationResult");
    Py_DECREF(pQxelarator);
    if (!pSimulationResult) {
        throw std::runtime_error("Cannot import qxelarator.SimulationResult Python class");
    }

    auto* pSimulationResultInstance = PyObject_CallObject(pSimulationResult, NULL);
    Py_DECREF(pSimulationResult);
    if (!pSimulationResultInstance) {
        throw std::runtime_error("Cannot create Python qxelarator.SimulationResult instance");
    }

    auto* results = getPythonMeasurementRegisterStatistics(*quantumState);
    if (results) {
        // FIXME: it's unclear in the docs whether this steals the ref or if Py_DECREF needs to be called.
        PyObject_SetAttrString(pSimulationResultInstance, "results", results);
    }

    auto* state = getPythonState(*quantumState);
    if (state) { // FIXME
        PyObject_SetAttrString(pSimulationResultInstance, "state", state);
    } else {
        auto* densityMatrix = getPythonDensityMatrix(quantumState.get());
        if (densityMatrix) {
            PyObject_SetAttrString(pSimulationResultInstance, "densityMatrix", densityMatrix);
        }
    }

    return pSimulationResultInstance;
}

qx::Operations getOperationsFromPython(PyObject* input) {
    if (!PyDict_Check(input)) {
        throw std::runtime_error("Operations argument needs to be a Python dictionary");
    }

    if (PyDict_Size(input) == 0) {
        throw std::runtime_error("Python Operations dictionary argument is empty");
    }

    qx::Operations result;
    
    PyObject *key, *value;
    Py_ssize_t pos = 0;

    while (PyDict_Next(input, &pos, &key, &value)) {
        if (!PyUnicode_Check(key)) {
            throw std::runtime_error("Operations argument dictionary keys must be Python strings");
        }

        if (PyUnicode_GET_LENGTH(key) == 0) {
            throw std::runtime_error("Operations dictionary keys cannot be empty Python string");
        }

        auto name = qxelarator::getPythonString(key);

        if (!PyTuple_Check(value)) {
            throw std::runtime_error("Operations dictionary values must be Python tuples");
        }
        
        if (PyTuple_Size(value) != 2) {
            throw std::runtime_error("Operations dictionary value tuples must be Python tuples of size 2");
        }
        
        auto signature = qxelarator::getSignatureFromPython(PyTuple_GetItem(value, 0));

        auto* pGenerationFunction = PyTuple_GetItem(value, 1);
        if (!pGenerationFunction) {
            throw std::runtime_error("Cannot retrieve Python Kraus operators generation function");
        }

        Py_INCREF(pGenerationFunction);

        if (!PyFunction_Check(pGenerationFunction)) {
            throw std::runtime_error("Operations dictionary values must be Python callables/lambdas");
        }

        auto cppGenerationFunction = [name, pGenerationFunction](qx::Operations::StaticOperands const& operands) {
            auto* arguments = PyTuple_New(operands.size());
            if (!arguments) {
                throw std::runtime_error("Internal error");
            }

            for (std::uint64_t operandIndex = 0; operandIndex < operands.size(); ++operandIndex) {
                auto const& operand = operands[operandIndex];
                if (std::holds_alternative<double>(operand)) {
                    auto* pyDouble = PyFloat_FromDouble(std::get<double>(operand));
                    if (!pyDouble) {
                        throw std::runtime_error("Internal error");
                    }
                    PyTuple_SetItem(arguments, operandIndex, pyDouble); // Steals the ref to pyDouble.
                } else {
                    assert(std::holds_alternative<std::int64_t>(operand));
                    auto* pyLong = PyLong_FromLong(std::get<std::int64_t>(operand));
                    if (!pyLong) {
                        throw std::runtime_error("Internal error");
                    }
                    PyTuple_SetItem(arguments, operandIndex, pyLong); // Steals the ref to pyLong.
                }
            }

            auto* res = PyObject_CallObject(pGenerationFunction, arguments);
            if (!res) {
                auto argsString = qxelarator::getPythonString(PyObject_Str(arguments));
                Py_DECREF(arguments);
                PyErr_Clear();
                throw std::runtime_error(std::string("Call to Python function for operation ") + *name + " with arguments " + argsString.value_or("<unknown>") + " failed");
            }
            Py_DECREF(arguments);

            if (!PyList_Check(res)) {
                throw std::runtime_error("Kraus operators generation function needs to return a list of numpy.ndarray");
            }

            qx::Operations::KrausOperators cppKrausOperators;

            for (std::uint64_t krausOperatorIndex = 0; krausOperatorIndex < PyList_GET_SIZE(res); ++krausOperatorIndex) {
                auto* krausOperator = PyList_GetItem(res, krausOperatorIndex);
                if (!krausOperator) {
                    throw std::runtime_error("Cannot retrieve Kraus operator from Python function");
                }

                if (!PyArray_Check(krausOperator)) {
                    throw std::runtime_error("Kraus operators need to be instances of numpy.ndarray");
                }

                cppKrausOperators.push_back(qxelarator::numpyArrayToQXSquareMatrix(krausOperator));
            }


            return cppKrausOperators;
        };

        result.add(*name, *signature, cppGenerationFunction);
    }

    return result;
}

} // namespace qxelarator