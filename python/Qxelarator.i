%define DOCSTRING
"`qxelarator` Python interface to QX simulator as an accelerator."
%enddef

%module(docstring=DOCSTRING) qxelarator

%include "std_string.i"

%typemap(in) std::optional<std::uint_fast64_t> {
    if($input == Py_None) {
        $1 = std::nullopt;
    } else if (PyLong_Check($input)) {
        $1 = PyLong_AsLong($input); // FIXME
    } else {
        return qxelarator::getPythonSimulationError("Seed argument must be a Python integer");
    }
}

%{
#include "Qxelarator.hpp"
#include "Helpers.hpp"
%}

%typemap(in) qx::Operations {
try {
    import_array();
    $1 = qxelarator::getOperationsFromPython($input);
} catch (std::exception const& e) {
    return qxelarator::getPythonSimulationError(e.what());
}
}

%exception {
    try {
        $action
    } catch (std::exception const& e) {
        return qxelarator::getPythonSimulationError(e.what());
    }
}

// Map the output of execute_string/execute_file to a simple Python class for user-friendliness.
%typemap(out) std::unique_ptr<qx::core::MixedStateBase> {
try {
    import_array();
    $result = qxelarator::getPythonSimulationResult(std::move($1));
} catch (std::exception const& e) {
    return qxelarator::getPythonSimulationError(e.what());
}
}

// Include the header file with above prototypes
%include "Qxelarator.hpp"
