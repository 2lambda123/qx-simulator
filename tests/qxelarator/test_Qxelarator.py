import unittest
import os
import qxelarator
import numpy as np
import math

class QxelaratorTest(unittest.TestCase):
    def test_execute_file(self):
        cQasmFileName = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'bell_pair.qasm')
        simulationResult = qxelarator.execute_file(cQasmFileName)

        self.assertIsInstance(simulationResult, qxelarator.SimulationResult)
        self.assertEqual(simulationResult.results, {"00": 1.0})
        self.assertEqual(simulationResult.state, {"11": complex(1., 0.)})

    def test_execute_string_fails_returns_error(self):
        cQasmString = \
"""
version 1.0

qubits 1

.myCircuit
    x q[0],
    measure q[0], b[0]
"""
        simulationError = qxelarator.execute_string(cQasmString)

        self.assertIsInstance(simulationError, qxelarator.SimulationError)
        self.assertEqual(simulationError.message, "<unknown>:7:12: syntax error, unexpected NEWLINE")

    def test_execute_string(self):
        cQasmString = \
"""
version 1.0

qubits 1

.myCircuit
    x q[0]
    measure q[0], b[0]
"""
        simulationResult = qxelarator.execute_string(cQasmString)

        self.assertIsInstance(simulationResult, qxelarator.SimulationResult)
        self.assertEqual(simulationResult.results, {"1": 1.0})
        self.assertEqual(simulationResult.state, {"1": complex(1., 0.)})
        self.assertIsNone(simulationResult.densityMatrix)
    
    def test_density_matrix(self):
        cQasmString = \
"""
version 1.0
qubits 1

h q[0]
measure q[0], b[0]
"""
        simulationResult = qxelarator.execute_string(cQasmString)

        self.assertIsInstance(simulationResult.densityMatrix, np.ndarray)
        self.assertTrue(np.allclose(simulationResult.densityMatrix, np.array([[0.5, 0], [0, 0.5]])))

        cQasmString = \
"""
version 1.0
qubits 1

h q[0]
"""
        simulationResult = qxelarator.execute_string(cQasmString)
        self.assertIsNone(simulationResult.densityMatrix)

        cQasmString = \
"""
version 1.0
qubits 100

h q[0]
measure q[0], b[0]
"""
        simulationResult = qxelarator.execute_string(cQasmString)
        self.assertIsNone(simulationResult.state)
        self.assertEqual(simulationResult.densityMatrix, "Number of qubits is too large to output the density matrix; support for sparse density matrix is still to be added")
    
    def test_seed_deterministic(self):
        cQasmString = \
"""
version 1.0
qubits 1

.myCircuit
    h q[0]
    measure q[0], b[0]
"""
        simulationResult = qxelarator.execute_string(cQasmString, seed = 123)
        self.assertEqual(list(simulationResult.results.keys()), ["1", "0"])
        self.assertAlmostEqual(simulationResult.results["0"], 0.5)
        self.assertAlmostEqual(simulationResult.results["1"], 0.5)

    def test_execute_string_custom_operations(self):
        cQasmString = \
"""
version 1.0
qubits 1

.myCircuit
    mycustom q[0]
"""

        operations = {"mycustom": ("q", [np.array([[0, 1], [1, 0]])])}
        output = qxelarator.execute_string(cQasmString, operations = operations)
        self.assertIsInstance(output, qxelarator.SimulationError)
        self.assertEqual(output.message, "Operations dictionary values must be Python callables/lambdas")

        operations = {"mycustom": ("dq", lambda : [np.array([[0, 1], [1, 0]])])}
        output = qxelarator.execute_string(cQasmString, operations = operations)
        self.assertIsInstance(output, qxelarator.SimulationError)
        self.assertEqual(output.message, "Requested operation mycustom with signature (Qubit), but the registered operation with that name has signature (Double, Qubit)")

        operations = {"mycustom": ("q", lambda : [np.array([[0, 1], [1, 1]])])}
        output = qxelarator.execute_string(cQasmString, operations = operations)
        self.assertIsInstance(output, qxelarator.SimulationError)
        self.assertEqual(output.message, "Kraus operators for operation mycustom are not non-trace-increasing")

        operations = {"mycustom": ("q", lambda : [np.array([[0, 1]])])}
        output = qxelarator.execute_string(cQasmString, operations = operations)
        self.assertIsInstance(output, qxelarator.SimulationError)
        self.assertEqual(output.message, "Not a square matrix")

        operations = {"mycustom": ("q", lambda j : [np.array([[0, 1], [1, 0]])])}
        output = qxelarator.execute_string(cQasmString, operations = operations)
        self.assertIsInstance(output, qxelarator.SimulationError)
        self.assertEqual(output.message, "Call to Python function for operation mycustom with arguments () failed")

        operations = {"mycustom": ("qd", lambda : [np.array([[0, 1], [1, 0]])])}
        output = qxelarator.execute_string(cQasmString, operations = operations)
        self.assertIsInstance(output, qxelarator.SimulationError)
        self.assertEqual(output.message, "Requested operation mycustom with signature (Qubit), but the registered operation with that name has signature (Qubit, Double)")

        operations = {"mycustom": ("q", lambda : [np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])])}
        output = qxelarator.execute_string(cQasmString, operations = operations)
        self.assertIsInstance(output, qxelarator.SimulationError)
        self.assertEqual(output.message, "Kraus operators for operation mycustom have size 4, which do not match the number of operands (1)")

        operations = {"mycustom": ("q", lambda : [np.array([[0, 1], [1, 0]])])}
        output = qxelarator.execute_string(cQasmString, operations = operations)
        self.assertIsInstance(output, qxelarator.SimulationResult)
        self.assertEqual(output.state, {"1": 1.0 + 0.0j})

        operations = {"mycustom": ("q", lambda : [np.array([[1, 0], [0, 0]]), np.array([[0, 0], [0, 1]])])}
        output = qxelarator.execute_string(cQasmString, operations = operations)
        self.assertIsInstance(output, qxelarator.SimulationResult)
        self.assertEqual(output.state, {"0": 1.0 + 0.0j})
        

    
if __name__ == '__main__':
    unittest.main()
