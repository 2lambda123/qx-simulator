class SimulationResult:
    def __init__(self):
        self.results = None
        self.state = None
        self.densityMatrix = None

    def __repr__(self):
        if self.densityMatrix is not None:
            return f"""Density matrix:
{self.densityMatrix}

Measurement register probabilities:
{self.results}"""
        else:
            return f"""State:
{self.state}

Measurement register probabilities:
{self.results}"""

class SimulationError:
    def __init__(self, message):
        self.message = message

    def __repr__(self):
        return f"Quantum simulation error: {self.message}"

from .qxelarator import *