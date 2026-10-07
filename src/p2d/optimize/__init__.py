"""Constrained codon optimisation (v0.2).

``optimize_coding`` is a DFA-constrained DP; see ``core.py``.
"""

from .automaton import Automaton, Match, expand_pattern, reverse_complement
from .core import CodonStrategy, OptimizationRequest, OptimizationResult, optimize_coding

__all__ = [
    "Automaton",
    "CodonStrategy",
    "Match",
    "OptimizationRequest",
    "OptimizationResult",
    "expand_pattern",
    "optimize_coding",
    "reverse_complement",
]
