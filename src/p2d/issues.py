"""Diagnostics shared by the frame engine and the optimiser.

Kept out of ``frame.py`` so ``p2d.optimize`` can report issues without importing
the module that calls it.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class Severity(str, Enum):
    ERROR = "error"
    WARNING = "warning"
    INFO = "info"


@dataclass(frozen=True)
class Issue:
    severity: Severity
    code: str
    message: str

    def __str__(self) -> str:
        return f"[{self.severity.value}] {self.code}: {self.message}"
