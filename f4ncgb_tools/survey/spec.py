"""Algebra specification for a survey point."""
from __future__ import annotations
from dataclasses import dataclass

# untwisting name -> (zipper_in_ms --untwist value, --framing-power)
UNTWIST = {
    "linear": ("bigelow", 1),    # p_2 = sigma - t      ; p* = L + t
    "jordan": ("jordan", 1),     # p_2 = (sigma - t)^2  ; p* = (L + t)^2
    "framing": ("bigelow", 2),   # p_2 = sigma^2 - t    ; p* = L^2 - t
}


@dataclass(frozen=True)
class Spec:
    untwist: str                 # key of UNTWIST
    n: int                       # number of strands
    max_untwist_level: int = 2   # 2 = two-strand untwisting only (the survey default)

    @property
    def emit_untwist(self) -> str:
        return UNTWIST[self.untwist][0]

    @property
    def framing_power(self) -> int:
        return UNTWIST[self.untwist][1]

    @property
    def name(self) -> str:
        return f"{self.untwist}-n{self.n}-mul{self.max_untwist_level}"
