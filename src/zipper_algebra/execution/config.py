"""Configuration for one solver-backed computation, without workspace defaults."""
from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path


DEFAULT_FLAGS = ("-e", "0", "-v", "2", "-r", "0")


@dataclass(frozen=True)
class RunConfig:
    """Paths and limits passed explicitly to the runner.

    ``command_prefix`` permits an interpreter followed by executable source files;
    the contents of those files are part of the solver identity. It is used by
    the deterministic process-boundary tests, not by the normal f4ncgb command.
    """

    binary: Path
    output_root: Path
    cache_root: Path | None = None
    timeout_seconds: float = 1800
    max_iterations: int = 100000
    max_degree: int | None = None
    flags: tuple[str, ...] = DEFAULT_FLAGS
    solver_revision: str = "unrecorded"
    monomial_order: str = "degree-lex"
    command_prefix: tuple[Path, ...] = ()

    def __post_init__(self):
        for name in ("binary", "output_root", "cache_root"):
            value = getattr(self, name)
            if value is None:
                if name != "cache_root":
                    raise ValueError(f"{name} must be an explicit path")
            else:
                object.__setattr__(self, name, Path(value).resolve())
        if not self.binary.is_file():
            raise ValueError("binary must name an existing executable file")
        if (isinstance(self.timeout_seconds, bool)
                or not isinstance(self.timeout_seconds, (int, float))
                or not math.isfinite(self.timeout_seconds) or self.timeout_seconds <= 0):
            raise ValueError("timeout_seconds must be a finite positive number")
        for name in ("max_iterations", "max_degree"):
            value = getattr(self, name)
            if value is not None and (type(value) is not int or value < 1):
                raise ValueError(f"{name} must be a positive integer")
        if self.max_iterations is None:
            raise ValueError("max_iterations must be a positive integer")
        if (not isinstance(self.flags, tuple)
                or any(not isinstance(arg, str) or not arg for arg in self.flags)):
            raise ValueError("flags must be a tuple of nonempty strings")
        reserved = {"-i", "--input", "-o", "--output", "-p", "--proof",
                    "-m", "--maxiter", "-d", "--maxdeg"}
        if any(arg.split("=", 1)[0] in reserved
               or any(arg.startswith(short) for short in ("-i", "-o", "-p", "-m", "-d"))
               for arg in self.flags):
            raise ValueError("input, output, proof and solver limits are owned by RunConfig")
        if not isinstance(self.solver_revision, str) or not self.solver_revision:
            raise ValueError("solver_revision must be a nonempty description")
        if self.monomial_order != "degree-lex":
            raise ValueError("this runner supports only the degree-lex monomial order")
        if not isinstance(self.command_prefix, tuple):
            raise ValueError("command_prefix must be a tuple of executable source paths")
        prefix = tuple(Path(path).resolve() for path in self.command_prefix)
        if any(not path.is_file() for path in prefix):
            raise ValueError("each command_prefix path must be an existing file")
        object.__setattr__(self, "command_prefix", prefix)

    def arguments(self) -> tuple[str, ...]:
        """Solver flags including explicit iteration and optional degree limits."""
        degree = () if self.max_degree is None else ("-d", str(self.max_degree))
        return (*self.flags, "-m", str(self.max_iterations), *degree)
