"""Solver proof validation, byte identities and source consistency checks."""
from __future__ import annotations

from collections.abc import Mapping
from graphlib import CycleError, TopologicalSorter
import hashlib
import json
from pathlib import Path
import re

from ..cofactor.proof import Basis, Proof
from ..coefficients import Coef
from ..polynomial import Poly
from ..provenance import current_source_files, PREFIX


SCHEMA = "zipper-solver-run-v1"
ARTIFACT_NAMES = ("input.ms", "basis.txt", "proof.txt", "stdout.txt", "stderr.txt")


class ArtifactError(ValueError):
    """A recorded computation does not match its declared inputs or proof."""


def validate_source_snapshot(producer):
    """Refuse publication if the recorded package source bytes changed."""
    expected = {name: identity for name, identity in producer["source_files"].items()
                if name.startswith(PREFIX)}
    try:
        current = current_source_files(expected)
    except (OSError, ValueError) as error:
        raise ArtifactError("zipper-algebra source files changed or became unreadable during the computation") from error
    if current != expected:
        raise ArtifactError("zipper-algebra source files changed during the computation")


def prepare_output_root(root: Path):
    """Keep generated output local while preserving an existing ignore policy."""
    if (root / ".git").exists():
        raise ArtifactError("output directory must not be a Git checkout root")
    existed = root.is_dir()
    root.mkdir(parents=True, exist_ok=True)
    if existed:
        absolute = root.resolve()
        checkout = next((parent for parent in absolute.parents
                         if (parent / ".git").is_file()
                         or (parent / ".git" / "HEAD").is_file()), None)
        if checkout is not None and absolute.relative_to(checkout).parts[0] not in (
                "runs", ".research-data", ".cache"):
            return
    marker = root / ".gitignore"
    try:
        with marker.open("x") as stream:
            stream.write("*\n!.gitignore\n")
    except FileExistsError:
        if not marker.is_file():
            raise




def canonical_json(value) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()


def digest_bytes(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def file_identity(path: Path) -> dict:
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as stream:
        while chunk := stream.read(131072):
            size += len(chunk)
            digest.update(chunk)
    return {"bytes": size, "sha256": digest.hexdigest()}


def _check_syntax(text: str, names, *, proof: bool):
    """Consume complete solver terms before using the historical tolerant parser."""
    variable = "(?:" + "|".join(re.escape(name) for name in names) + ")"
    scalar = r"\d+(?:/\d+)?"
    word = variable + r"(?:\*" + variable + ")*"
    if proof:
        term = r"[+-]?(?:" + scalar + r"\*)?(?:" + variable + r"\*)*\[i?\d+\](?:\*" + variable + ")*"
    else:
        term = r"[+-]?(?:" + scalar + r"(?:\*" + word + ")?|" + word + ")"
    pattern = re.compile(term)
    for raw in text.splitlines():
        line = raw.strip()
        if line.endswith(","):
            line = line[:-1]
        tokens = line.split()
        pending = False
        for index, token in enumerate(tokens):
            if token in ("+", "-"):
                if pending:
                    raise ValueError("repeated standalone sign in solver output")
                pending = True
            else:
                if pending and token.startswith(("+", "-")):
                    raise ValueError("repeated sign in solver output")
                if index and not pending and not token.startswith(("+", "-")):
                    raise ValueError("missing sign between solver output terms")
                if pattern.fullmatch(token) is None:
                    raise ValueError("malformed solver output term")
                pending = False
        if pending:
            raise ValueError("incomplete solver output term")


def validate_proofs(source, inputs, names) -> int:
    """Check every output identity against the actual ordered input relations.

    Local equalities on the acyclic derivation graph establish that the output
    polynomials belong to the input ideal, without expanding shared proofs.
    It does not establish that they form a complete Gröbner basis.
    """
    try:
        def read(name):
            if isinstance(source, Mapping):
                data = source[name]
                if not isinstance(data, bytes):
                    raise ValueError("solver artifact contents must be bytes")
                return data.decode("utf-8")
            return (Path(source) / name).read_text()
        basis_text, proof_text = read("basis.txt"), read("proof.txt")
        _check_syntax(basis_text, names, proof=False)
        _check_syntax(proof_text, names, proof=True)
        basis = Basis.from_text(basis_text, 0, names)
        proof = Proof.from_text(proof_text, names)
        if not len(basis) or len(basis) != len(proof):
            raise ValueError("basis and proof must contain the same positive number of entries")
        for terms in proof.terms:
            for _coefficient, _left, kind, index, _right in terms:
                bound = len(inputs) if kind == "input" else len(basis)
                if not 0 <= index < bound:
                    raise ValueError("a raw proof reference does not name an existing relation")
        try:
            TopologicalSorter({index: proof.refs(index) for index in range(len(basis))}).prepare()
        except CycleError as error:
            raise ValueError("cyclic proof references") from error
        for index in range(len(basis)):
            reconstructed = Poly(0)
            for coefficient, left, kind, reference, right in proof.terms[index]:
                polynomial = inputs[reference] if kind == "input" else basis.polys[reference]
                scalar = Coef(coefficient, 0)
                for word, value in polynomial.terms.items():
                    reconstructed.add_term(left + word + right, scalar * value)
            if (reconstructed - basis.polys[index]).terms:
                raise ValueError(f"basis entry {index} fails its input-ideal identity")
    except (OSError, ValueError, KeyError, IndexError, ZeroDivisionError) as error:
        raise ArtifactError(f"invalid solver proof: {error}") from error
    return len(basis)
