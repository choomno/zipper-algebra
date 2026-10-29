"""Prepare ignored directories for optional experiment outputs."""

from pathlib import Path
import hashlib
import os
import tempfile


def prepare_output_dir(path: Path | str) -> Path:
    """Create an output directory and preserve any existing ignore rules."""
    directory = Path(path)
    if (directory / ".git").exists():
        raise ValueError("output directory must not be a Git checkout root")
    existed = directory.is_dir()
    directory.mkdir(parents=True, exist_ok=True)
    if existed:
        absolute = directory.resolve()
        checkout = next((parent for parent in absolute.parents
                         if (parent / ".git").is_file()
                         or (parent / ".git" / "HEAD").is_file()), None)
        if checkout is not None and absolute.relative_to(checkout).parts[0] not in (
                "runs", ".research-data", ".cache"):
            return directory
    try:
        with (directory / ".gitignore").open("x", encoding="utf-8") as marker:
            marker.write("*\n!.gitignore\n")
    except FileExistsError:
        pass
    return directory


def capture_experiment_provenance(source, *helpers):
    """Snapshot the package, recipe family and explicitly shared source helpers."""
    from zipper_algebra.provenance import capture_provenance
    source = Path(source)
    sources = [*source.parent.glob("*.py"), Path(__file__),
               Path(__file__).with_name("__init__.py"), *helpers]
    return capture_provenance(extra_sources=sources)


def input_identity(path):
    """Identify an input before a transformation can replace its bytes."""
    from zipper_algebra.provenance import file_identity
    path = Path(path)
    return {"name": path.name, **file_identity(path)}


def read_input(path):
    """Return input bytes and their identity from the same read."""
    path = Path(path)
    content = path.read_bytes()
    return content, {"name": path.name, "bytes": len(content),
                     "sha256": hashlib.sha256(content).hexdigest()}


def save_pickle_output(path, payload, *, producer, inputs=()):
    """Keep the legacy pickle payload unchanged and write its provenance beside it."""
    import pickle
    from zipper_algebra.provenance import write_output_provenance
    path = Path(path)
    stream = tempfile.NamedTemporaryFile(mode="wb", dir=path.parent,
                                         prefix="." + path.name + ".", delete=False)
    temporary = Path(stream.name)
    try:
        with stream:
            pickle.dump(payload, stream)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)
    write_output_provenance(path, producer=producer, inputs=inputs)


def tag_solver_input(path, *, producer):
    """Tag an emitted solver input and retain its identity before execution."""
    from zipper_algebra.provenance import write_output_provenance
    inputs = [input_identity(path)]
    write_output_provenance(path, producer=producer)
    return inputs


def tag_solver_basis(path, *, producer, inputs):
    """Tag an existing solver output without changing its polynomial format."""
    from zipper_algebra.provenance import write_output_provenance
    if Path(path).is_file():
        write_output_provenance(path, producer=producer, inputs=inputs)
