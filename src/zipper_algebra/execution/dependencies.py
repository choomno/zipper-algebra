"""Identify the installed native reconstruction code used by a computation."""
from __future__ import annotations

import hashlib
import importlib
import os
from pathlib import Path
import stat


def _mapped_files():
    """Read Linux's loaded file identities, retaining separate inode generations."""
    loaded: dict[str, set[tuple[int, int]]] = {}
    for line in Path("/proc/self/maps").read_text().splitlines():
        fields = line.split(maxsplit=5)
        if len(fields) != 6 or not fields[5].startswith("/"):
            continue
        path = fields[5]
        major, minor = (int(piece, 16) for piece in fields[3].split(":"))
        loaded.setdefault(path, set()).add((os.makedev(major, minor), int(fields[4])))
    return loaded


def _file_record(path, mapped=None):
    """Hash one stable file and, for native objects, verify its loaded inode."""
    path = Path(path)
    with path.open("rb") as stream:
        before = os.fstat(stream.fileno())
        if not stat.S_ISREG(before.st_mode):
            raise RuntimeError("reconstruction dependency is not a regular file")
        # Overlay mounts can report a different st_dev through stat than the
        # underlying mapping. Require one mapped generation and matching inode;
        # the on-disk device is still included in the hash-time stability check.
        if mapped is not None and (len(mapped) != 1 or next(iter(mapped))[1] != before.st_ino):
            raise RuntimeError("reconstruction dependency differs from the loaded native object")
        digest = hashlib.file_digest(stream, "sha256").hexdigest()
        after = os.fstat(stream.fileno())
    final = path.stat()
    def identity(value):
        return value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns
    if identity(before) != identity(after) or identity(before) != identity(final):
        raise RuntimeError("reconstruction dependency changed while recording its identity")
    return {"name": path.name, "bytes": before.st_size, "sha256": digest}


def reconstruction_provenance() -> dict:
    """Record exact binding and loaded-library bytes without installation paths.

    Linux's process maps identify the library actually loaded by the CFFI
    extension, rather than a library guessed from an installation prefix.
    Native inode checks reject a replaced or deleted on-disk library. Capture
    again before publication and compare records to detect intervening changes.
    Version and fork metadata are declarations by the installed Python binding;
    no source revision is inferred from the required installation pin. A file
    edited in place before this first capture cannot be compared with the bytes
    originally loaded by Python or the dynamic linker; this is an installed-file
    identity record, not a memory-image attestation.
    """
    try:
        binding = importlib.import_module("fflow")
        extension = importlib.import_module("_cffi_fflow")
    except ImportError as error:
        raise ImportError("discovery requires the FiniteFlow Python bindings; "
                          "see docs/RECONSTRUCTION.md") from error
    version = getattr(binding, "__version__", None)
    fork = getattr(binding, "__fork__", None)
    api = getattr(binding, "__fork_api_version__", None)
    if (not isinstance(version, str) or not version or not isinstance(fork, str) or not fork
            or type(api) is not int or api < 1):
        raise RuntimeError("the installed FiniteFlow binding does not identify the maintained fork API")
    binding_path = getattr(binding, "__file__", None)
    extension_path = getattr(extension, "__file__", None)
    if (not isinstance(binding_path, str) or not binding_path
            or not isinstance(extension_path, str) or not extension_path):
        raise RuntimeError("reconstruction dependencies do not identify their installed files")
    try:
        loaded = _mapped_files()
        libraries = [path for path in loaded
                     if Path(path).name == "libfflow.so" or Path(path).name.startswith("libfflow.so.")
                     or Path(path).name == "libfflow.so (deleted)"]
        if len(libraries) != 1 or libraries[0].endswith(" (deleted)"):
            raise RuntimeError("expected exactly one available loaded FiniteFlow library")
        # Extension __file__ can be a symbolic path. Process maps normally use
        # its resolved filename, so compare using that same absolute filename.
        extension_path = str(Path(extension_path).resolve())
        if extension_path not in loaded:
            raise RuntimeError("the FiniteFlow CFFI extension is not present in the loaded native objects")
        library_path = libraries[0]
        files = {"binding": _file_record(binding_path),
                 "extension": _file_record(extension_path, loaded[extension_path]),
                 "library": _file_record(library_path, loaded[library_path])}
    except OSError as error:
        raise RuntimeError("could not read the installed reconstruction dependencies") from error
    return {"schema": "zipper-reconstruction-dependencies-v1", "implementation": "FiniteFlow",
            "version": version, "fork": fork, "fork_api_version": api, "files": files}
