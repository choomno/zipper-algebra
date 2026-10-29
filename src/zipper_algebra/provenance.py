"""Deterministic source snapshots and output metadata, without source contents.

The snapshot identifies files on disk when captured. It is neither an authorship
signature nor a copy of uncommitted source, and does not introspect code already
loaded into a long-lived Python interpreter. Git dirty state describes Git's
tracked/declared-source changes. A positive omitted_source_count identifies a
partial source snapshot because ignored candidates were excluded; their names
and contents are not recorded.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import subprocess
import tempfile
import tomllib


SCHEMA = "zipper-source-provenance-v1"
BUILD_SCHEMA = "zipper-build-provenance-v1"
OUTPUT_SCHEMA = "zipper-output-provenance-v1"
PREFIX = "src/zipper_algebra/"
BUILD_SOURCES = ("MANIFEST.in", "pyproject.toml", "setup.py")


def _canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def file_identity(path):
    """Hash the exact bytes of one file."""
    content = Path(path).read_bytes()
    return {"bytes": len(content), "sha256": hashlib.sha256(content).hexdigest()}


def _digest(value):
    return hashlib.sha256(_canonical(value)).hexdigest()


def _relative(name):
    if not isinstance(name, str):
        raise ValueError("provenance paths must be relative public file names")
    path = PurePosixPath(name)
    if (not name or path.is_absolute() or ".." in path.parts
            or "\\" in name or ":" in name or any(ord(c) < 32 for c in name)):
        raise ValueError("provenance paths must be relative public file names")
    return name


def _git_relative(name):
    """Preserve legal Git names as escaped JSON, including colons and newlines."""
    if (not isinstance(name, str) or not name or PurePosixPath(name).is_absolute()
            or ".." in PurePosixPath(name).parts or "\0" in name):
        raise ValueError("provenance Git paths must be relative file names")
    return name


def _plain_name(name):
    _relative(name)
    return PurePosixPath(name).name == name and name not in (".", "..")


def _hash(value):
    return isinstance(value, str) and len(value) == 64 and all(c in "0123456789abcdef" for c in value)


def _identity(value):
    return (isinstance(value, dict) and set(value) == {"bytes", "sha256"}
            and type(value["bytes"]) is int and value["bytes"] >= 0 and _hash(value["sha256"]))


def _validate_git(value):
    if not isinstance(value, dict) or set(value) != {"commit", "dirty", "changes"}:
        raise ValueError("invalid provenance Git state")
    commit = value["commit"]
    if commit is not None and (not isinstance(commit, str) or len(commit) != 40
                               or any(c not in "0123456789abcdef" for c in commit)):
        raise ValueError("invalid provenance Git commit")
    if value["dirty"] is not None and type(value["dirty"]) is not bool:
        raise ValueError("invalid provenance dirty state")
    if not isinstance(value["changes"], list):
        raise ValueError("invalid provenance change list")
    for change in value["changes"]:
        if (not isinstance(change, dict)
                or not {"path", "index", "worktree", "bytes", "sha256"} <= set(change)
                or set(change) - {"path", "index", "worktree", "bytes", "sha256", "previous_path"}
                or any(not isinstance(change[k], str) or len(change[k]) != 1
                       for k in ("index", "worktree"))):
            raise ValueError("invalid provenance change")
        _git_relative(change["path"])
        if "previous_path" in change:
            _git_relative(change["previous_path"])
        if (change["sha256"], change["bytes"]) != (None, None) and not _identity(
                {"sha256": change["sha256"], "bytes": change["bytes"]}):
            raise ValueError("invalid provenance changed-file hash")


def validate_provenance(value):
    """Validate metadata shape and internal hashes, without authenticating it."""
    fields = {"schema", "origin", "git", "source_files", "source_sha256",
              "changes_since_build", "extra_source_git", "omitted_source_count"}
    if (not isinstance(value, dict) or set(value) != fields or value["schema"] != SCHEMA
            or value["origin"] not in ("checkout", "build", "unversioned")):
        raise ValueError("invalid source provenance schema")
    _validate_git(value["git"])
    if type(value["omitted_source_count"]) is not int or value["omitted_source_count"] < 0:
        raise ValueError("omitted source count must be a nonnegative integer")
    if not isinstance(value["source_files"], dict) or not value["source_files"]:
        raise ValueError("source provenance needs source files")
    for name, identity in value["source_files"].items():
        _relative(name)
        computation = name.endswith(".py") and (name.startswith(PREFIX) or name.startswith("experiments/"))
        if not (computation or name in BUILD_SOURCES) or not _identity(identity):
            raise ValueError("invalid provenance source identity")
    if value["source_sha256"] != _digest(value["source_files"]):
        raise ValueError("source provenance digest mismatch")
    if not isinstance(value["changes_since_build"], list):
        raise ValueError("invalid post-build change list")
    for change in value["changes_since_build"]:
        if (not isinstance(change, dict) or set(change) != {"path", "status"}
                or change["status"] not in ("added", "modified", "deleted")):
            raise ValueError("invalid post-build change")
        _relative(change["path"])
    if not isinstance(value["extra_source_git"], list):
        raise ValueError("invalid additional-source Git state")
    for entry in value["extra_source_git"]:
        if not isinstance(entry, dict) or set(entry) != {"paths", "git"} or not isinstance(entry["paths"], list):
            raise ValueError("invalid additional-source origin")
        for name in entry["paths"]:
            _relative(name)
            if name not in value["source_files"] and {"path": name, "status": "deleted"} not in value["changes_since_build"]:
                raise ValueError("additional-source origin names an absent source")
        _validate_git(entry["git"])
    return value


def _project(root):
    marker = root / "pyproject.toml"
    if not marker.is_file():
        return False
    return tomllib.loads(marker.read_text(encoding="utf-8")).get("project", {}).get("name") == "zipper-algebra"


def _git(root, *arguments, input=None):
    return subprocess.check_output(["git", "-C", str(root), *arguments], input=input,
                                   stderr=subprocess.DEVNULL,
                                   env=dict(os.environ, GIT_OPTIONAL_LOCKS="0"))


def _ignored_sources(root, names):
    if root is None or not (root / ".git").exists() or not names:
        return set()
    try:
        data = _git(root, "check-ignore", "-z", "--stdin",
                    input=b"\0".join(os.fsencode(name) for name in names) + b"\0")
    except (OSError, subprocess.CalledProcessError):
        # Git may be absent; the separately recorded Git state remains unknown.
        # Exit 1 also means the query succeeded with no ignored candidates.
        return set()
    return {os.fsdecode(name) for name in data.split(b"\0") if name}


def _unknown_git():
    return {"commit": None, "dirty": None, "changes": []}


def _git_state(root, sources):
    if root is None or not (root / ".git").exists():
        return _unknown_git()
    try:
        top = Path(os.fsdecode(_git(root, "rev-parse", "--show-toplevel")).strip()).resolve()
        if top != root:
            return _unknown_git()
        commit = _git(root, "rev-parse", "--verify", "HEAD").decode().strip()
        status = _git(root, "status", "--porcelain=v1", "-z", "--untracked-files=no").split(b"\0")
        tracked = {os.fsdecode(name) for name in _git(root, "ls-files", "-z").split(b"\0") if name}
    except (OSError, subprocess.CalledProcessError):
        return _unknown_git()
    changes = []
    index = 0
    while index < len(status) and status[index]:
        item = os.fsdecode(status[index]); index += 1
        name = _git_relative(item[3:])
        change = {"path": name, "index": item[0], "worktree": item[1]}
        if "R" in item[:2] or "C" in item[:2]:
            change["previous_path"] = _git_relative(os.fsdecode(status[index])); index += 1
        current = root / name
        identity = file_identity(current) if current.is_file() and current.resolve().is_relative_to(root) else {"bytes": None, "sha256": None}
        changes.append({**change, **identity})
    for name in sorted(set(sources) - tracked):
        changes.append({"path": name, "index": "?", "worktree": "?", **file_identity(root / name)})
    return {"commit": commit, "dirty": bool(changes), "changes": sorted(changes, key=lambda item: item["path"])}


def _extra_sources(paths, project, baseline):
    roots, omitted = {}, set()
    for item in paths:
        path = Path(item).resolve()
        if path.suffix != ".py":
            raise ValueError("additional provenance sources must be Python files")
        root = project
        if root is None:
            # Only inspect ancestors of the explicitly declared recipe, not CWD.
            root = next((parent for parent in path.parents
                         if (parent / "experiments").is_dir() and _project(parent)), None)
        if root is None or not path.is_relative_to(root):
            raise ValueError("additional source is outside the declared project")
        name = _relative(path.relative_to(root).as_posix())
        if not (name.startswith("experiments/") or name.startswith(PREFIX)):
            raise ValueError("additional source is outside public computation sources")
        if name in _ignored_sources(root, [name]):
            omitted.add((root, name))
            continue
        if not path.is_file() and (baseline is None or name not in baseline["source_files"]):
            raise ValueError("additional provenance sources must exist or have a build baseline")
        roots.setdefault(root, {})[name] = file_identity(path) if path.is_file() else None
    if len(roots) > 1:
        raise ValueError("additional sources must belong to one project")
    return roots, omitted


def current_source_files(names=None):
    """Hash package source files without Git; optionally rehash only declared names."""
    package = Path(__file__).resolve().parent
    if names is None:
        paths = sorted(package.rglob("*.py"))
    else:
        paths = []
        for name in names:
            _relative(name)
            if not name.startswith(PREFIX) or not name.endswith(".py"):
                raise ValueError("source rehash names must identify package Python files")
            paths.append(package / name[len(PREFIX):])
    sources = {}
    for path in paths:
        if not path.resolve().is_relative_to(package):
            raise ValueError("package source path leaves its source tree")
        if path.is_file():
            sources[PREFIX + path.relative_to(package).as_posix()] = file_identity(path)
    return sources


def _build_origin(package, project):
    build = package / "_build_provenance.json"
    if (project is not None and (project / ".git").exists()) or not build.is_file():
        return None
    metadata = json.loads(build.read_text(encoding="utf-8"))
    if not isinstance(metadata, dict) or set(metadata) != {"schema", "origin"} or metadata["schema"] != BUILD_SCHEMA:
        raise ValueError("invalid build provenance schema")
    return validate_provenance(metadata["origin"])


def capture_provenance(*, extra_sources=()):
    """Capture source bytes and revision at operation entry, without private notes."""
    package = Path(__file__).resolve().parent
    candidate = package.parent.parent
    project = candidate if package.parent.name == "src" and _project(candidate) else None
    sources = current_source_files()
    ignored = _ignored_sources(project, sources)
    omitted = {(project, name) for name in ignored}
    for name in ignored:
        sources.pop(name, None)
    if not sources:
        raise ValueError("source provenance needs source files")
    if project is not None:
        for name in BUILD_SOURCES:
            path = project / name
            if path.is_file():
                if not path.resolve().is_relative_to(project):
                    raise ValueError("build source path leaves its source tree")
                sources[name] = file_identity(path)
    baseline = _build_origin(package, project)
    extras, omitted_extras = _extra_sources(extra_sources, project, baseline)
    omitted.update(omitted_extras)
    requested = {name for entries in extras.values() for name in entries}
    for entries in extras.values():
        sources.update({name: identity for name, identity in entries.items() if identity is not None})
    state = _git_state(project, sources)
    origin, changes = "checkout" if state["commit"] is not None else "unversioned", []
    if baseline is not None:
        state, origin = baseline["git"], "build"
        before = baseline["source_files"]
        compared = {name for name in before.keys() | sources.keys() if name.startswith(PREFIX)} | requested
        if project is not None:
            compared.update(BUILD_SOURCES)
        changes = [{"path": name, "status": "added" if name not in before else "deleted" if name not in sources else "modified"}
                   for name in sorted(compared) if before.get(name) != sources.get(name)]
    extra_git = [{"paths": sorted(entries), "git": _git_state(root, {name: identity for name, identity in entries.items() if identity is not None})}
                 for root, entries in extras.items() if root != project or not (root / ".git").exists()]
    result = {"schema": SCHEMA, "origin": origin, "git": state,
              "source_files": sources, "source_sha256": _digest(sources),
              "changes_since_build": changes, "extra_source_git": extra_git,
              "omitted_source_count": len(omitted)}
    return validate_provenance(result)


def _write_atomic(path, value):
    path = Path(path)
    stream = tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=path.parent,
                                         prefix="." + path.name + ".", delete=False)
    temporary = Path(stream.name)
    try:
        with stream:
            stream.write(json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _input_identities(inputs):
    result = []
    for item in inputs:
        if isinstance(item, dict):
            if (set(item) != {"name", "bytes", "sha256"} or not _plain_name(item["name"])
                    or not _identity({key: item[key] for key in ("bytes", "sha256")})):
                raise ValueError("invalid prehashed provenance input")
            result.append(dict(item))
        else:
            path = Path(item)
            result.append({"name": path.name, **file_identity(path)})
    return result


def _output_record(directory, names, producer, inputs, extra_sources):
    files = {}
    for name in names:
        if not _plain_name(name) or name == "provenance.json":
            raise ValueError("bundle provenance requires distinct plain output file names")
        if name in files:
            raise ValueError("bundle provenance requires distinct plain output file names")
        files[name] = file_identity(directory / name)
    if not files:
        raise ValueError("output provenance requires at least one output")
    provenance = capture_provenance(extra_sources=extra_sources) if producer is None else validate_provenance(producer)
    return {"schema": OUTPUT_SCHEMA, "producer": provenance, "outputs": files, "inputs": _input_identities(inputs)}


def write_output_provenance(path, *, producer=None, inputs=(), extra_sources=()):
    """Write a byte-bound sidecar for an existing output, preserving its format."""
    path = Path(path)
    record = _output_record(path.parent, [path.name], producer, inputs, extra_sources)
    destination = path.with_name(path.name + ".provenance.json")
    _write_atomic(destination, record)
    return record


def write_bundle_provenance(directory, names, *, producer=None, inputs=(), extra_sources=()):
    """Write one provenance record binding all named outputs of a generation."""
    directory = Path(directory)
    record = _output_record(directory, names, producer, inputs, extra_sources)
    _write_atomic(directory / "provenance.json", record)
    return record
