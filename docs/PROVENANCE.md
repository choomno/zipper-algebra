# Source provenance of computation outputs

New Python computation outputs record the zipper-algebra source revision and
the exact source-file bytes observed by the producing process. The solver's
revision and executable hash remain separate from the Python implementation.

## Where the records live

| Output | Provenance location |
| --- | --- |
| Discovery result | `result.json`, under `producer`. |
| Numeric solver run | Its database record, under `request.producer`; the record also identifies every input/output file by hash. `export-sample` restores it as `record.json`. |
| Final discovery bundle | `provenance.json`, covering the five mathematical products by filename, size and SHA-256. The result's generation hash includes this metadata file. |
| Bound skeleton, lift or certificate JSON | The `provenance` field in the version-2 envelope. Earlier JSON envelope versions are rejected. |
| Standalone emitted solver input or saved certificate, atoms, lift or skeleton | An adjacent `<filename>.provenance.json` sidecar. |
| Composition and verification JSON reports | `producer`, with hashes of the supplied input files. |
| Standalone GAP certificate script | A comment header records the producing source and input-certificate hash. |
| Experimental Python file outputs | Adjacent sidecars, including source and parent-input hashes where the calculation reads saved data. |

The solver-input, certificate-text and experimental pickle formats are unchanged.
Keep a standalone file together with its sidecar. A managed bundle carries one
shared record instead of a separate sidecar for every intermediate file.

Description reports requested with `--json` also identify the implementation.
Plain console output, arbitrary shell redirects and compiled documentation are outside these Python file-output hooks.

## What is captured

`producer` is a `zipper-source-provenance-v1` record:

- `git.commit`: the full source commit, or `null` when unavailable.
- `git.dirty`: whether the captured checkout has recorded local changes;
  `null` means that Git state could not be determined.
- `git.changes`: staged and working-tree status, logical paths and current
  byte hashes. Renames also record the preceding path. Deleted files have no
  current byte hash. Untracked implementation and declared build-source files
  are included.
- `source_files`: byte lengths and SHA-256 hashes of the maintained Python
  sources and, in a source checkout, the known build files. Experimental
  writers additionally identify their recipe family and shared helpers.
- `source_sha256`: a deterministic digest of that source-file map.
- `omitted_source_count`: how many ignored source candidates were excluded.
  A positive count identifies a partial snapshot without naming those files.
- `origin`: `checkout`, `build` or `unversioned`.
- `changes_since_build`: added, modified or deleted package, declared recipe
  and build sources relative to the embedded build origin.
- `extra_source_git`: a separate source origin when a checkout recipe uses an
  installed package from another revision.

Capture does not include source contents, patches, conversation records,
remote URLs or absolute machine paths. Git-ignored working files are excluded
when the checkout's Git ignore policy is available. A dirty-source
record identifies the changed bytes but does not preserve those bytes for
reconstruction. Commit the intended implementation before a publishable run,
or preserve its source snapshot separately.

Git cleanliness describes the checkout's status. It does not mean that ignored
sources were captured; inspect `omitted_source_count` as well.

The source list is deliberately conservative: it covers all maintained package
modules, including ones a particular operation might not import. It is not a
minimal dependency analysis. The separate reproduction inventory also includes
tests and build support; those files are not all computation inputs.

The older experimental sidecars identify Python sources and supplied inputs.
Their raw solver-output hooks do not yet record the complete external solver
configuration. Maintained discovery records include the executable digest,
declared revision, arguments and selected I/O route.

## Installed distributions

Wheels and source distributions embed `_build_provenance.json`. Outside a source
checkout, outputs retain that original build's commit and Git state while
rehashing the installed Python files. Interpret `git.dirty` as the build
checkout's state and inspect `changes_since_build` for later installed-file
changes. The build-file hashes remain available in the embedded build record.

An extracted, unchanged source distribution preserves its original provenance
when rebuilt, including when it lies inside an unrelated Git checkout. A source
archive whose recorded computation or build sources were edited is rejected for rebuilding with the
old provenance; make those changes in a Git checkout to obtain a new source
origin. Unversioned exports without embedded metadata record unknown Git state
and their actual source hashes.

Release builds from a Git checkout require a successful ignore-policy check.
Ignored symlinked payloads are omitted; other symlinked payloads are rejected.
The build stops if it cannot classify local files for release.

## Cache and publication behaviour

The solver request includes its producer snapshot. A different commit, recorded
local state or source-file contents produces a different cache key. Existing
requests without provenance are not silently relabelled or reused. Cached
records retain their original producer. The solver executable digest, declared
source revision and selected interface (`files` or `sqlite`) are part of the
request identity. File-based and native SQLite runs therefore have distinct
cache keys even when their mathematical outputs agree. Both routes retain the
same input, basis, proof and process-output bytes for inspection with
`export-sample`; see the [reproduction guide](X3_REPRODUCTION.md#read-the-result).

Discovery captures one snapshot for its samples and final products. The package
sources are rehashed before publication. A detected source-file change aborts
publication and yields an inconclusive result; it cannot retain an earlier
successful result as the current selection. These checks observe files on disk;
they do not inspect code previously loaded into a long-lived interpreter.

Managed generations are published atomically. A standalone data file and its
sidecar are two files, so their pair cannot be replaced atomically. Check the
sidecar's output hash before associating its metadata with a file. For example:

```python
import json
from pathlib import Path
from zipper_algebra.provenance import file_identity, validate_provenance

sidecar = Path("recomposed.txt.provenance.json")
record = json.loads(sidecar.read_text())
validate_provenance(record["producer"])
for name, expected in record["outputs"].items():
    assert file_identity(sidecar.parent / name) == expected
```

The hashes detect mismatched bytes. Provenance is a process's record of what it
observed, not a signature or a mathematical verification result. Certificate
verification still checks the complete algebraic identity separately.
