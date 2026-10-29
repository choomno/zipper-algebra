# f4ncgb setup and parser patch

Fresh solver computations use the external [f4ncgb project](https://gitlab.sai.jku.at/f4ncgb/f4ncgb).
The solver is developed by Maximilian Heisinger and Clemens Hofstadler; cite
their [CASC 2025 paper](https://doi.org/10.1007/978-3-032-09645-6_5),
*f4ncgb: High Performance Gröbner Basis Computations in Free Algebras*,
Lecture Notes in Computer Science **16235**, pp. 79–97.

Checking the saved certificate, recomposing it from the saved symbolic lift,
and replaying its retained proof require neither f4ncgb nor this patch; see the
[reproduction guide](../../docs/X3_REPRODUCTION.md).

The bundled [parser patch](msolve-parser-fix.patch) targets upstream commit
`f152bb0d040d7947d4296f9d68ca849ee135629c` (v1.1.3). This patch baseline is
separate from the original research binary, whose recovered settings describe
v1.0.1 without establishing its exact executable bytes.
That revision was still upstream `main` when checked on 27 October 2026.
The production parser changes extend recovered local commit
`0301fb910fdfed3270cf95fe2742d34d0683cca7` with single-variable CRLF header
handling. The test harness has also been strengthened with isolated temporary
files, exact LF/CRLF block comparisons and declared/undeclared exponent-zero
controls. The bundled parser and tests match
[the maintained fix commit](https://github.com/choomno/f4ncgb/commit/789504d7c7da6af6946faecb73421b037766a5b8).
The [parser-fix PR](https://github.com/choomno/f4ncgb/pull/1) is merged in the personal
GitHub copy; the official source remains on JKU’s GitLab.

The patch fixes CRLF handling at the final input line and in single-variable
and multi-block headers, and rejects missing or undeclared variable names when they are read, including
undeclared names raised to exponent zero. Older f4ncgb versions could silently
introduce unintended variables from malformed input. Correct LF input with
declared variables does not inherently require this patch. It is a parser fix;
it does not establish completion or correctness of a Gröbner computation.

## Apply to the recorded revision

Run from the zipper-algebra checkout; choose a destination outside this repository:

```sh
patch_file="$PWD/patches/f4ncgb/msolve-parser-fix.patch"
git clone https://gitlab.sai.jku.at/f4ncgb/f4ncgb.git ../f4ncgb-source
git -C ../f4ncgb-source checkout --detach f152bb0d040d7947d4296f9d68ca849ee135629c
git -C ../f4ncgb-source apply --check "$patch_file"
git -C ../f4ncgb-source apply "$patch_file"
```

The patch changes `src/parser.cpp`, adds parser regressions in `src/test_parse.cpp`,
and updates upstream `NEWS.md` and `README.md`. The upstream MIT
[license](LICENSE) is included beside the patch.

## Build and test

The recorded upstream build needs a C++20 compiler, CMake, Boost (at least 1.80),
FLINT and GMP; mimalloc is optional. Consult the pinned upstream README for its
supported dependency versions and platform details. After installing them:

```sh
cmake -S ../f4ncgb-source -B ../f4ncgb-source/build-release -DCMAKE_BUILD_TYPE=Release
cmake --build ../f4ncgb-source/build-release
../f4ncgb-source/build-release/f4ncgb_test
../f4ncgb-source/build-release/f4ncgb --version
```

Build the default target: the target named `f4ncgb` is the library, while the
executable target is `f4ncgb_bin`.

All 31 tests passed in the recorded
[Ubuntu 24.04 GitHub build](https://github.com/choomno/f4ncgb/actions/runs/40625473978)
at `789504d7c7da6af6946faecb73421b037766a5b8` and locally on Fedora 44. The
historical 23 accepted-input output comparisons used iteration caps and do not
assert completed Gröbner bases.

Pass the resulting executable explicitly to `zipper-algebra discover --solver`.
The installed command takes output/cache locations and limits directly; it does
not use the old research workspace configuration. The reproduction guide gives
the complete command and its smaller three-strand example, as well as
solver-free certificate checking. The patch is reviewed in the personal
GitHub copy; no upstream issue or merge request is being submitted.


## Optional direct SQLite I/O

Discovery stores sample data in SQLite with either solver interface. The default
`--solver-io files` works with the parser-patched revision above: each active
solver call uses temporary input, basis and proof files, then their bytes are
stored in the database and the temporary files are removed.

The optional `--solver-io sqlite` route requires the SQLite extension in the
maintained f4ncgb copy. This extension is separate from the bundled parser patch.
The solver reads the input and writes the basis and nonexpanded proof directly
in the database. Python still checks every proof identity before recording a
validated sample. The same presentation, arithmetic and certificate checks apply
to both routes.

The [SQLite extension](https://github.com/choomno/f4ncgb/pull/2) is pinned to
`cb37d6570a42cf2f786cc96180b18456ae6d1bcf`. Clone that source into a
separate checkout; the original upstream remote does not contain this extension:

```sh
git clone https://github.com/choomno/f4ncgb.git ../f4ncgb-sqlite-source
git -C ../f4ncgb-sqlite-source fetch origin cb37d6570a42cf2f786cc96180b18456ae6d1bcf
git -C ../f4ncgb-sqlite-source checkout --detach cb37d6570a42cf2f786cc96180b18456ae6d1bcf
```

Install the additional SQLite development dependency (`sqlite-devel` on Fedora or `libsqlite3-dev` on
Debian/Ubuntu), then build with the feature enabled:

```sh
cmake -S ../f4ncgb-sqlite-source -B ../f4ncgb-sqlite-source/build-sqlite \
  -DCMAKE_BUILD_TYPE=Release -DENABLE_SQLITE=ON
cmake --build ../f4ncgb-sqlite-source/build-sqlite
../f4ncgb-sqlite-source/build-sqlite/f4ncgb_test
../f4ncgb-sqlite-source/build-sqlite/f4ncgb --features
```

The feature report must contain `"sqlite": true`. Builds with
`-DENABLE_SQLITE=OFF` (the default) retain file I/O without linking SQLite.
To use the enabled executable, add `--solver-io sqlite` to the discovery command,
set `--solver` to this executable and set `--solver-revision` to its actual source
commit (`git -C ../f4ncgb-sqlite-source rev-parse HEAD`). Do not use the older parser-only
revision label for the extended binary.

The adapter removes intermediate algebra files; SQLite may create temporary
journal files during transactions. Completed samples use the same database
schema and `export-sample` inspection command with either solver interface.
