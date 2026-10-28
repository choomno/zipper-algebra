# f4ncgb parser patch

Fresh solver computations use the external [f4ncgb project](https://gitlab.sai.jku.at/f4ncgb/f4ncgb).
Checking the saved certificate, recomposing it from the saved symbolic lift,
and replaying its retained proof require neither f4ncgb nor this patch; see the
[reproduction guide](../../docs/X3_REPRODUCTION.md).

The recovered research environment used the [parser patch](msolve-parser-fix.patch)
against upstream commit `f152bb0d040d7947d4296f9d68ca849ee135629c` (v1.1.3).
That revision was still upstream `main` when checked on 27 October 2026.
The patch payload matches recovered local commit
`0301fb910fdfed3270cf95fe2742d34d0683cca7`.

The patch fixes CRLF handling at the final input line and in block headers,
and rejects missing or undeclared variable names when they are read, including
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

The recovered record reports 31 sanitizer tests and 23 accepted-input output
comparisons. During the current cleanup the patch was checked against the exact
base source; the compiler and CMake are absent from this environment, so those
historical build results have not been repeated here. The 23 output comparisons
used iteration caps and do not assert completed Gröbner bases.

The historical solver driver also needs workspace configuration and path repairs;
applying this patch alone does not make that driver portable. Its supported
alternative for auditing the result is the solver-free route in the reproduction
guide. The upstream issue and merge request have not been filed as part of this
cleanup.
