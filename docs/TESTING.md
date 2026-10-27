# Tests and continuous integration

Use Python 3.11 or newer. Install the package before running the source tests:

```sh
python -m pip install -e .
python -m unittest discover -s tests -v
python -m unittest tools.cofactor.tests_proof_skeleton -v
```

The maintained tests need no solver, personal configuration, sibling research
folders, populated cache or network access after installation. They cover exact
and modular arithmetic, basis parsing, presentation contracts, proof expansion,
certificate validation, contextual composition, repeated reduction steps, exact
interpolation, held-out validation, inconclusive outcomes, package imports and
CLI exit codes. Tests use temporary directories for output.

The historical comparisons use the archived implementation and immutable
fixtures under `tests/fixtures/legacy_x3/`. Their hashes and original records are
listed in the adjacent manifest. The retained large free-algebra certificate is
a compatibility fixture; tests do not identify it as the latest mathematical proof.

## Distribution checks

Build both formats, then inspect their contents:

```sh
python -m pip install build
python -m build
python ci/check_distribution.py dist
```

In a fresh environment, install the built wheel and run
`python ci/check_installed.py` from the checkout. That script imports the package
and runs commands from a separate temporary working directory in isolated Python
processes. It rejects an editable installation for this check. It distinguishes
a valid certificate, a false identity and malformed input.

`python ci/check_sdist.py dist` extracts the source distribution, rebuilds its
wheel and runs the source regressions. Run it in a disposable environment: it
replaces that environment's installed development package with the rebuilt wheel.

The wheel contains only the maintained package and distribution metadata. The
source archive includes public documentation, tests and selected regression
fixtures. Private histories, recovered manuscripts, machine configuration,
caches and experimental result collections are excluded.

## GitHub Actions

The workflow runs on pushes, pull requests and manual dispatch. Its matrix
checks Linux on Python 3.11, 3.12, 3.13 and 3.14, plus Windows on Python 3.13.
It builds the wheel and source archive, tests the installed wheel, runs the
regressions and rebuilds the source archive. A separate job checks the optional
Sympy factorization extra. Documentation checks validate local links and compile
the mathematical specification.

Jobs have time limits, read-only repository permissions and no publication step.
The workflow does not require project secrets or an algebra solver. A local test
pass and a GitHub-hosted run are separate evidence; the validation record states
which was actually observed.

## Experimental integrations

The older `python -m tools.tests`, solver sampling, dimension and trace scripts
retain separate dependency and completion requirements. They are outside the
maintained CI boundary. Large discovery computations need explicit tool versions
and completion evidence before their mathematical conclusions can be accepted.

See the [mathematical specification](mathematical-specification/main.tex) for
the coefficient domains, assumptions and conclusions of the maintained operations.
