# Mathematical reference

Read the [compiled PDF](mathematical-specification.pdf), or edit its
[LaTeX source](main.tex).

From the repository root, compile with `python ci/build_docs.py`. The build
requires `pdflatex` with the standard AMS, geometry and hyperref packages.
Auxiliary files are temporary; the resulting PDF is kept here for readers.
GitHub Actions also provides its independently compiled copy as a workflow
artifact. This is the maintained software reference; the root-level
reconstruction manuscript is a separate historical document.
