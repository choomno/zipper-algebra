# Mathematical reference

Read the [compiled PDF](mathematical-specification.pdf), or edit its
[LaTeX source](main.tex).

From the repository root, compile with `python ci/build_docs.py`. The build
requires `pdflatex` with the standard AMS, geometry and hyperref packages.
Auxiliary files are temporary; the resulting PDF is kept here for readers.
GitHub Actions also provides its independently compiled copy as a workflow
artifact. The reference begins with Bigelow’s problem and the thesis’s polynomial
untwisting family and finite-generation theorem, then specifies the implemented
algebra and certificate operations. Numbered references are included in the PDF.
Thesis citations identify the recovered manuscript version; the thesis itself
is kept outside the software distribution.
