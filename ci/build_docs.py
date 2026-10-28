"""Compile the maintained mathematical reference with temporary TeX auxiliaries."""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "docs/mathematical-specification/main.tex"
OUTPUT = SOURCE.parent / "mathematical-specification.pdf"


def source_epoch():
    """The source revision date is unaffected by later PDF-only commits."""
    try:
        result = subprocess.run(
            ["git", "log", "-1", "--format=%ct", "--", str(SOURCE.relative_to(ROOT))],
            cwd=ROOT, text=True, capture_output=True,
        )
    except OSError:
        return None
    value = result.stdout.strip()
    return value if result.returncode == 0 and value.isdecimal() else None


def compile_reference(output=OUTPUT, engine="pdflatex", engine_args=()):
    output = Path(output).resolve()
    if output.suffix.lower() != ".pdf":
        raise ValueError("the output path must end in .pdf")
    environment = dict(os.environ)
    epoch = source_epoch()
    if epoch is not None:
        environment.setdefault("SOURCE_DATE_EPOCH", epoch)
        environment.setdefault("FORCE_SOURCE_DATE", "1")
    with tempfile.TemporaryDirectory(prefix="zipper-math-reference-") as directory:
        command = [engine, *engine_args, "-no-shell-escape", "-halt-on-error",
                   "-interaction=nonstopmode", "-file-line-error",
                   f"-output-directory={directory}", str(SOURCE.relative_to(ROOT))]
        for number in (1, 2):
            result = subprocess.run(command, cwd=ROOT, env=environment, text=True,
                                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            if result.returncode:
                raise RuntimeError(f"LaTeX pass {number} failed:\n{result.stdout}")
        log = (Path(directory) / "main.log").read_text(encoding="utf-8", errors="replace")
        problems = (r"undefined references", r"Reference .* undefined",
                    r"Label\(s\) may have changed", r"Rerun to get", r"Overfull \\[hv]box")
        if any(re.search(pattern, log, flags=re.I) for pattern in problems):
            raise RuntimeError(f"The compiled reference has unresolved layout/reference warnings:\n{log}")
        pdf = Path(directory) / "main.pdf"
        if not pdf.is_file() or not pdf.read_bytes().startswith(b"%PDF-"):
            raise RuntimeError("LaTeX did not produce a PDF")
        output.parent.mkdir(parents=True, exist_ok=True)
        descriptor, pending = tempfile.mkstemp(prefix=".compiled-", suffix=".pdf", dir=output.parent)
        os.close(descriptor)
        pending = Path(pending)
        try:
            shutil.copyfile(pdf, pending)
            pending.chmod(0o644)
            pending.replace(output)
        finally:
            pending.unlink(missing_ok=True)
    print(f"Compiled mathematical reference: {output}")
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--engine", default="pdflatex")
    parser.add_argument("--engine-arg", action="append", default=[])
    args = parser.parse_args()
    try:
        compile_reference(args.output, args.engine, args.engine_arg)
    except (OSError, RuntimeError, ValueError) as error:
        parser.exit(1, f"Could not compile the mathematical reference: {error}\n")
