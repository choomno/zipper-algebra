"""Prepare ignored directories for optional experiment outputs."""

from pathlib import Path


def prepare_output_dir(path: Path | str) -> Path:
    """Create an output directory and preserve any existing ignore rules."""
    directory = Path(path)
    if (directory / ".git").exists():
        raise ValueError("output directory must not be a Git checkout root")
    directory.mkdir(parents=True, exist_ok=True)
    try:
        with (directory / ".gitignore").open("x", encoding="utf-8") as marker:
            marker.write("*\n!.gitignore\n")
    except FileExistsError:
        pass
    return directory
