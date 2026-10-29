"""Check local Markdown links in current public usage and reference documentation."""
from pathlib import Path
import re
from urllib.parse import unquote, urlsplit


if __name__ == "__main__":
    root = Path(__file__).resolve().parents[1]
    files = ([root / name for name in ("README.md", "INDEX.md", "WORKFLOW.md")]
             + list((root / "docs").rglob("*.md"))
             + [path for family in ("experiments", "examples", "patches", "singular_run")
                for path in (root / family).rglob("README.md")])
    checked = 0
    for source in files:
        text = re.sub(r"```.*?```", "", source.read_text(encoding="utf-8"), flags=re.S)
        for match in re.finditer(r"\[[^\]]*\]\(([^\s)]+)\)", text):
            url = urlsplit(match.group(1))
            if url.scheme or url.netloc or not url.path:
                continue
            target = source.parent / unquote(url.path)
            assert target.exists(), f"{source.relative_to(root)}: missing {match.group(1)}"
            checked += 1
    print(f"{checked} local documentation links resolve")
