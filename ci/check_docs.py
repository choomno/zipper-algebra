"""Check local Markdown links in current public usage and reference documentation."""
from pathlib import Path
import re
from urllib.parse import unquote, urlsplit


if __name__ == "__main__":
    root = Path(__file__).resolve().parents[1]
    files = (list(root.glob("*.md")) + list((root / "docs").rglob("*.md"))
             + list((root / "tools").rglob("README.md")))
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
