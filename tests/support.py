"""Repository-owned scientific fixtures used by source-distribution tests."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LEGACY_X3 = ROOT / "tests" / "fixtures" / "legacy_x3"
