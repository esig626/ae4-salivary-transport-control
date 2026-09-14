"""Execute the bounded Task 30 NHE1 repair verification."""
from __future__ import annotations

import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from modern_full_model.task30_nhe1_repair import run


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, sort_keys=True))
