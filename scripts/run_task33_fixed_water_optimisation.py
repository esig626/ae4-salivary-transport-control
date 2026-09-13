"""Run one of the two checkpoint-separated Task 33 stages."""
from __future__ import annotations

from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from modern_full_model.task33_fixed_water_optimisation import main


if __name__ == "__main__":
    main()
