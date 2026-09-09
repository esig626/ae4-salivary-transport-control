"""Run the isolated WT-only NABS-S4 by H2.5 hydraulic diagnostic hierarchy."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence

from .absolute_hydraulic_diagnostic import run_absolute_hydraulic_diagnostic


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-directory",
        type=Path,
        default=Path("results/13B_modern_full_model"),
    )
    parser.add_argument(
        "--native-output-directory",
        type=Path,
        default=Path("results/13B_modern_full_model"),
        help="read-only directory containing the committed native WT freeze",
    )
    parser.add_argument("--workers", type=int, default=9)
    parser.add_argument("--seed", type=int, default=13601)
    parser.add_argument("--max-nfev", type=int, default=3000)
    args = parser.parse_args(argv)
    if args.workers < 1:
        parser.error("--workers must be positive")
    if args.max_nfev < 1:
        parser.error("--max-nfev must be positive")
    summary = run_absolute_hydraulic_diagnostic(
        output=args.output_directory,
        native_output=args.native_output_directory,
        workers=args.workers,
        seed=args.seed,
        max_nfev=args.max_nfev,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

