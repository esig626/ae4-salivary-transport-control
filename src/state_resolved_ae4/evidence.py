"""Mechanical evidence-freeze and holdout boundary for Task 13."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
from pathlib import Path
from typing import Iterable, Mapping


EXPECTED_FREEZE_HASHES: Mapping[str, str] = {
    "analysis/13_state_resolved_ae4/evidence_freeze.md": (
        "48022d5f69c4459edf6c9cceb914dd5ce092ffe73ffaffdbf9ceedc578005631"
    ),
    "analysis/13_state_resolved_ae4/question_tree.md": (
        "34343e076d46b00ce6eadc7a419688747b13917afce576611a2149c79aba019e"
    ),
}


CAL_T_TARGET_IDS = frozenset(
    {
        "E16-01",
        "E16-02",
        "E16-03",
        "E16-04",
        "E16-05",
        "E16-06",
        "E16-07",
        "E16-09",
        "E16-12",
        "E21-01",
        "E21-02",
        "E21-03",
        "E21-04",
        "E21-05",
        "E25-01",
        "E25-02",
        "E25-03",
        "E25-04",
    }
)

CAL_WT_TARGET_IDS = frozenset(
    {
        "E15-01",
        "E15-05",
        "E15-07_WT",
        "E15-08_WT",
        "E15-09_WT",
    }
)

STAGE_A_HELDOUT_IDS = frozenset(
    {
        "E15-02",
        "E15-03",
        "E15-06",
        "E15-07_KO",
        "E15-08_KO",
        "E15-09_KO",
        "E15-10",
        "E15-13",
        "E15-16",
    }
)

STRICT_SECRETION_HOLDOUT_IDS = frozenset({"E15-02", "E15-03"})


@dataclass(frozen=True)
class FreezeBoundary:
    """Hashes that must be verified before any calibration is run."""

    evidence_hash: str
    question_tree_hash: str
    status: str
    allowed_target_ids: tuple[str, ...]
    forbidden_target_ids: tuple[str, ...]

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_evidence_freeze(repository_root: str | Path) -> FreezeBoundary:
    """Verify the exact curator-frozen evidence before calibration.

    This deliberately fails closed if either document changes.  A scientific
    amendment needs a new reviewed hash rather than silent calibration drift.
    """

    root = Path(repository_root)
    actual: dict[str, str] = {}
    for relative, expected in EXPECTED_FREEZE_HASHES.items():
        path = root / relative
        if not path.is_file():
            raise RuntimeError(f"required frozen evidence file is absent: {relative}")
        value = _sha256(path)
        actual[relative] = value
        if value != expected:
            raise RuntimeError(
                f"evidence freeze hash changed for {relative}: {value} != {expected}"
            )
    evidence_text = (root / "analysis/13_state_resolved_ae4/evidence_freeze.md").read_text(
        encoding="utf-8"
    )
    question_text = (root / "analysis/13_state_resolved_ae4/question_tree.md").read_text(
        encoding="utf-8"
    )
    if "Status: FROZEN before Task 13 fitting" not in evidence_text:
        raise RuntimeError("evidence document does not declare the frozen boundary")
    if "Frozen with the evidence table before fitting" not in question_text:
        raise RuntimeError("question tree does not declare the frozen boundary")
    return FreezeBoundary(
        evidence_hash=actual["analysis/13_state_resolved_ae4/evidence_freeze.md"],
        question_tree_hash=actual["analysis/13_state_resolved_ae4/question_tree.md"],
        status="verified_before_calibration",
        allowed_target_ids=tuple(sorted(CAL_T_TARGET_IDS | CAL_WT_TARGET_IDS)),
        forbidden_target_ids=tuple(sorted(STAGE_A_HELDOUT_IDS)),
    )


def assert_calibration_targets_allowed(target_ids: Iterable[str]) -> tuple[str, ...]:
    """Reject any calibration API call containing a held-out or unknown ID."""

    target_tuple = tuple(target_ids)
    allowed = CAL_T_TARGET_IDS | CAL_WT_TARGET_IDS
    forbidden = [target for target in target_tuple if target not in allowed]
    if forbidden:
        raise PermissionError(
            "calibration target is not CAL-T/CAL-WT under the frozen evidence: "
            + ", ".join(sorted(forbidden))
        )
    return target_tuple

def assert_stage_b_ionic_targets_only(target_ids: Iterable[str]) -> tuple[str, ...]:
    """Enforce the lead's narrow Stage-B release: resting Cl and pH only."""

    target_tuple = tuple(target_ids)
    permitted = {"E15-06", "E15-07_KO"}
    forbidden = set(target_tuple) - permitted
    if forbidden:
        raise PermissionError(
            "Stage B may localize only with AE4-null resting Cl/pH; prohibited: "
            + ", ".join(sorted(forbidden))
        )
    return target_tuple
