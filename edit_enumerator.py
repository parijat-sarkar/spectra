"""
Enumerate partial-edit outcomes for a protospacer.

For ABE: every A in the window is a potential edit site (A->G on the guide
strand). For CBE: every C (C->T).

Given N editable bases in the window, we enumerate all 2^N - 1 non-empty
subsets. Each subset is one "outcome" where exactly those sites are edited
and the others are not. "Full edit" (all N sites) is one of the 2^N - 1.
"No edit" is excluded because it's the reference.

For a guide with 0 editable bases in the window, no outcomes are returned.
"""

from __future__ import annotations

from itertools import combinations
from typing import Iterable


EDITOR_SPECS = {
    "ABE": {"target_base": "A", "product_base": "G", "edit_type": "A-G"},
    "CBE": {"target_base": "C", "product_base": "T", "edit_type": "C-T"},
}


def editable_positions(protospacer: str, editor: str, window: tuple[int, int]) -> list[int]:
    """
    Return 1-based protospacer positions where the editable base appears
    inside the window.

    window = (start, end), both inclusive, 1-based.
    """
    spec = EDITOR_SPECS[editor]
    w_start, w_end = window
    target = spec["target_base"]
    return [
        p for p in range(w_start, w_end + 1)
        if 1 <= p <= len(protospacer) and protospacer[p - 1] == target
    ]


def enumerate_partial_outcomes(
    protospacer: str,
    editor: str,
    window: tuple[int, int],
    max_edits: int | None = None,
) -> list[list[int]]:
    """
    Return a list of edit-position subsets. Each subset is a list of 1-based
    protospacer positions that are edited simultaneously.

    If max_edits is set, only subsets with <= max_edits positions are returned.
    Subsets are emitted in ascending order of subset size (1, 2, ..., N),
    and within each size in combinatorial order.

    The full-edit outcome (all N positions) is included. The no-edit
    outcome (empty subset) is excluded.
    """
    positions = editable_positions(protospacer, editor, window)
    n = len(positions)
    if n == 0:
        return []
    max_k = n if max_edits is None else min(n, max_edits)
    outcomes: list[list[int]] = []
    for k in range(1, max_k + 1):
        for combo in combinations(positions, k):
            outcomes.append(list(combo))
    return outcomes


def outcome_cardinality(n_editable: int, max_edits: int | None = None) -> int:
    """Number of outcomes for an N-editable-site guide."""
    if n_editable == 0:
        return 0
    max_k = n_editable if max_edits is None else min(n_editable, max_edits)
    from math import comb
    return sum(comb(n_editable, k) for k in range(1, max_k + 1))
