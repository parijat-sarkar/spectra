"""
PAM scanning and protospacer enumeration on a genomic + strand sequence.

A "guide" here is a dict with keys:
    sequence          : 20 nt protospacer on the guide strand (5'->3')
    pam               : 3 nt PAM on the guide strand
    context           : 30 nt context: 4 upstream + 20 protospacer + 3 PAM + 3 downstream
    start_global      : 1-based lowest + strand coord of the protospacer
    orientation       : 'sense' or 'antisense'

Coordinates note:
    - Input sequence `plus_seq` is the + strand, 1-based genomic coord
      `seq_start` corresponds to plus_seq[0].
    - `start_global` is the lowest + strand coord of the protospacer,
      matching Beagle's "sgRNA Sequence Start Pos. (global)" column.
"""

from __future__ import annotations

import re

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "N": "[ACGT]",
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]",
    "B": "[CGT]", "D": "[AGT]", "H": "[ACT]", "V": "[ACG]",
}

_COMP = str.maketrans("ACGTN", "TGCAN")


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def compile_pam(pam_pattern: str) -> re.Pattern:
    expanded = "".join(IUPAC[c.upper()] for c in pam_pattern)
    return re.compile(expanded)


def pam_to_rc_pattern(pam_pattern: str) -> str:
    """Return the + strand string pattern whose RC matches pam_pattern."""
    return reverse_complement(pam_pattern.upper()).replace("T", "T")  # noop but explicit


def find_guides(
    plus_seq: str,
    seq_start: int,
    pam_pattern: str = "NGG",
    protospacer_len: int = 20,
) -> list[dict]:
    """
    Scan both strands of plus_seq for PAM matches and return all protospacers.

    plus_seq      : + strand DNA sequence, upper-case.
    seq_start     : 1-based genomic coord of plus_seq[0].
    pam_pattern   : PAM on the guide strand, e.g., "NGG", "NGN", "NG".
    protospacer_len: length of protospacer (20 for SpCas9).

    Returns guides with the invariant that start_global is always the
    lowest + strand coordinate of the protospacer, regardless of orientation.
    """
    pam_pattern = pam_pattern.upper()
    pam_len = len(pam_pattern)
    plus_re = compile_pam(pam_pattern)
    minus_re = compile_pam(reverse_complement(pam_pattern))

    guides: list[dict] = []
    L = len(plus_seq)

    # --- SENSE strand ---
    # protospacer at [i : i+P], PAM at [i+P : i+P+pam_len]
    # require: i >= 4 (for 4nt upstream context) and i+P+pam_len+3 <= L (for 3nt downstream context)
    for i in range(4, L - protospacer_len - pam_len - 3 + 1):
        pam = plus_seq[i + protospacer_len : i + protospacer_len + pam_len]
        if not plus_re.fullmatch(pam):
            continue
        protospacer = plus_seq[i : i + protospacer_len]
        context = plus_seq[i - 4 : i + protospacer_len + pam_len + 3]
        if len(context) != 4 + protospacer_len + pam_len + 3:
            continue  # guard against boundary cases
        guides.append({
            "sequence": protospacer,
            "pam": pam,
            "context": context,
            "start_global": seq_start + i,  # 1-based lowest + strand coord of protospacer
            "orientation": "sense",
        })

    # --- ANTISENSE strand ---
    # PAM (on + strand, reverse-complemented) at [i : i+pam_len],
    # protospacer (on + strand, will be RC'd) at [i+pam_len : i+pam_len+P].
    # Require: i >= 3 (for 3nt downstream context on - strand = upstream on +)
    # and i+pam_len+P+4 <= L (for 4nt upstream context on - strand = downstream on +)
    for i in range(3, L - pam_len - protospacer_len - 4 + 1):
        pam_on_plus = plus_seq[i : i + pam_len]
        if not minus_re.fullmatch(pam_on_plus):
            continue
        proto_on_plus = plus_seq[i + pam_len : i + pam_len + protospacer_len]
        # context on - strand = RC of plus[i-3 : i+pam_len+P+4]
        plus_context = plus_seq[i - 3 : i + pam_len + protospacer_len + 4]
        if len(plus_context) != 4 + protospacer_len + pam_len + 3:
            continue
        context = reverse_complement(plus_context)
        guides.append({
            "sequence": reverse_complement(proto_on_plus),
            "pam": reverse_complement(pam_on_plus),
            "context": context,
            "start_global": seq_start + i + pam_len,  # lowest + strand coord of protospacer
            "orientation": "antisense",
        })

    return guides


def protospacer_pos_to_plus_coord(start_global: int, orientation: str, pos_1based: int, protospacer_len: int = 20) -> int:
    """
    Map a 1-based protospacer position (1 = 5' end of the guide) to its + strand coord.
    """
    if orientation == "sense":
        return start_global + pos_1based - 1
    elif orientation == "antisense":
        return start_global + protospacer_len - pos_1based
    else:
        raise ValueError(f"unknown orientation: {orientation}")
