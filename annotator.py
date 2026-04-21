"""
Transcript-aware annotation engine.

Given a transcript bundle (from ensembl_client.get_transcript_bundle) and a
list of genomic edits (on the + strand), compute:

    - per-edit HGVS coding coordinate (c.123, c.1005+7, c.-25, c.*7)
    - per-edit nucleotide change in HGVS (transcript-sense) notation
    - per-edit domain classification (CDS / UTR / Intron / Splice-donor/acceptor)
    - combined amino-acid change for all edits that land in CDS codons,
      applied simultaneously so partial edits in the same codon interact

Strand handling:
    - bundle['strand']  is +1 (forward) or -1 (reverse) for the transcript
    - HGVS nucleotide changes are reported on the transcript sense strand:
      for -1 transcripts we complement the + strand ref/alt bases.
    - Amino acid translation uses the transcript sense coding sequence.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable

# Standard genetic code
CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}
AA_THREE = {
    'A':'Ala','R':'Arg','N':'Asn','D':'Asp','C':'Cys','Q':'Gln','E':'Glu',
    'G':'Gly','H':'His','I':'Ile','L':'Leu','K':'Lys','M':'Met','F':'Phe',
    'P':'Pro','S':'Ser','T':'Thr','W':'Trp','Y':'Tyr','V':'Val','*':'Ter',
}

_COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def comp(base: str) -> str:
    return base.translate(_COMP)


# ---------------------------------------------------------------------------
# Transcript index

@dataclass
class CDSSegment:
    """One contiguous stretch of CDS within a single exon, in transcript order."""
    g_5prime: int      # + strand coord of the 5'-most CDS base of the segment
    g_3prime: int      # + strand coord of the 3'-most CDS base of the segment
    c_start: int       # 1-based c. position of g_5prime
    c_end: int         # 1-based c. position of g_3prime
    exon_tx_idx: int   # 0-based exon index in transcript order


@dataclass
class ExonTx:
    """Exon in transcript order."""
    idx: int            # 0-based in transcript order
    g_5prime: int
    g_3prime: int
    strand: int         # +1 or -1


@dataclass
class TranscriptIndex:
    strand: int
    exons_tx: list[ExonTx]
    cds_segments: list[CDSSegment]
    total_cds_len: int
    cds_start_g: int    # genomic + strand coord of the first CDS base (start codon)
    cds_end_g: int      # genomic + strand coord of the last CDS base (before stop)
    # Note: Ensembl's Translation.start/end are + strand coords of first/last codon bases.
    # On - strand, Translation.start > Translation.end is possible but usually start<=end as
    # genomic coords; we always treat them as min/max.


def build_transcript_index(bundle: dict) -> TranscriptIndex:
    strand = bundle["strand"]
    exons = sorted(bundle["exons"], key=lambda e: e["start"])
    if strand == -1:
        exons = list(reversed(exons))  # now in transcript 5'->3' order

    exons_tx: list[ExonTx] = []
    for i, e in enumerate(exons):
        if strand == 1:
            exons_tx.append(ExonTx(i, e["start"], e["end"], strand))
        else:
            exons_tx.append(ExonTx(i, e["end"], e["start"], strand))

    cds_segments: list[CDSSegment] = []
    total = 0
    cds_g_min = bundle.get("cds_start")
    cds_g_max = bundle.get("cds_end")
    if cds_g_min is not None and cds_g_max is not None:
        if cds_g_min > cds_g_max:
            cds_g_min, cds_g_max = cds_g_max, cds_g_min
        for ex in exons_tx:
            # Exon genomic range on + strand:
            g_lo = min(ex.g_5prime, ex.g_3prime)
            g_hi = max(ex.g_5prime, ex.g_3prime)
            seg_lo = max(g_lo, cds_g_min)
            seg_hi = min(g_hi, cds_g_max)
            if seg_lo > seg_hi:
                continue
            length = seg_hi - seg_lo + 1
            if strand == 1:
                g5, g3 = seg_lo, seg_hi
            else:
                g5, g3 = seg_hi, seg_lo
            cds_segments.append(
                CDSSegment(
                    g_5prime=g5,
                    g_3prime=g3,
                    c_start=total + 1,
                    c_end=total + length,
                    exon_tx_idx=ex.idx,
                )
            )
            total += length

    return TranscriptIndex(
        strand=strand,
        exons_tx=exons_tx,
        cds_segments=cds_segments,
        total_cds_len=total,
        cds_start_g=cds_g_min if cds_g_min else 0,
        cds_end_g=cds_g_max if cds_g_max else 0,
    )


# ---------------------------------------------------------------------------
# Per-position annotation

@dataclass
class PosAnno:
    """Per-genomic-position annotation."""
    g_pos: int
    domain: str       # 'CDS', 'UTR5', 'UTR3', 'Intron', 'Outside'
    is_splice_donor: bool = False       # intronic +1/+2
    is_splice_acceptor: bool = False    # intronic -1/-2
    c_coord: str | None = None          # HGVS c. coord like '123', '1005+7', '-25', '*7'
    cds_pos: int | None = None          # 1-based if in CDS, else None
    codon_idx: int | None = None        # 1-based amino acid index
    frame: int | None = None            # 0/1/2 position within codon


def _tx_dist(a: int, b: int, strand: int) -> int:
    """Return transcript-direction distance b - a (positive if b is downstream of a)."""
    return (b - a) if strand == 1 else (a - b)


def annotate_position(idx: TranscriptIndex, g_pos: int) -> PosAnno:
    strand = idx.strand

    # 1) First check if in any exon
    in_exon_idx = None
    for ex in idx.exons_tx:
        lo, hi = min(ex.g_5prime, ex.g_3prime), max(ex.g_5prime, ex.g_3prime)
        if lo <= g_pos <= hi:
            in_exon_idx = ex.idx
            break

    if in_exon_idx is not None:
        # Exonic. Is it in CDS?
        for seg in idx.cds_segments:
            lo, hi = min(seg.g_5prime, seg.g_3prime), max(seg.g_5prime, seg.g_3prime)
            if lo <= g_pos <= hi:
                # In CDS
                # Compute c. position
                dist = _tx_dist(seg.g_5prime, g_pos, strand)
                c_pos = seg.c_start + dist
                codon_idx = (c_pos - 1) // 3 + 1
                frame = (c_pos - 1) % 3
                return PosAnno(
                    g_pos=g_pos,
                    domain="CDS",
                    c_coord=str(c_pos),
                    cds_pos=c_pos,
                    codon_idx=codon_idx,
                    frame=frame,
                )
        # Exonic but not CDS -> UTR
        return _annotate_utr(idx, g_pos, in_exon_idx)

    # 2) Check if in any intron (between two adjacent exons in tx order)
    for i in range(len(idx.exons_tx) - 1):
        ex_prev = idx.exons_tx[i]
        ex_next = idx.exons_tx[i + 1]
        # Intron is between ex_prev's 3' and ex_next's 5'
        intron_lo = min(ex_prev.g_3prime, ex_next.g_5prime) + 1
        intron_hi = max(ex_prev.g_3prime, ex_next.g_5prime) - 1
        # Because exons are in tx order, on + strand ex_prev.g_3prime < ex_next.g_5prime;
        # on - strand ex_prev.g_3prime > ex_next.g_5prime. In both cases the intron's
        # + strand coords are (min+1) .. (max-1).
        if intron_lo <= g_pos <= intron_hi:
            return _annotate_intron(idx, g_pos, ex_prev, ex_next)

    # 3) Outside the transcript
    return PosAnno(g_pos=g_pos, domain="Outside")


def _annotate_intron(idx: TranscriptIndex, g_pos: int, ex_prev: ExonTx, ex_next: ExonTx) -> PosAnno:
    strand = idx.strand
    # Distance from 5' end of intron = tx-direction distance from ex_prev.g_3prime to g_pos
    # (5' end of intron is just past ex_prev's 3' boundary in tx direction)
    plus_5 = _tx_dist(ex_prev.g_3prime, g_pos, strand)   # 1..intron_len
    plus_from_donor = plus_5
    # Distance from 3' end of intron (from g_pos to ex_next.g_5prime)
    plus_to_acceptor = _tx_dist(g_pos, ex_next.g_5prime, strand)

    # Splice site flags (standard ±1/±2 GT..AG)
    is_donor = plus_from_donor in (1, 2)
    is_acceptor = plus_to_acceptor in (1, 2)

    # HGVS: choose +k or -k based on which end is closer (tie -> +k)
    if plus_from_donor <= plus_to_acceptor:
        suffix = f"+{plus_from_donor}"
        ref_c = _genomic_to_c_at_exon_boundary(idx, ex_prev, end="3prime")
    else:
        suffix = f"-{plus_to_acceptor}"
        ref_c = _genomic_to_c_at_exon_boundary(idx, ex_next, end="5prime")

    c_coord = f"{ref_c}{suffix}"
    return PosAnno(
        g_pos=g_pos,
        domain="Intron",
        is_splice_donor=is_donor,
        is_splice_acceptor=is_acceptor,
        c_coord=c_coord,
    )


def _genomic_to_c_at_exon_boundary(idx: TranscriptIndex, ex: ExonTx, end: str) -> str:
    """
    Return the HGVS c. coordinate (as a string) of an exon's boundary base.
    end = '3prime' (3' end of exon in tx direction) or '5prime' (5' end).
    """
    g_pos = ex.g_3prime if end == "3prime" else ex.g_5prime
    # Try to find in CDS first
    for seg in idx.cds_segments:
        if seg.exon_tx_idx != ex.idx:
            continue
        lo, hi = min(seg.g_5prime, seg.g_3prime), max(seg.g_5prime, seg.g_3prime)
        if lo <= g_pos <= hi:
            dist = _tx_dist(seg.g_5prime, g_pos, idx.strand)
            return str(seg.c_start + dist)
    # Otherwise the boundary is in a UTR: fall back to UTR encoding
    anno = _annotate_utr(idx, g_pos, ex.idx)
    return anno.c_coord or "?"


def _annotate_utr(idx: TranscriptIndex, g_pos: int, exon_tx_idx: int) -> PosAnno:
    """Compute c.-N or c.*N HGVS coord for a non-CDS exonic position."""
    strand = idx.strand
    if not idx.cds_segments:
        return PosAnno(g_pos=g_pos, domain="UTR5", c_coord="?")

    first_cds = idx.cds_segments[0]
    last_cds = idx.cds_segments[-1]

    # Is the position upstream of the start codon (5' UTR) or downstream of stop (3' UTR)?
    # Compare tx-direction distance.
    d_to_start = _tx_dist(g_pos, first_cds.g_5prime, strand)  # positive if g_pos is 5' of start
    d_to_stop = _tx_dist(last_cds.g_3prime, g_pos, strand)    # positive if g_pos is 3' of stop

    if d_to_start > 0:
        # 5' UTR: count transcript-5' bases (including exon hops) between g_pos and start codon
        n = _count_tx_bases_between(idx, g_pos, first_cds.g_5prime) - 1  # exclusive of start codon
        return PosAnno(g_pos=g_pos, domain="UTR5", c_coord=f"-{n}")
    elif d_to_stop > 0:
        # 3' UTR
        n = _count_tx_bases_between(idx, last_cds.g_3prime, g_pos) - 1
        return PosAnno(g_pos=g_pos, domain="UTR3", c_coord=f"*{n}")
    else:
        # Shouldn't happen; fall through as CDS
        return PosAnno(g_pos=g_pos, domain="CDS", c_coord="?")


def _count_tx_bases_between(idx: TranscriptIndex, g_a: int, g_b: int) -> int:
    """
    Count the number of exonic (transcript) bases from g_a (inclusive) to g_b
    (inclusive), following the transcript in its 5'->3' direction. Assumes
    g_a is 5' of (or equal to) g_b in transcript order, and both lie within
    exons.
    """
    strand = idx.strand
    total = 0
    for ex in idx.exons_tx:
        lo, hi = min(ex.g_5prime, ex.g_3prime), max(ex.g_5prime, ex.g_3prime)
        seg_lo = max(lo, min(g_a, g_b))
        seg_hi = min(hi, max(g_a, g_b))
        if seg_lo > seg_hi:
            continue
        # Only count if this exon is between the two in tx order
        # (approximation: clamp to exon range)
        total += seg_hi - seg_lo + 1
    return total


# ---------------------------------------------------------------------------
# Edit application and per-guide consequence

@dataclass
class AnnotatedEdit:
    g_pos: int
    ref_plus: str              # + strand ref base
    alt_plus: str              # + strand alt base
    ref_tx: str                # transcript-sense ref (= comp if -strand)
    alt_tx: str                # transcript-sense alt
    domain: str
    c_coord: str | None
    hgvs_nt: str | None        # e.g., '1005+7A>G' or '128A>T' or '-25T>C'
    is_splice_donor: bool
    is_splice_acceptor: bool
    cds_pos: int | None
    codon_idx: int | None
    frame: int | None
    category: str              # per-edit tag used for Beagle-style Mutation Category string


def annotate_edits(
    idx: TranscriptIndex,
    edits: list[tuple[int, str, str]],
    cds_seq: str | None = None,
) -> tuple[list[AnnotatedEdit], dict]:
    """
    Given a transcript index and a list of (g_pos, ref_plus, alt_plus) edits,
    return per-edit annotations and a combined-codon summary.

    cds_seq is the transcript-sense CDS nucleotide sequence (ATG...stop). If
    provided, amino acid changes are computed for CDS edits. If None, AA
    changes are returned as '(NC)'.

    Returns:
        per_edit: list[AnnotatedEdit]
        summary : {
            'categories': list[str] (one per edit, in input order),
            'aa_edits': list[str]   (one per edit; same length as edits; '(NC)' for non-coding),
            'combined_aa_changes': list[str]  (one per distinct affected codon, combined),
            'mutation_category': str (comma-joined categories, Beagle-style),
        }
    """
    per_edit: list[AnnotatedEdit] = []
    cds_changes: dict[int, dict] = {}   # codon_idx -> {'frame_to_alt': {frame: alt_tx}, 'ref_codon': str, ...}

    for g_pos, ref_plus, alt_plus in edits:
        anno = annotate_position(idx, g_pos)
        ref_tx = ref_plus if idx.strand == 1 else comp(ref_plus)
        alt_tx = alt_plus if idx.strand == 1 else comp(alt_plus)

        hgvs_nt = None
        if anno.c_coord is not None and anno.c_coord != "?":
            hgvs_nt = f"{anno.c_coord}{ref_tx}>{alt_tx}"

        # Category logic
        if anno.domain == "CDS":
            category = "CDS"  # placeholder, set after AA translation
        elif anno.domain == "Intron":
            if anno.is_splice_donor:
                category = "Splice-donor"
            elif anno.is_splice_acceptor:
                category = "Splice-acceptor"
            else:
                category = "Intron"
        elif anno.domain in ("UTR5", "UTR3"):
            category = "UTR"
        else:
            category = "Outside"

        per_edit.append(AnnotatedEdit(
            g_pos=g_pos,
            ref_plus=ref_plus,
            alt_plus=alt_plus,
            ref_tx=ref_tx,
            alt_tx=alt_tx,
            domain=anno.domain,
            c_coord=anno.c_coord,
            hgvs_nt=hgvs_nt,
            is_splice_donor=anno.is_splice_donor,
            is_splice_acceptor=anno.is_splice_acceptor,
            cds_pos=anno.cds_pos,
            codon_idx=anno.codon_idx,
            frame=anno.frame,
            category=category,
        ))

    # --- Combined codon translation for CDS edits ---
    aa_edits = ["(NC)"] * len(per_edit)
    combined_aa_changes: list[str] = []

    if cds_seq:
        # Group CDS edits by codon
        cds_edits_by_codon: dict[int, list[int]] = {}
        for i, e in enumerate(per_edit):
            if e.domain == "CDS" and e.codon_idx is not None:
                cds_edits_by_codon.setdefault(e.codon_idx, []).append(i)

        for codon_idx, edit_indices in sorted(cds_edits_by_codon.items()):
            start_c = (codon_idx - 1) * 3
            if start_c + 3 > len(cds_seq):
                for i in edit_indices:
                    aa_edits[i] = "(?)"
                continue
            ref_codon = cds_seq[start_c : start_c + 3]
            alt_codon = list(ref_codon)
            for i in edit_indices:
                e = per_edit[i]
                if e.frame is None or not (0 <= e.frame <= 2):
                    continue
                alt_codon[e.frame] = e.alt_tx
            alt_codon_s = "".join(alt_codon)
            ref_aa = CODON_TABLE.get(ref_codon, "?")
            alt_aa = CODON_TABLE.get(alt_codon_s, "?")

            aa_str = f"p.{AA_THREE.get(ref_aa, ref_aa)}{codon_idx}{AA_THREE.get(alt_aa, alt_aa)}"
            combined_aa_changes.append(aa_str)

            if ref_aa == alt_aa:
                cat = "Silent"
            elif alt_aa == "*":
                cat = "Nonsense"
            elif ref_aa == "*":
                cat = "Stop-loss"
            elif ref_aa == "M" and codon_idx == 1:
                cat = "Start-loss"
            else:
                cat = "Missense"
            for i in edit_indices:
                per_edit[i].category = cat
                aa_edits[i] = aa_str

    categories = [e.category for e in per_edit]
    mutation_category = ", ".join(categories) if categories else ""

    return per_edit, {
        "categories": categories,
        "aa_edits": aa_edits,
        "combined_aa_changes": combined_aa_changes,
        "mutation_category": mutation_category,
    }


def build_cds_sequence(bundle: dict, idx: TranscriptIndex) -> str:
    """Build the transcript-sense CDS nucleotide sequence from the bundle's genomic seq."""
    seq = bundle["seq"]
    seq_start = bundle["seq_start"]
    strand = idx.strand
    cds_nt = []
    for seg in idx.cds_segments:
        g_lo = min(seg.g_5prime, seg.g_3prime)
        g_hi = max(seg.g_5prime, seg.g_3prime)
        s = seq[g_lo - seq_start : g_hi - seq_start + 1]
        if strand == -1:
            s = revcomp(s)
        cds_nt.append(s)
    return "".join(cds_nt)
