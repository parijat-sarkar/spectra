"""
Top-level pipeline: ENST in, full guide design row table out (with all
partial-edit combinations enumerated).

Public entry point:
    generate_table(ensembl_id, editor, pam, window, flank=50, max_edits=None)

Returns a dict with:
    columns : list of column names (24 standard columns)
    rows    : list of list[str] of table rows
    meta    : {ensembl_id, gene_symbol, gene_id, assembly, strand,
               n_guides, n_rows, chromosome, species}

Notes:
    - "Edit Window" input uses 1-based 1..20 scheme
      (e.g. (3, 10) covers positions 3 through 10).
    - PAM is one of: NGG, NGN, NG (NG will be treated as NG-N, so NGN).
    - Editor is "ABE" (A->G) or "CBE" (C->T).
    - Enumerates all 2^N-1 partial-edit combinations and collapses to unique AAs.
"""

from __future__ import annotations

from typing import Iterable

from ensembl_client import get_transcript_bundle
from guide_finder import find_guides, protospacer_pos_to_plus_coord
from edit_enumerator import enumerate_partial_outcomes, editable_positions, EDITOR_SPECS
from annotator import build_transcript_index, annotate_edits, build_cds_sequence


ENZYME_LABEL = {
    "NGG": "SpyoCas9",
    "NGN": "SpyoCas9NG",
    "NG": "SpyoCas9NG",
}

COLUMNS = [
    "Input",
    "CRISPR Enzyme",
    "Edit Type",
    "Edit Window",
    "Target Taxon",
    "Target Assembly",
    "Target Genome Sequence",
    "Target Gene ID",
    "Target Gene Symbol",
    "Target Gene Strand",
    "Target Transcript ID",
    "Target Domain",
    "sgRNA Sequence",
    "sgRNA Context Sequence",
    "PAM Sequence",
    "sgRNA Sequence Start Pos. (global)",
    "sgRNA Orientation",
    "Nucleotide Edits (global)",
    "Guide Edits",
    "Nucleotide Edits",
    "Amino Acid Edits",
    "Mutation Category",
    "Constraint Violations",
]


# Severity order for "Worst Mutation Category" (highest first wins).
_CATEGORY_SEVERITY = [
    "Nonsense",
    "Start-loss",
    "Splice-donor",
    "Splice-acceptor",
    "Missense",
    "Silent",
    "UTR",
    "Intron",
    "Outside",
    "",
]


def _worst_category(cats: list[str]) -> str:
    for want in _CATEGORY_SEVERITY:
        if want in cats:
            return want
    return cats[0] if cats else ""


def _flatten_categories(per_combo: list[str]) -> list[str]:
    """A per-combo mutation-category entry may itself be a comma-separated
    list (when a single combo spans multiple domains, e.g. 'Missense, Silent').
    Break those up so severity-ranking operates on atomic tokens."""
    flat: list[str] = []
    for entry in per_combo:
        for tok in (t.strip() for t in entry.split(",")):
            if tok:
                flat.append(tok)
    return flat


def _constraint_violations(sequence: str) -> str:
    """Flag Beagle-style constraints: poly-T runs >=4, strong GC extremes, BsmBI sites."""
    viols = []
    # Poly-T run of 4 or more (Pol III termination signal)
    import re
    m = re.search(r"T{4,}", sequence)
    if m:
        viols.append(f"poly(T):{m.group(0)}")
    # BsmBI site
    if "CGTCTC" in sequence or "GAGACG" in sequence:
        viols.append("BsmBI")
    # BbsI site
    if "GAAGAC" in sequence or "GTCTTC" in sequence:
        viols.append("BbsI")
    return ", ".join(viols) if viols else ""


def _pam_normalize(pam: str) -> str:
    pam = pam.upper()
    if pam == "NG":
        return "NGN"
    return pam


def _merge_intervals(raw: list[tuple[int, int]]) -> list[tuple[int, int]]:
    raw = sorted(raw, key=lambda x: x[0])
    merged: list[tuple[int, int]] = []
    for lo, hi in raw:
        if merged and lo <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], hi))
        else:
            merged.append((lo, hi))
    return merged


def _build_allowed_regions(idx, exons: list[dict], intron_flank: int) -> list[tuple[int, int]]:
    """
    Allowed guide region = CDS ± intron_flank, merged.

    Using CDS (not full exons) automatically drops 5'UTR / 3'UTR from scope
    while still extending `intron_flank` nt into each flanking intron AND
    `intron_flank` nt into 5'UTR/3'UTR past the start/stop codons.

    If the transcript has no CDS (non-coding), fall back to exons ± flank.
    """
    if idx.cds_segments:
        raw = [
            (max(1, min(s.g_5prime, s.g_3prime) - intron_flank),
             max(s.g_5prime, s.g_3prime) + intron_flank)
            for s in idx.cds_segments
        ]
    else:
        raw = [(max(1, e["start"] - intron_flank), e["end"] + intron_flank) for e in exons]
    return _merge_intervals(raw)


def _guide_in_allowed_region(start_global: int, protospacer_len: int,
                             allowed: list[tuple[int, int]]) -> bool:
    """True if the protospacer [start_global, start_global + P - 1] lies
    entirely within one allowed region."""
    lo = start_global
    hi = start_global + protospacer_len - 1
    # binary search could help; linear is fine for dozens of regions
    for a, b in allowed:
        if a <= lo and hi <= b:
            return True
        if b < lo:
            continue
        if a > hi:
            break
    return False


def generate_rows(
    bundle: dict,
    editor: str,
    pam: str,
    window: tuple[int, int],
    max_edits: int | None = None,
    intron_flank: int = 20,
) -> Iterable[list[str]]:
    """Yield Beagle-style rows (one per guide, with all 2^N-1 combo outcomes collapsed to unique AAs).

    For each guide with N editable positions in the window:
    - Enumerate all 2^N-1 non-empty subsets (combos)
    - For each combo, annotate all edits simultaneously (codon interactions captured)
    - Collapse to unique amino acid outcomes (same AA from different combos = one entry)
    - One row per guide lists all unique AAs, their categories, and worst category overall

    Guides are restricted to the union of (exon ± intron_flank) regions.
    Set intron_flank=0 to scan exons only; set to a large number (e.g. 10000)
    to scan the whole transcript.
    """
    pam = _pam_normalize(pam)
    if editor not in EDITOR_SPECS:
        raise ValueError(f"Unknown editor: {editor}")
    spec = EDITOR_SPECS[editor]

    idx = build_transcript_index(bundle)
    cds_seq = build_cds_sequence(bundle, idx) if idx.cds_segments else None

    allowed = _build_allowed_regions(idx, bundle["exons"], intron_flank)

    all_guides = find_guides(
        plus_seq=bundle["seq"],
        seq_start=bundle["seq_start"],
        pam_pattern=pam,
        protospacer_len=20,
    )
    guides = [g for g in all_guides
              if _guide_in_allowed_region(g["start_global"], 20, allowed)]

    input_tag = bundle["ensembl_id"]
    enzyme = ENZYME_LABEL.get(pam, f"SpyoCas9({pam})")
    edit_type = spec["edit_type"]
    window_str = f"{window[0]}..{window[1]}"
    taxon = bundle["species"]
    if taxon == "homo_sapiens" or taxon == "human":
        taxon_tag = "9606"
    else:
        taxon_tag = taxon
    assembly = bundle["assembly"]
    chrom_acc = f"{bundle['chrom']}"  # raw chromosome name
    gene_id = bundle["gene_id"] or ""
    gene_symbol = bundle["gene_symbol"] or ""
    strand_sym = "+" if bundle["strand"] == 1 else "-"
    transcript_id = bundle["ensembl_id"]

    # One row per guide. Enumerate all 2^N-1 combo outcomes, collapse to unique AAs.
    # For each combo, annotate all edits simultaneously so codon interactions are captured.
    for g in guides:
        ep = editable_positions(g["sequence"], editor, window)
        if not ep:
            continue

        # Enumerate all 2^N-1 combos of editable positions
        combos = enumerate_partial_outcomes(
            g["sequence"], editor, window, max_edits=max_edits
        )
        if not combos:
            continue

        # Collect all unique AA outcomes across all combos
        # Structure: aa_change -> {"cat": str, "nuc_global": str, "guide_edits": str, "nuc_hgvs": str}
        unique_aas: dict[str, dict] = {}
        all_domains: set[str] = set()

        for combo in combos:
            # Build list of edits for this combo (all positions in combo)
            edits = []
            edit_info = []  # Track (protospacer_pos, g_coord, ref, alt) for each edit

            for p in combo:
                g_coord = protospacer_pos_to_plus_coord(
                    g["start_global"], g["orientation"], p, protospacer_len=20
                )
                seq_offset = g_coord - bundle["seq_start"]
                ref_plus = bundle["seq"][seq_offset]
                if g["orientation"] == "sense":
                    alt_plus = spec["product_base"]
                else:
                    _C = {"A":"T","T":"A","C":"G","G":"C","N":"N"}
                    alt_plus = _C[spec["product_base"]]
                edits.append((g_coord, ref_plus, alt_plus))
                edit_info.append((p, g_coord, ref_plus, alt_plus))

            # Annotate all edits in this combo simultaneously (codon-aware)
            per_edit, summary = annotate_edits(idx, edits, cds_seq=cds_seq)

            # Build column values for this combo
            nuc_global_parts = []
            guide_edit_parts = []
            nuc_hgvs_parts = []

            for i, (p, g_coord, ref, alt) in enumerate(edit_info):
                if i < len(per_edit):
                    e = per_edit[i]
                    all_domains.add(e.domain)

                    # Nucleotide Edits (global): genomic coordinates like 23532211T>C
                    nuc_global_parts.append(f"{g_coord}{ref}>{alt}")

                    # Guide Edits: which base in guide like A_5, C_7
                    guide_edit_parts.append(f"{ref}_{p}")

                    # Nucleotide Edits (HGVS): transcript coordinates like 3828A>G
                    # Use the c_coord but extract just the position part
                    c_pos = e.c_coord.split(':')[-1] if ':' in e.c_coord else e.c_coord
                    nuc_hgvs_parts.append(f"{c_pos}{ref}>{alt}")

            # Get AA outcomes
            aa_outcomes = summary.get("combined_aa_changes", summary.get("aa_edits", []))
            cats = summary.get("categories", [])

            has_aa_outcome = False
            for i, aa_str in enumerate(aa_outcomes):
                if aa_str and aa_str not in ("(NC)", "(?)"):
                    has_aa_outcome = True
                    cat = cats[i] if i < len(cats) else ""

                    # Remove "p." prefix if present
                    aa_clean = aa_str.replace("p.", "") if aa_str else ""

                    if aa_clean not in unique_aas:
                        unique_aas[aa_clean] = {
                            "cat": cat,
                            "nuc_global": "; ".join(nuc_global_parts),
                            "guide_edits": ", ".join(guide_edit_parts),
                            "nuc_hgvs": "; ".join(nuc_hgvs_parts)
                        }
                    else:
                        # Keep worst category for this AA, merge edit lists
                        existing_cat = unique_aas[aa_clean]["cat"]
                        for severity_cat in _CATEGORY_SEVERITY:
                            if severity_cat == existing_cat:
                                break
                            if severity_cat == cat:
                                unique_aas[aa_clean]["cat"] = cat
                                break

            # If no AA outcome but we have edits, record them
            if not has_aa_outcome and nuc_global_parts:
                key = "(Non-coding edit)"
                cat = summary.get("categories", [""])[0] if summary.get("categories") else ""
                if key not in unique_aas:
                    unique_aas[key] = {
                        "cat": cat,
                        "nuc_global": "; ".join(nuc_global_parts),
                        "guide_edits": ", ".join(guide_edit_parts),
                        "nuc_hgvs": "; ".join(nuc_hgvs_parts)
                    }

        if "CDS" in all_domains:
            target_domain = "CDS"
        elif "UTR5" in all_domains or "UTR3" in all_domains:
            target_domain = "UTR"
        elif "Intron" in all_domains:
            target_domain = "CDS"
        else:
            target_domain = "Outside"

        # Format output
        aa_list = sorted(unique_aas.keys())
        cat_list = [unique_aas[aa]["cat"] for aa in aa_list]
        nuc_global_list = [unique_aas[aa].get("nuc_global", "") for aa in aa_list]
        guide_edits_list = [unique_aas[aa].get("guide_edits", "") for aa in aa_list]
        nuc_hgvs_list = [unique_aas[aa].get("nuc_hgvs", "") for aa in aa_list]
        worst_cat = _worst_category(cat_list)

        row = [
            input_tag,
            enzyme,
            edit_type,
            window_str,
            taxon_tag,
            assembly,
            chrom_acc,
            gene_id,
            gene_symbol,
            strand_sym,
            transcript_id,
            target_domain,
            g["sequence"],
            g["context"],
            g["pam"],
            str(g["start_global"]),
            g["orientation"],
            "; ".join(nuc_global_list),  # Nucleotide Edits (global) - genomic coords
            "; ".join(guide_edits_list),  # Guide Edits - which bases in guide
            "; ".join(nuc_hgvs_list),  # Nucleotide Edits (HGVS) - transcript coords
            "; ".join(aa_list),  # Amino Acid Edits (all unique AAs from all combos)
            "; ".join(cat_list),  # Mutation Categories (per AA)
            _constraint_violations(g["sequence"]),
        ]
        yield row


def generate_table(
    ensembl_id: str,
    editor: str,
    pam: str,
    window: tuple[int, int] = (3, 10),
    flank: int = 50,
    max_edits: int | None = None,
    intron_flank: int = 20,
) -> dict:
    bundle = get_transcript_bundle(ensembl_id, flank=max(flank, intron_flank + 50))
    rows = list(generate_rows(bundle, editor, pam, window,
                              max_edits=max_edits, intron_flank=intron_flank))
    meta = {
        "ensembl_id": ensembl_id,
        "gene_symbol": bundle.get("gene_symbol"),
        "gene_id": bundle.get("gene_id"),
        "assembly": bundle.get("assembly"),
        "strand": "+" if bundle["strand"] == 1 else "-",
        "chromosome": bundle["chrom"],
        "species": bundle["species"],
        "n_rows": len(rows),
    }
    return {"columns": COLUMNS, "rows": rows, "meta": meta}
