"""
Regression tests for the column-alignment bugs.

Run:  python test_alignment.py

Bug 1 (category misalignment)
    `combined_aa_changes` is indexed per affected CODON (ascending codon order);
    `categories` is indexed per EDIT (caller input order). The old code zipped
    them by position, so on antisense guides -- where protospacer order is the
    reverse of codon order -- every category landed on the wrong amino acid.
    A Trp->Ter next to a silent change came out as "Nonsense; Nonsense".

Bug 2 (non-coding edits dropped)
    Every intronic/UTR outcome across all 2^N-1 combos was stored under the
    single key "(Non-coding edit)", so only the first survived; and any combo
    that also contained a coding edit never registered its non-coding edits at
    all. Splice-site and intronic consequences silently vanished from the CSV.

Bug 3 (edit columns not parallel)
    Guide Edits / Nucleotide Edits carried the whole combo's edits against each
    amino acid, so those cells had a different arity from Amino Acid Edits and
    could not be split and zipped.
"""

from annotator import build_transcript_index, annotate_edits, build_cds_sequence
from spectra_core import generate_rows, COLUMNS, SEP

FAILURES = []


def check(label, got, want):
    ok = got == want
    print(f"  {'PASS' if ok else 'FAIL'}  {label}")
    if not ok:
        print(f"        got  {got!r}\n        want {want!r}")
        FAILURES.append(label)


# --------------------------------------------------------------------------
# Bug 1: categories must be parallel to combined_aa_changes
# --------------------------------------------------------------------------
def test_category_alignment():
    print("Bug 1 -- category/AA alignment")
    cds = "ATGGCTACTCAAGCTGATTTGATGGAGTTGGACATGGCCATGGAACCAGACAGAAAATGGCTGTAA"
    bundle = {
        "strand": 1,
        "exons": [{"start": 1001, "end": 1000 + len(cds)}],
        "cds_start": 1001,
        "cds_end": 1000 + len(cds),
        "seq": cds,
        "seq_start": 1001,
    }
    idx = build_transcript_index(bundle)
    seq = build_cds_sequence(bundle, idx)
    check("codon 20 is TGG (Trp)", seq[57:60], "TGG")
    check("codon 21 is CTG (Leu)", seq[60:63], "CTG")

    # Edits handed over in codon-DESCENDING order, as an antisense guide does.
    #   c.63 G>A -> CTG->CTA  Leu21Leu  Silent
    #   c.60 G>A -> TGG->TGA  Trp20Ter  Nonsense
    _, s = annotate_edits(idx, [(1063, "G", "A"), (1060, "G", "A")], cds_seq=seq)

    check("combined_aa_changes", s["combined_aa_changes"], ["p.Trp20Ter", "p.Leu21Leu"])
    check("combined_categories", s["combined_categories"], ["Nonsense", "Silent"])
    check("lengths match", len(s["combined_categories"]), len(s["combined_aa_changes"]))
    # The per-edit list is still input-ordered -- that is fine, it is a different axis.
    check("per-edit categories stay input-ordered", s["categories"], ["Silent", "Nonsense"])

    # A silent change must never be reported as Nonsense.
    for aa, cat in zip(s["combined_aa_changes"], s["combined_categories"]):
        if aa.endswith("Ter"):
            check(f"{aa} is Nonsense", cat, "Nonsense")
        elif aa[2:5] == aa[-3:]:
            check(f"{aa} is Silent", cat, "Silent")


# --------------------------------------------------------------------------
# Bug 2 + 3: full pipeline on a two-exon transcript with an intron
# --------------------------------------------------------------------------
def make_two_exon_bundle():
    """
    + strand, 2 exons, CDS spans both, 100 nt intron.
      Exon1 1001..1060   Intron 1061..1160   Exon2 1161..1260
    Intron starts GT and ends AG so the splice sites are well formed.
    """
    exon1 = "ATGGCTACTCAAGCTGATTTGATGGAGTTGGACATGGCCATGGAACCAGACAGAAAATGG"   # 60 nt
    intron = "GT" + "C" * 96 + "AG"                                            # 100 nt
    exon2 = "CTGTGGCCAACCAGACAGAAAATGGCTGTGGCCAACCAGACAGAAAATGGCTGTGGCCAACCAGACAGAAAATGGCTGTGGCCAACCAGACATGGCTGTAA"
    seq = exon1 + intron + exon2
    return {
        "strand": 1,
        "exons": [{"start": 1001, "end": 1060}, {"start": 1161, "end": 1160 + len(exon2)}],
        "cds_start": 1001,
        "cds_end": 1160 + len(exon2),
        "seq": seq,
        "seq_start": 1001,
        "ensembl_id": "ENSTTEST",
        "species": "homo_sapiens",
        "assembly": "TEST",
        "chrom": "1",
        "gene_id": "ENSGTEST",
        "gene_symbol": "TEST",
    }


def test_row_columns_parallel():
    print("\nBugs 2 & 3 -- every edit column parallel, non-coding edits kept")
    bundle = make_two_exon_bundle()
    rows = list(generate_rows(bundle, "CBE", "NGN", (3, 10), intron_flank=20))
    check("pipeline produced rows", bool(rows), True)

    ci = {c: i for i, c in enumerate(COLUMNS)}
    parallel = ["Nucleotide Edits (global)", "Guide Edits", "Nucleotide Edits",
                "Amino Acid Edits", "Mutation Category",
                "Edit Combination", "Num Edits in Combination"]

    bad_arity = 0
    bad_total = 0
    bad_guide_base = 0
    noncoding_rows = 0
    for r in rows:
        lens = {c: len(r[ci[c]].split(SEP)) for c in parallel}
        if len(set(lens.values())) != 1:
            bad_arity += 1
        n = lens["Amino Acid Edits"]
        if r[ci["Total Combinations for Guide"]] != str(n):
            bad_total += 1
        # CBE: every Guide Edit must name a C on the protospacer
        for tok in r[ci["Guide Edits"]].split(SEP):
            for part in tok.split(", "):
                if part and not part.startswith("C_"):
                    bad_guide_base += 1
        if "Intron" in r[ci["Mutation Category"]] or "Splice" in r[ci["Mutation Category"]]:
            noncoding_rows += 1

    check("all edit columns have equal arity", bad_arity, 0)
    check("Total Combinations matches arity", bad_total, 0)
    check("CBE guide edits all named C_n", bad_guide_base, 0)
    check("intronic/splice outcomes survive", noncoding_rows > 0, True)

    # Distinct intronic edits must not collapse into one bucket.
    for r in rows:
        cats = r[ci["Mutation Category"]].split(SEP)
        aas = r[ci["Amino Acid Edits"]].split(SEP)
        nts = r[ci["Nucleotide Edits"]].split(SEP)
        nc = [nts[i] for i, c in enumerate(cats) if c in ("Intron", "Splice-donor", "Splice-acceptor")]
        if len(nc) > 1:
            check("distinct intronic edits kept separate", len(set(nc)), len(nc))
            break

    # Categories must agree with the amino acid strings they sit next to.
    mismatched = 0
    for r in rows:
        for aa, cat in zip(r[ci["Amino Acid Edits"]].split(SEP),
                           r[ci["Mutation Category"]].split(SEP)):
            if len(aa) > 6 and aa[:3].isalpha() and aa[-3:].isalpha():
                want = "Nonsense" if aa.endswith("Ter") else (
                    "Silent" if aa[:3] == aa[-3:] else None)
                if want and cat != want and not (aa[:3] == "Met" and cat == "Start-loss"):
                    mismatched += 1
    check("no category contradicts its amino acid change", mismatched, 0)


if __name__ == "__main__":
    test_category_alignment()
    test_row_columns_parallel()
    print()
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S): " + ", ".join(FAILURES))
        raise SystemExit(1)
    print("All alignment regression tests passed.")
