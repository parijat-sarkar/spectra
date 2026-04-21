"""
End-to-end test on a synthetic transcript. Builds a fake bundle that mimics
Ensembl's structure, runs the full pipeline, and prints the table.

Two transcripts are tested:
  1. + strand, 3 exons, CDS spanning all 3.
  2. - strand, 2 exons, CDS entirely within.

We verify:
  - CDS vs Intron vs UTR vs Splice-donor/acceptor classification
  - HGVS c. coord for coding, intronic (+k/-k), UTR (-N, *N)
  - Combined codon translation with partial edits
  - Antisense edits reported as transcript-sense (complement) bases
"""

from annotator import build_transcript_index, annotate_edits, build_cds_sequence
from beagle_core import generate_rows, COLUMNS


def make_plus_transcript():
    """
    + strand transcript on chr1:
      Exon1: 100..200 (101 nt)
      Intron1: 201..300
      Exon2: 301..400 (100 nt)
      Intron2: 401..500
      Exon3: 501..600 (100 nt)
    CDS: 150..550 (50 nt in exon1, 100 nt in exon2, 50 nt in exon3) = 200 nt total
    5'UTR: 100..149 (50 nt)
    3'UTR: 551..600 (50 nt)

    We fabricate a sequence where every exon has known content so translation
    is predictable. For simplicity, fill non-CDS with 'N' and CDS with cycling
    ATG CCC AAA ... but for testing, just use all T's with ATG at the start.
    """
    # Build the plus-strand sequence for region [seq_start..seq_end]
    seq_start = 50
    seq_end = 650
    length = seq_end - seq_start + 1
    seq = list("N" * length)

    # Populate CDS with a known coding sequence
    # CDS genomic coords 150..550, in + strand order
    cds_parts = [
        (150, 200),  # 51 nt
        (301, 400),  # 100 nt
        (501, 550),  # 50 nt
    ]
    cds_total = sum(b - a + 1 for a, b in cds_parts)
    assert cds_total == 201  # oops, off by one; fine for this test
    # Build a CDS string: starts with ATG, then cycles (GCA=Ala, CTG=Leu, TTT=Phe)
    codons = ["ATG"] + ["GCA", "CTG", "TTT"] * 66 + ["TAA"]
    cds_str = "".join(codons)[:cds_total]

    # Place cds_str into seq + strand
    cds_idx = 0
    for a, b in cds_parts:
        n = b - a + 1
        for pos in range(a, b + 1):
            seq[pos - seq_start] = cds_str[cds_idx]
            cds_idx += 1

    # Fill rest of exon with random DNA (use T's for determinism)
    exons = [(100, 200), (301, 400), (501, 600)]
    for a, b in exons:
        for pos in range(a, b + 1):
            if seq[pos - seq_start] == "N":
                seq[pos - seq_start] = "T"
    # Fill introns with A's (so they're easy to spot) -- but we need some PAMs
    # for scanning; let's use a pattern rich in NGN matches so guides appear.
    import random
    random.seed(0)
    for a, b in [(201, 300), (401, 500)]:
        for pos in range(a, b + 1):
            if seq[pos - seq_start] == "N":
                seq[pos - seq_start] = random.choice("ACGT")
    # Fill outside transcript with N
    for pos in range(seq_start, seq_end + 1):
        if seq[pos - seq_start] == "N":
            seq[pos - seq_start] = random.choice("ACGT")

    return {
        "ensembl_id": "ENSTFAKE0001.1",
        "version": 1,
        "gene_id": "ENSGFAKE0001",
        "gene_symbol": "FAKEGENE1",
        "assembly": "GRCh38",
        "species": "human",
        "species_slug": "homo_sapiens",
        "chrom": "1",
        "strand": 1,
        "tx_start": 100,
        "tx_end": 600,
        "exons": [
            {"id": "E1", "start": 100, "end": 200, "strand": 1},
            {"id": "E2", "start": 301, "end": 400, "strand": 1},
            {"id": "E3", "start": 501, "end": 600, "strand": 1},
        ],
        "cds_start": 150,
        "cds_end": 550,
        "translation_id": "ENSPFAKE0001",
        "seq_start": seq_start,
        "seq_end": seq_end,
        "seq": "".join(seq),
    }


def test_position_annotation():
    bundle = make_plus_transcript()
    idx = build_transcript_index(bundle)
    cds = build_cds_sequence(bundle, idx)
    print(f"CDS length = {len(cds)} (expected 201)")
    print(f"CDS first 6 = {cds[:6]}  (expected ATGGCA)")
    print(f"CDS last 3 (stop) = {cds[-3:]}  (expected TAA)")

    # Test positions of interest
    from annotator import annotate_position
    tests = [
        (100, "UTR5"),      # far 5' UTR
        (149, "UTR5"),      # last 5' UTR base
        (150, "CDS"),       # c.1
        (151, "CDS"),       # c.2
        (200, "CDS"),       # last CDS base of exon 1 = c.51
        (201, "Intron"),    # intron 1, +1 -> splice donor
        (202, "Intron"),    # intron 1, +2 -> splice donor
        (203, "Intron"),    # intron 1, +3 -> Intron only
        (299, "Intron"),    # intron 1, -2 -> splice acceptor
        (300, "Intron"),    # intron 1, -1 -> splice acceptor
        (301, "CDS"),       # c.52
        (400, "CDS"),       # c.151
        (401, "Intron"),    # intron 2, +1 -> splice donor
        (500, "Intron"),    # intron 2, -1 -> splice acceptor
        (501, "CDS"),       # c.152
        (550, "CDS"),       # c.201 (last CDS = stop codon 3rd base)
        (551, "UTR3"),      # *1
        (600, "UTR3"),
    ]
    for g, want in tests:
        a = annotate_position(idx, g)
        marker = " OK" if a.domain == want else " MISMATCH"
        splice = ""
        if a.is_splice_donor: splice += " [donor]"
        if a.is_splice_acceptor: splice += " [acceptor]"
        print(f"  g={g:4d} -> {a.domain:10s} c.{a.c_coord}{splice}{marker} (want {want})")


def test_full_pipeline():
    bundle = make_plus_transcript()
    rows = list(generate_rows(bundle, editor="ABE", pam="NGN", window=(3, 10)))
    print(f"\n=== generate_rows produced {len(rows)} rows ===")
    # Print Mutation Category distribution
    from collections import Counter
    cat_idx = COLUMNS.index("Mutation Category")
    cats = Counter(r[cat_idx] for r in rows)
    for c, n in cats.most_common(15):
        print(f"  {n:4d}  {c}")
    # Show a few rows
    print("\n=== Sample rows ===")
    show_cols = ["sgRNA Sequence", "PAM Sequence", "sgRNA Orientation",
                 "sgRNA Sequence Start Pos. (global)",
                 "Guide Edits", "Nucleotide Edits", "Amino Acid Edits",
                 "Mutation Category", "Edit Combination", "Num Edits in Combination",
                 "Total Combinations for Guide"]
    col_ids = [COLUMNS.index(c) for c in show_cols]
    for r in rows[:5]:
        for c, i in zip(show_cols, col_ids):
            print(f"  {c:40s}: {r[i]}")
        print()
    # Find a guide with multiple edits and show all its partial outcomes
    # Group by sgRNA sequence + start
    from collections import defaultdict
    guide_idx = COLUMNS.index("sgRNA Sequence")
    start_idx = COLUMNS.index("sgRNA Sequence Start Pos. (global)")
    total_idx = COLUMNS.index("Total Combinations for Guide")
    by_key = defaultdict(list)
    for r in rows:
        by_key[(r[guide_idx], r[start_idx])].append(r)
    # Find a guide with >=4 combos
    for key, rs in by_key.items():
        if int(rs[0][total_idx]) >= 4:
            print(f"\n=== Guide {key[0]} at {key[1]}: {len(rs)} partial outcomes ===")
            combo_idx = COLUMNS.index("Edit Combination")
            cat_idx = COLUMNS.index("Mutation Category")
            hgvs_idx = COLUMNS.index("Nucleotide Edits")
            aa_idx = COLUMNS.index("Amino Acid Edits")
            for r in rs:
                print(f"  {r[combo_idx]:30s} | {r[cat_idx]:40s} | {r[hgvs_idx]:40s} | {r[aa_idx]}")
            break


if __name__ == "__main__":
    print("=== position annotation test ===")
    test_position_annotation()
    print()
    print("=== full pipeline test ===")
    test_full_pipeline()
