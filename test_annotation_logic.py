"""
Test annotation logic: verify that annotate_edits handles multi-edit combos correctly.
Simulates your TTT→multiple outcomes scenario.
"""

from annotator import annotate_edits, build_transcript_index, CODON_TABLE, AA_THREE

# Build a minimal mock transcript index
class MockCDSSegment:
    def __init__(self):
        self.g_5prime = 1
        self.g_3prime = 100
        self.c_start = 1
        self.c_end = 100
        self.exon_tx_idx = 0

class MockExonTx:
    def __init__(self):
        self.idx = 0
        self.g_5prime = 1
        self.g_3prime = 100
        self.strand = 1

class MockTranscriptIndex:
    def __init__(self):
        self.strand = 1
        self.exons_tx = [MockExonTx()]
        self.cds_segments = [MockCDSSegment()]
        self.total_cds_len = 100
        self.cds_start_g = 1
        self.cds_end_g = 100

# Simulate codon TTT (Phe) at genomic positions 10,11,12
# CDS position would be codon_idx = 4 (positions 10-12 map to codon 4)
idx = MockTranscriptIndex()

# Create a fake CDS sequence: ...XXXTTTYYYY...
# Position 10-12 are the TTT codon (codon_idx 4 in 0-based = positions 10-12 in 1-based)
# Actually, let's simplify: create a short CDS with TTT at the right place
# CDS: positions 1-3=codon1, 4-6=codon2, 7-9=codon3, 10-12=codon4 (TTT)
cds_seq = "ATGGGTAAATTTGGG"  # Codons: ATG(M), GGT(G), AAA(K), TTT(F), GGG(G)
#           1-3   4-6   7-9  10-12 13-15

print("CDS sequence:", cds_seq)
print("Codon 4 (positions 10-12): TTT → Phe\n")

# Test 1: Single edit at position 1 (T→C): CTT (Leu)
print("=" * 60)
print("TEST 1: Single edit pos 1 (T→C) → CTT = Leu")
edits = [(10, 'T', 'C')]  # genomic pos 10, ref T, alt C
per_edit, summary = annotate_edits(idx, edits, cds_seq=cds_seq)
print(f"AA changes: {summary['combined_aa_changes']}")
print(f"Categories: {summary['categories']}")
print()

# Test 2: Single edit at position 2 (T→C): TCT (Ser)
print("=" * 60)
print("TEST 2: Single edit pos 2 (T→C) → TCT = Ser")
edits = [(11, 'T', 'C')]  # genomic pos 11, ref T, alt C
per_edit, summary = annotate_edits(idx, edits, cds_seq=cds_seq)
print(f"AA changes: {summary['combined_aa_changes']}")
print(f"Categories: {summary['categories']}")
print()

# Test 3: Single edit at position 3 (T→C): TTC (Phe, SILENT!)
print("=" * 60)
print("TEST 3: Single edit pos 3 (T→C) → TTC = Phe (SILENT)")
edits = [(12, 'T', 'C')]  # genomic pos 12, ref T, alt C
per_edit, summary = annotate_edits(idx, edits, cds_seq=cds_seq)
print(f"AA changes: {summary['combined_aa_changes']}")
print(f"Categories: {summary['categories']}")
print()

# Test 4: Combo edits pos 1+2 (T→C both): CCT (Pro)
print("=" * 60)
print("TEST 4: Combo edits pos 1+2 (T→C both) → CCT = Pro")
edits = [(10, 'T', 'C'), (11, 'T', 'C')]  # both positions edited
per_edit, summary = annotate_edits(idx, edits, cds_seq=cds_seq)
print(f"AA changes: {summary['combined_aa_changes']}")
print(f"Categories: {summary['categories']}")
print()

# Test 5: Combo edits pos 1+2+3 (T→C all): TTT→CCC (Pro)
print("=" * 60)
print("TEST 5: Combo edits pos 1+2+3 (T→C all) → CCC = Pro")
edits = [(10, 'T', 'C'), (11, 'T', 'C'), (12, 'T', 'C')]  # all three
per_edit, summary = annotate_edits(idx, edits, cds_seq=cds_seq)
print(f"AA changes: {summary['combined_aa_changes']}")
print(f"Categories: {summary['categories']}")
print()

print("=" * 60)
print("✓ If you see all 5 outcomes above, annotation is working!")
print("✓ Notice: TTC is SILENT (Phe→Phe), others are MISSENSE except maybe Pro")
