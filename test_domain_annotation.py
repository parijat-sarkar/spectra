"""
Test domain + position annotations in output.
Verify that Note column now shows: "Splice-donor:+2; Arg294Arg" etc.
"""

from annotator import annotate_position, build_transcript_index

# Mock transcript with intron
class MockExonTx:
    def __init__(self, idx, g_5prime, g_3prime, strand):
        self.idx = idx
        self.g_5prime = g_5prime
        self.g_3prime = g_3prime
        self.strand = strand

class MockCDSSegment:
    def __init__(self, g_5prime, g_3prime, c_start, c_end):
        self.g_5prime = g_5prime
        self.g_3prime = g_3prime
        self.c_start = c_start
        self.c_end = c_end
        self.exon_tx_idx = 0

class MockTranscriptIndex:
    def __init__(self):
        # Two exons: Exon1 (1-100), Intron (101-200), Exon2 (201-300)
        self.strand = 1
        self.exons_tx = [
            MockExonTx(0, 1, 100, 1),
            MockExonTx(1, 201, 300, 1),
        ]
        # CDS spans both exons
        self.cds_segments = [
            MockCDSSegment(1, 100, 1, 100),
            MockCDSSegment(201, 300, 101, 200),
        ]
        self.total_cds_len = 200
        self.cds_start_g = 1
        self.cds_end_g = 300

idx = MockTranscriptIndex()

print("=" * 60)
print("DOMAIN ANNOTATION TEST")
print("=" * 60)
print()

# Test 1: CDS position
print("TEST 1: Position in CDS (should show 'CDS')")
anno = annotate_position(idx, 50)  # In first exon CDS
print(f"  Position 50 domain: {anno.domain}")
print(f"  C-coordinate: {anno.c_coord}")
print()

# Test 2: Intron position (donor)
print("TEST 2: Splice donor position (should show 'Splice-donor:+2')")
anno = annotate_position(idx, 102)  # +2 from exon end (position 100+2)
print(f"  Position 102 domain: {anno.domain}")
print(f"  Is splice donor: {anno.is_splice_donor}")
print(f"  C-coordinate: {anno.c_coord}")
print()

# Test 3: Intron position (acceptor)
print("TEST 3: Splice acceptor position (should show 'Splice-acceptor:-2')")
anno = annotate_position(idx, 199)  # -2 from next exon (position 201-2)
print(f"  Position 199 domain: {anno.domain}")
print(f"  Is splice acceptor: {anno.is_splice_acceptor}")
print(f"  C-coordinate: {anno.c_coord}")
print()

print("=" * 60)
print("✓ If you see splice-donor and splice-acceptor flags above,")
print("  domain annotation is working correctly.")
print("=" * 60)
