"""
Quick test: verify combo enumeration logic with your TTT→Phe example.
"""

from edit_enumerator import enumerate_partial_outcomes, editable_positions

# Example: CCC (Pro) with CBE editor (C→T)
# CCC has 3 editable C's → 2^3 - 1 = 7 combos

protospacer = "ATGCCCGGTAA"  # Contains CCC at positions 4,5,6
editor = "CBE"
window = (4, 6)  # 1-based, the CCC codon

ep = editable_positions(protospacer, editor, window)
print(f"Editable positions in window {window}: {ep}")
# Expected: [4, 5, 6] (all three C's)

combos = enumerate_partial_outcomes(protospacer, editor, window)
print(f"\nAll {len(combos)} combo outcomes:")
for i, combo in enumerate(combos, 1):
    print(f"  {i}. Positions {combo}")

# Manually verify:
# Pos 4,5,6 are the three C's
# 2^3 - 1 = 7 combos:
# {4} → ref_codon=CCC, alt at pos1 → TCC
# {5} → ref_codon=CCC, alt at pos2 → CTC
# {6} → ref_codon=CCC, alt at pos3 → CCT
# {4,5} → TTC
# {4,6} → TCT
# {5,6} → CTT
# {4,5,6} → TTT

print("\nExpected 7 combos total. Got:", len(combos))
print("✓ Combo enumeration works!" if len(combos) == 7 else "✗ Mismatch!")
