# Beagle+

A local base-editor guide designer that does what the Broad's Beagle does
**plus** enumerates every partial-edit outcome in the editing window.

Beagle reports one row per guide, showing what happens when *all* editable
bases in the window fire at once. If a guide has 4 A's in the window, Beagle
gives you one outcome (A→G×4 simultaneously). Beagle+ still emits **one row
per guide**, but collapses all 15 non-empty subsets of those 4 positions
into that single row — every combination's HGVS nucleotide change, combined
codon translation, and mutation category is present in the relevant cell,
separated by ` | ` between combinations (and `, ` within a combination when
multiple bases fire).

This matters when a single guide can theoretically produce many different
protein variants (missense, silent, bystander combinations) depending on
which bases the editor actually hits in a given cell.

### Cell format at a glance

For a guide with 3 editable A's at protospacer positions 3, 7, 10 (7 combos):

```
Edit Combination          : A_3 | A_7 | A_10 | A_3 + A_7 | A_3 + A_10 | A_7 + A_10 | A_3 + A_7 + A_10
Nucleotide Edits (global) : 5A>G | 9A>G | 12A>G | 5A>G, 9A>G | 5A>G, 12A>G | 9A>G, 12A>G | 5A>G, 9A>G, 12A>G
Mutation Category         : Missense | Missense | Silent  | Missense    | Missense     | Missense     | Missense
Num Edits in Combination  : 1, 1, 1, 2, 2, 2, 3
Total Combinations        : 7
Worst Mutation Category   : Missense
Unique Mutation Categories: Missense, Silent
```

## What it does

- Takes an Ensembl transcript ID (ENST…).
- Fetches exons, CDS, and genomic sequence from Ensembl REST.
- Scans both strands for NGG (SpCas9) or NGN (SpCas9-NG / SpG) PAMs.
- By default restricts guides to **exons ± 20 nt** into each intron (matches
  Beagle's behavior). Adjust via the "Intron flank" field on the form, or
  `intron_flank` in the JSON API. Set to 0 for exon-only; set very large to
  scan whole introns (warning: 50k+ guides on long genes).
- Enumerates every 2^N − 1 non-empty subset of editable bases in the
  user-specified window.
- For each outcome, computes:
  - Genomic ref>alt edits on the + strand,
  - Transcript-sense HGVS c. notation (including intronic `+k/-k` and UTR
    `-N/*N`),
  - Combined codon translation for edits in CDS,
  - Mutation category (Missense / Silent / Nonsense / Splice-donor /
    Splice-acceptor / Intron / UTR / Outside),
  - Constraint flags (poly-T, BsmBI, BbsI).

## Columns

The output table has the 24 columns Beagle emits, in the same order, **plus**
three extension columns at the end:

1. Input, 2. CRISPR Enzyme, 3. Edit Type, 4. Edit Window,
5. Target Taxon, 6. Target Assembly, 7. Target Genome Sequence,
8. Target Gene ID, 9. Target Gene Symbol, 10. Target Gene Strand,
11. Target Transcript ID, 12. Target Domain,
13. sgRNA Sequence, 14. sgRNA Context Sequence, 15. PAM Sequence,
16. sgRNA Sequence Start Pos. (global), 17. sgRNA Orientation,
18. Nucleotide Edits (global), 19. Guide Edits, 20. Nucleotide Edits,
21. Amino Acid Edits, 22. Mutation Category, 23. Constraint Violations,
24. Note,

**New:** 25. Edit Combination, 26. Num Edits in Combination,
27. Total Combinations for Guide,
28. Worst Mutation Category (most severe outcome across combinations —
    sortable; severity ranks Nonsense > Start-loss > Splice-donor >
    Splice-acceptor > Missense > Silent > UTR > Intron > Outside),
29. Unique Mutation Categories (deduplicated atomic categories across
    combinations, e.g. `Missense, Silent`).

## Setup

```bash
cd beagle_plus
python3 -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

## Run

```bash
python app.py
```

Then open <http://localhost:5000> in your browser. Paste an ENST ID
(e.g. `ENST00000349945.7`), pick ABE/CBE, NGG/NGN, an edit window
(default 3..10), and submit.

The first 500 rows render in the browser. Download the full CSV for
everything.

## Validate against a Beagle export

```bash
python validate_vs_beagle.py beagle-export.xlsx ENST00000349945.7 \
       --editor ABE --pam NGN --window 3..10
```

This checks:
1. Every unique guide in the Beagle export appears in Beagle+'s output.
2. The full-edit row for each guide matches Beagle's Nucleotide Edits
   (global).
3. How many partial-edit rows Beagle+ adds on top of Beagle.

## Known limitations (read this before trusting the output)

1. **No probability weighting.** Every partial outcome is reported as
   equally possible. In reality, BE-Hive predicts some combinations are
   10× more likely than others. If you need to prioritize, wire in
   `maxwshen/be_predict_bystander`. This is the single biggest gap.
2. **Editor-specific sequence biases are ignored.** CBE with rAPOBEC1
   prefers TC context; Beagle+ does not downweight non-TC C's. Treat all
   outcomes as "theoretically possible," not "likely."
3. **Single Ensembl release, single assembly.** The Ensembl REST API
   returns whatever release the server currently points to. For
   reproducibility, record the assembly from the output and pin the
   Ensembl release via their URL scheme (`https://e111.rest.ensembl.org`
   for release 111) if needed.
4. **ENST only.** If you provide a gene ID (ENSG), the app fails.
   Different transcripts of the same gene have different exon structures;
   the user has to choose which transcript.
5. **Splice sites are ±1/±2 only.** Near-splice edits (e.g., branch point
   or +3 to +6 positions) are classified as plain "Intron." Edit those
   thresholds in `annotator.py` if you need extended splice windows.
6. **Combinatorial blowup.** A guide with 8 editable bases in the window
   produces 255 partial outcomes. Across 2000 guides × 256 outcomes =
   512k rows. Use the "Max edits per outcome" field to cap at e.g. 3.
7. **Intron flank default is 20 nt.** If you're hunting branch-point
   mutations, intronic enhancer edits, or anything deeper than 20 nt into
   an intron, raise "Intron flank". If you want exon-only, set it to 0.
7. **Off-target scoring is not performed.** Use CRISPick, Cas-OFFinder,
   or CHOPCHOP in addition to this tool.

## Architecture

```
ensembl_client.py   -- Ensembl REST wrapper (lookup, sequence fetch)
guide_finder.py     -- PAM scanning, protospacer + context extraction
edit_enumerator.py  -- 2^N partial edit subset enumeration
annotator.py        -- Transcript index + HGVS + AA translation
beagle_core.py      -- Orchestration, Beagle column assembly
app.py              -- Flask routes (form + results + CSV download + JSON API)
templates/          -- HTML templates for index + results
validate_vs_beagle.py -- Cross-check output against a Beagle export
test_synthetic.py   -- End-to-end test on a synthetic transcript
```

## JSON API

POST to `/api/generate` with JSON body:

```json
{
  "ensembl_id": "ENST00000349945.7",
  "editor": "ABE",
  "pam": "NGN",
  "window": [3, 10],
  "flank": 50,
  "max_edits": 3
}
```

Returns `{ "meta": {...}, "columns": [...], "rows": [[...], ...] }`.

## License / credits

Independent implementation. Produces output columns modeled after the
Broad Institute's Beagle portal, which is **not** open source. Base-editing
outcome science is from Arbab et al. (*Cell*, 2020) and the BE-Hive
repositories (`maxwshen/be_predict_bystander` and `maxwshen/be_predict_efficiency`).
