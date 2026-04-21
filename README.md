# Spectra

Base editor guide designer that enumerates **all** partial-edit combinations.

## What it does

- Enumerates all 2^N-1 partial-edit outcomes for base editors (ABE/CBE)
- Shows domain + position annotations (Splice-donor, Intron, CDS, etc.)
- Collapses to unique amino acid outcomes
- Generates CSV export for further analysis

## Quick Start

### Local

```bash
pip install -r requirements.txt
python app.py
```

Then visit: http://localhost:5000

### Online

Visit: https://spectra-32j8.onrender.com

(You'll get the actual URL from Render after deployment completes)

## Input

- ENST ID (Ensembl transcript)
- Editor: ABE or CBE
- PAM: NGG or NGN
- Edit window (e.g., 3-10)

## Output

Interactive table + CSV download with:
- All unique AA outcomes across all partial combos
- Domain annotations (Splice-donor:+2, Intron:+6, etc.)
- Mutation categories per outcome

## Technologies

- Python 3.10+
- Flask
- Ensembl REST API
