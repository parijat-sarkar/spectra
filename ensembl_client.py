"""
Ensembl REST client.

Fetches transcript metadata (exons, CDS, strand, assembly) and genomic
sequence given an ENST ID. All results are cached in memory for the
lifetime of the process.
"""

from __future__ import annotations

import os
import time
from functools import lru_cache
from typing import Any

import requests

# Override with ENSEMBL_REST_URL to pin a release (e.g. https://e111.rest.ensembl.org)
# or to switch to a different mirror when the main endpoint is unhealthy.
ENSEMBL_REST = os.environ.get("ENSEMBL_REST_URL", "https://rest.ensembl.org").rstrip("/")

# Per-request socket timeout. Short on purpose: a hung Ensembl connection used
# to block a worker for a full minute at a time.
TIMEOUT = float(os.environ.get("ENSEMBL_TIMEOUT", 20))
MAX_RETRIES = int(os.environ.get("ENSEMBL_MAX_RETRIES", 5))
RETRY_DELAY = 1.0  # seconds, doubled each attempt
MAX_RETRY_DELAY = 8.0

# Hard ceiling on the time one logical fetch may spend retrying. This must stay
# comfortably below gunicorn's --timeout so we raise a catchable EnsemblError
# (which the app renders as a readable message) instead of letting the arbiter
# SIGABRT the worker -- SystemExit is a BaseException and escapes the app's
# `except Exception`, which is what produced bare 500s with no error page.
DEADLINE = float(os.environ.get("ENSEMBL_DEADLINE", 90))


class EnsemblError(RuntimeError):
    pass


def _get(path: str, params: dict | None = None, accept: str = "application/json") -> Any:
    """HTTP GET with exponential-backoff retry on 429 / 5xx, bounded by DEADLINE."""
    url = f"{ENSEMBL_REST}{path}"
    headers = {"Accept": accept}
    last_err = None
    started = time.monotonic()
    delay = RETRY_DELAY

    for attempt in range(MAX_RETRIES):
        remaining = DEADLINE - (time.monotonic() - started)
        if remaining <= 0:
            break
        try:
            r = requests.get(url, params=params or {}, headers=headers,
                             timeout=min(TIMEOUT, remaining))
            if r.status_code == 429 or r.status_code >= 500:
                last_err = EnsemblError(f"HTTP {r.status_code} for {url}")
                wait = float(r.headers.get("Retry-After", delay))
            elif not r.ok:
                # 4xx other than 429: a bad ID won't fix itself, fail immediately.
                raise EnsemblError(f"HTTP {r.status_code} for {url}: {r.text[:200]}")
            else:
                return r.json() if accept == "application/json" else r.text
        except requests.RequestException as e:
            last_err = e
            wait = delay

        if attempt < MAX_RETRIES - 1:
            time.sleep(min(wait, MAX_RETRY_DELAY, max(0.0, DEADLINE - (time.monotonic() - started))))
            delay = min(delay * 2, MAX_RETRY_DELAY)

    waited = time.monotonic() - started
    raise EnsemblError(
        f"Ensembl did not respond successfully after {MAX_RETRIES} attempts "
        f"over {waited:.0f}s ({last_err}). The REST service is likely having "
        f"a wobble -- wait a minute and retry, or set ENSEMBL_REST_URL to a "
        f"release-pinned mirror such as https://e111.rest.ensembl.org"
    )


@lru_cache(maxsize=256)
def get_transcript(ensembl_id: str) -> dict:
    """
    Fetch transcript lookup with expanded exon + Translation info.

    Accepts ENSTxxx or ENSTxxx.N form; the version suffix is stripped
    before querying, but the user's version is not enforced. Consumers
    should record assembly + release from the response.
    """
    base_id = ensembl_id.split(".", 1)[0]
    data = _get(f"/lookup/id/{base_id}", params={"expand": 1})
    if data.get("object_type") != "Transcript":
        raise EnsemblError(f"{ensembl_id} is not a transcript (got {data.get('object_type')})")
    return data


@lru_cache(maxsize=256)
def get_sequence(species: str, region: str) -> str:
    """
    Fetch a genomic region as an upper-case DNA string on the + strand.

    region format: "chrom:start..end" (1-based, inclusive). Always returns
    the + strand.
    """
    txt = _get(
        f"/sequence/region/{species}/{region}:1",
        accept="text/plain",
    )
    return txt.strip().upper()


def get_transcript_bundle(ensembl_id: str, flank: int = 50) -> dict:
    """
    Return a dict with everything needed downstream:
        ensembl_id, version, gene_id, gene_symbol, assembly, species,
        chrom, strand (+1/-1), tx_start, tx_end, exons (list of dicts,
        sorted by genomic start), cds_start, cds_end, translation_id,
        seq_start, seq_end, seq (plus-strand sequence from seq_start..seq_end)

    Coordinates are all 1-based inclusive. Exon list is ordered by
    genomic start (NOT by transcript order). For strand=-1 transcripts,
    transcript order is the reverse of the exon list.
    """
    t = get_transcript(ensembl_id)

    gene_id = None
    gene_symbol = None
    # Ensembl returns Parent as gene ID; gene symbol needs a second lookup.
    if "Parent" in t:
        gene_id = t["Parent"]
        try:
            g = get_transcript(gene_id) if gene_id.startswith("ENSG") else None
        except EnsemblError:
            g = None
        # Separate lookup for gene symbol
        try:
            g_data = _get(f"/lookup/id/{gene_id}")
            gene_symbol = g_data.get("display_name")
        except EnsemblError:
            pass

    exons = sorted(
        [
            {
                "id": e["id"],
                "start": e["start"],
                "end": e["end"],
                "strand": e.get("strand", t["strand"]),
            }
            for e in t.get("Exon", [])
        ],
        key=lambda e: e["start"],
    )
    cds_start = t["Translation"]["start"] if t.get("Translation") else None
    cds_end = t["Translation"]["end"] if t.get("Translation") else None
    translation_id = t["Translation"]["id"] if t.get("Translation") else None

    tx_start = t["start"]
    tx_end = t["end"]
    chrom = t["seq_region_name"]
    strand = t["strand"]  # +1 or -1
    assembly = t.get("assembly_name", "GRCh38")
    species = t.get("species", "human")

    # Translate Ensembl species to the REST species slug
    species_slug = {"human": "homo_sapiens"}.get(species, species)

    seq_start = max(1, tx_start - flank)
    seq_end = tx_end + flank
    region = f"{chrom}:{seq_start}..{seq_end}"
    seq = get_sequence(species_slug, region)

    return {
        "ensembl_id": ensembl_id,
        "version": t.get("version"),
        "gene_id": gene_id,
        "gene_symbol": gene_symbol,
        "assembly": assembly,
        "species": species,
        "species_slug": species_slug,
        "chrom": chrom,
        "strand": strand,
        "tx_start": tx_start,
        "tx_end": tx_end,
        "exons": exons,
        "cds_start": cds_start,
        "cds_end": cds_end,
        "translation_id": translation_id,
        "seq_start": seq_start,
        "seq_end": seq_end,
        "seq": seq,
    }
