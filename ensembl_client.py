"""
Ensembl REST client.

Fetches transcript metadata (exons, CDS, strand, assembly) and genomic
sequence given an ENST ID. All results are cached in memory for the
lifetime of the process.
"""

from __future__ import annotations

import time
from functools import lru_cache
from typing import Any

import requests

ENSEMBL_REST = "https://rest.ensembl.org"
TIMEOUT = 60
MAX_RETRIES = 3
RETRY_DELAY = 2  # seconds


class EnsemblError(RuntimeError):
    pass


def _get(path: str, params: dict | None = None, accept: str = "application/json") -> Any:
    """HTTP GET with basic retry on 429 / 5xx."""
    url = f"{ENSEMBL_REST}{path}"
    headers = {"Accept": accept}
    last_err = None
    for attempt in range(MAX_RETRIES):
        try:
            r = requests.get(url, params=params or {}, headers=headers, timeout=TIMEOUT)
            if r.status_code == 429 or r.status_code >= 500:
                retry_after = float(r.headers.get("Retry-After", RETRY_DELAY))
                time.sleep(retry_after)
                last_err = EnsemblError(f"HTTP {r.status_code} for {url}")
                continue
            if not r.ok:
                raise EnsemblError(f"HTTP {r.status_code} for {url}: {r.text[:200]}")
            if accept == "application/json":
                return r.json()
            return r.text
        except requests.RequestException as e:
            last_err = e
            time.sleep(RETRY_DELAY)
    raise EnsemblError(f"Failed after {MAX_RETRIES} retries: {last_err}")


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
