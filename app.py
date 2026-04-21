"""
Flask web app: Spectra (base-editor guide designer with full partial-edit enumeration).

Run:
    pip install -r requirements.txt
    python app.py

Then browse to http://localhost:5000
"""

from __future__ import annotations

import csv
import io
import time
from threading import Lock

from flask import Flask, render_template, request, send_file, abort, jsonify

from spectra_core import generate_table, COLUMNS
from ensembl_client import EnsemblError

app = Flask(__name__)

# Simple in-memory cache of the last N runs so the user can download CSVs.
_CACHE: dict[str, dict] = {}
_CACHE_ORDER: list[str] = []
_CACHE_LOCK = Lock()
_CACHE_MAX = 8


def _cache_put(key: str, value: dict) -> None:
    with _CACHE_LOCK:
        _CACHE[key] = value
        _CACHE_ORDER.append(key)
        while len(_CACHE_ORDER) > _CACHE_MAX:
            old = _CACHE_ORDER.pop(0)
            _CACHE.pop(old, None)


@app.route("/")
def index():
    return render_template("index.html", columns=COLUMNS)


@app.route("/generate", methods=["POST"])
def generate():
    try:
        ensembl_id = request.form.get("ensembl_id", "").strip()
        editor = request.form.get("editor", "ABE").strip().upper()
        pam = request.form.get("pam", "NGG").strip().upper()
        window_start = int(request.form.get("window_start", 3))
        window_end = int(request.form.get("window_end", 10))
        max_edits_raw = request.form.get("max_edits", "").strip()
        max_edits = int(max_edits_raw) if max_edits_raw else None
        flank = int(request.form.get("flank", 50))
        intron_flank = int(request.form.get("intron_flank", 20))

        if not ensembl_id:
            return render_template("index.html", columns=COLUMNS, error="ENST ID is required."), 400
        if editor not in ("ABE", "CBE"):
            return render_template("index.html", columns=COLUMNS, error="Editor must be ABE or CBE."), 400
        if pam not in ("NGG", "NGN", "NG"):
            return render_template("index.html", columns=COLUMNS, error="PAM must be NGG, NGN or NG."), 400
        if not (1 <= window_start <= window_end <= 20):
            return render_template("index.html", columns=COLUMNS, error="Edit window must satisfy 1 <= start <= end <= 20."), 400

        t0 = time.time()
        result = generate_table(
            ensembl_id=ensembl_id,
            editor=editor,
            pam=pam,
            window=(window_start, window_end),
            flank=flank,
            max_edits=max_edits,
            intron_flank=intron_flank,
        )
        elapsed = time.time() - t0

        # Cache under a simple key
        key = f"{ensembl_id}_{editor}_{pam}_{window_start}-{window_end}_{max_edits_raw}_{int(time.time())}"
        _cache_put(key, result)

        # Paginate for display (first 500 rows in HTML; full set in CSV)
        preview_limit = 500
        preview = result["rows"][:preview_limit]
        return render_template(
            "results.html",
            columns=result["columns"],
            rows=preview,
            meta=result["meta"],
            total_rows=len(result["rows"]),
            preview_limit=preview_limit,
            cache_key=key,
            elapsed=f"{elapsed:.1f}",
            form={
                "ensembl_id": ensembl_id,
                "editor": editor,
                "pam": pam,
                "window": f"{window_start}..{window_end}",
                "max_edits": max_edits_raw or "all",
                "intron_flank": intron_flank,
            },
        )
    except EnsemblError as e:
        return render_template("index.html", columns=COLUMNS, error=f"Ensembl error: {e}"), 502
    except Exception as e:
        return render_template("index.html", columns=COLUMNS, error=f"Error: {type(e).__name__}: {e}"), 500


@app.route("/download/<cache_key>.csv")
def download_csv(cache_key: str):
    # Thread-safe: copy cache data while holding lock, release before writing
    with _CACHE_LOCK:
        result = _CACHE.get(cache_key)
        if result is None:
            abort(404)
        # Deep copy to avoid holding lock during file write
        result = {
            "columns": result["columns"][:],
            "rows": [row[:] for row in result["rows"]]
        }

    # Write file outside lock
    buf = io.StringIO()
    w = csv.writer(buf)
    w.writerow(result["columns"])
    for row in result["rows"]:
        w.writerow(row)
    buf.seek(0)
    return send_file(
        io.BytesIO(buf.getvalue().encode("utf-8")),
        mimetype="text/csv",
        as_attachment=True,
        download_name=f"{cache_key}.csv",
    )


@app.route("/api/generate", methods=["POST"])
def api_generate():
    """JSON API for programmatic use."""
    data = request.get_json(force=True, silent=True) or {}
    try:
        result = generate_table(
            ensembl_id=data["ensembl_id"],
            editor=data.get("editor", "ABE"),
            pam=data.get("pam", "NGG"),
            window=tuple(data.get("window", [3, 10])),
            flank=int(data.get("flank", 50)),
            max_edits=data.get("max_edits"),
            intron_flank=int(data.get("intron_flank", 20)),
        )
        return jsonify({"meta": result["meta"],
                        "columns": result["columns"],
                        "n_rows": len(result["rows"]),
                        "rows": result["rows"]})
    except EnsemblError as e:
        return jsonify({"error": f"Ensembl: {e}"}), 502
    except Exception as e:
        return jsonify({"error": f"{type(e).__name__}: {e}"}), 500


if __name__ == "__main__":
    app.run(host="127.0.0.1", port=5000, debug=False)
