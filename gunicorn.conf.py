"""Gunicorn settings.

Gunicorn loads ./gunicorn.conf.py automatically, so these apply even when the
platform's start command is a bare `gunicorn app:app` with no flags.

Why the long timeout: a single request can spend minutes fetching a transcript
and its genomic sequence from Ensembl REST. Gunicorn's default 30s timeout
SIGABRTs the worker mid-request, which raises SystemExit -- a BaseException
that escapes the app's `except Exception`, so the user sees a bare 500 with no
explanation instead of a readable error page.

Keep this comfortably above ENSEMBL_DEADLINE in ensembl_client.py (90s) so the
client always gets to raise a catchable EnsemblError first.
"""

import os

bind = f"0.0.0.0:{os.environ.get('PORT', 10000)}"

# One worker: free tiers are memory-constrained and a whole gene's sequence is
# held in memory. Threads give concurrency without a second copy of the data.
workers = int(os.environ.get("WEB_CONCURRENCY", 1))
threads = int(os.environ.get("GUNICORN_THREADS", 4))
worker_class = "gthread"

timeout = int(os.environ.get("GUNICORN_TIMEOUT", 300))
graceful_timeout = 30
keepalive = 5

accesslog = "-"
errorlog = "-"
loglevel = "info"
