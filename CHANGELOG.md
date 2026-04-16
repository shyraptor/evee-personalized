# Changelog

## 1.1.0 — 2026-04-16

### Added

- **`verify` subcommand** — queries a single hardcoded well-known variant
  (rs1800896, IL10 promoter) to sanity-check API connectivity and the
  ClinVar→EVEE coordinate offset. No local files required. Useful before
  a long fetch run, and as a canary for Goodfire moving the API base URL.
- **Core-logic test suite** — `tests/test_core.py` covers `variant_id_for`,
  `zygosity`, `is_palindromic`, `_primary_sig`, `_is_ready`, `_is_queued`,
  `_passes_filters`, `_sort_key`, and `PATHO_RANK`, in addition to the
  existing parser tests.
- **`pyproject.toml`** — ruff config (`E`, `W`, `F`, `B`, `UP`), targeting
  Python 3.9. CI now runs `ruff check .` alongside the tests.
- **`__version__` constant** — exposed on the module and used in the
  `User-Agent` header so bumping is a one-line change.

### Changed

- **Response size cap (`MAX_RESPONSE_BYTES = 10 MB`)** in `http_get_json()`
  — defends against a misbehaving or hostile API base URL returning
  unbounded data.
- **Cache corruption counter** in `read_cache()` — corrupt / incomplete
  NDJSON lines are counted and reported on stderr (behavior unchanged —
  they're still skipped and re-queried on the next fetch).
- **README privacy section tightened** — clarifies that the *set* of
  variant-id lookups collectively encodes the user's carrier status at
  queried loci; this tool is not an anonymizing proxy.

### CI

- Bumped `actions/checkout` v4 → v5 and `actions/setup-python` v5 → v6
  (Node.js 24 runtime, ahead of June 2026 deprecation deadline).

## 1.0.0 — 2026-04-15

Initial release.

### Features

- **Multi-vendor parser** — auto-detects delimiter (tab vs comma) and header
  style (commented / uncommented / split alleles). Tested against synthetic
  fixtures for MyHeritage, 23andMe (standard + allele-split + Windows CRLF),
  AncestryDNA (standard + mitochondrial), FamilyTreeDNA (standard + Family
  Finder), and Living DNA.
- **ClinVar intersection** — uses NCBI's `variant_summary.txt.gz` (GRCh38)
  to pre-filter by rsid and attach classification metadata. Applies a
  blacklist of unclassified entries rather than a brittle whitelist.
- **EVEE fetch with async polling** — handles EVEE's queue-based analysis
  endpoint (HTTP 202 / `"status": "processing"`), polling on an adaptive
  schedule (3–15 s). Falls back to rsid search on 404 for multi-allelic
  sites. Resumable via NDJSON cache.
- **Carrier-only default** — only fetches/emits variants the user actually
  carries (direct allele match, no strand-flip guessing). Palindromic
  SNVs (A/T, C/G) are flagged in the output.
- **Non-benign default** — drops ClinVar benign / likely benign /
  benign+likely benign by default. Overridable with `--filter all`.
- **Coordinate offset** — converts ClinVar's 1-based `PositionVCF` to
  EVEE's 0-based `variant_id` by subtracting 1.
- **JSON + Markdown output** — `evee_results.json` (AI-chat-ready array)
  and `evee_results.md` (human-readable report), sorted by pathogenicity
  rank then gene.
- **Regression test suite** — 9 synthetic fixtures covering all supported
  formats, runnable via `python3 tests/test_parser.py`.
- **Zero dependencies** — Python 3.9+ standard library only.
