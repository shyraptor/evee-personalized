# Changelog

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
