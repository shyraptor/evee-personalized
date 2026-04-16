# evee-personalized

Fetch [EVEE](https://evee.goodfire.ai) variant-effect interpretations for the
variants in a consumer DNA raw-data export (MyHeritage, 23andMe, AncestryDNA,
FamilyTreeDNA, Living DNA, …), pre-filtered via ClinVar. Produces a JSON +
Markdown report of the variants you actually carry, with a mechanistic,
LLM-generated interpretation for each.

**This is not medical advice.** See [Caveats](#caveats) below.

## Built on top of EVEE

All variant interpretations in this tool come from **EVEE (Evo Variant Effect
Explorer)**, a research project from [Goodfire AI](https://goodfire.ai) that
combines embeddings from the Evo 2 genomic foundation model with a
Claude-Sonnet interpretation layer:

- **Web tool** — <https://evee.goodfire.ai>
- **Preprint** — Goodfire AI (2026), *EVEE: Interpretable variant effect
  prediction from genomic foundation model embeddings*, bioRxiv.
  <https://www.biorxiv.org/content/10.64898/2026.04.10.717844v1>
  (use the "Cite this article" button on that page for a formatted
  citation / BibTeX entry).
- **Manuscript / figure-reproduction repo** —
  <https://github.com/goodfire-ai/evee-manuscript>

This repository is an independent client that queries EVEE's public API. It
is not affiliated with or endorsed by Goodfire AI. If you publish anything
based on results from this tool, please cite the EVEE paper above.

## What it does

1. **build** — parses your raw-data CSV, intersects by rsid against
   NCBI's [ClinVar `variant_summary.txt.gz`](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz)
   (GRCh38 forward strand, classified SNVs), and writes a `candidates.json`
   of unique (rsid × allele) records.
2. **fetch** — queries EVEE's `/variants/{id}/analysis` endpoint for each
   candidate, polling EVEE's async generation queue as needed, and appends
   each result to a resumable `evee_cache.ndjson`.
3. **emit** — applies filters and produces `evee_results.json` (an
   AI-chat-ready array) and `evee_results.md` (a human-readable report).
4. **all** — runs all three in order.

By default, `fetch` and `emit` keep only:

- **variants you actually carry** — your genotype contains at least one alt
  allele — and
- **non-benign ClinVar classifications** — keeps P/LP, VUS, conflicting,
  drug response, risk factor, association; drops benign / likely benign /
  benign+likely benign.

Widen either axis with `--no-carriers-only` and/or `--filter all`.

## Requirements

- Python 3.9 or newer.
- No third-party dependencies — the pipeline is 100% standard library.

## Inputs you need

### 1. Your consumer DNA raw-data file

The parser handles every major direct-to-consumer format:

| Vendor | Format | Notes |
|---|---|---|
| **MyHeritage** | `.csv`, comma-separated | `RSID / CHROMOSOME / POSITION / RESULT`, quoted values, `--` no-call |
| **23andMe** | `.txt`, tab-separated | Commented header (`# rsid	…`), `--` no-call. Also supports the `allele1`/`allele2` split-column variant. |
| **AncestryDNA** | `.txt`, tab-separated | Uncommented header, split `allele1`/`allele2` columns, `0` no-call |
| **FamilyTreeDNA** | `.csv`, comma-separated | Standard format (like MyHeritage) or "Family Finder" variant (commented header, `name` + split alleles) |
| **Living DNA** | `.csv` extension but tab-separated | Commented header, same shape as 23andMe |

The parser auto-detects delimiter (tab vs comma) from the header line and
ignores the file extension — all you need to do is point `--input` at the
file. You can download your raw data from your testing provider's results
page (typically labelled "Download DNA data" or similar).

The raw-data file itself is **never uploaded anywhere by this tool** — all
parsing and filtering runs locally. The only network traffic is EVEE
`variant_id` lookups of the form `chr:pos:ref:alt`. With the default
`--carriers-only`, the set of variants queried *is* a subset of your
genotype at ClinVar-classified loci: collectively, those queries encode
your pathogenic / likely-pathogenic / VUS carrier status at every queried
site. EVEE sees and caches each lookup server-side. This tool is not an
anonymizing proxy; treat it as a side-channel that reveals your genetic
makeup at the queried loci to Goodfire's backend.

If you have a file in some other layout, the parser will work as long as it
has columns matching one of the alias sets: `rsid / rs_id / rs / name`,
`chromosome / chrom / chr`, `position / pos / bp`, and either
`result / genotype / gt` *or* `allele1 / allele2`.

The `tests/fixtures/` directory contains a small synthetic file of every
supported format if you want to see the exact shapes.

### 2. ClinVar's `variant_summary.txt.gz`

Download the latest snapshot from NCBI (no login required; ~400–500 MB
compressed):

```bash
curl -O https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
```

## Usage

```bash
# Sanity-check that the EVEE API is reachable and our coordinate offset is
# still correct. No local files needed — it queries a single well-known
# variant and verifies EVEE returns a result at the expected coordinates.
python3 evee_pipeline.py verify

# One-shot, default filters (non_benign + carriers-only)
python3 evee_pipeline.py all

# Or step-by-step
python3 evee_pipeline.py build
python3 evee_pipeline.py fetch
python3 evee_pipeline.py emit

# Smoke test: first 20 filtered candidates
python3 evee_pipeline.py fetch --sample 20

# Include ClinVar-benign variants too
python3 evee_pipeline.py fetch --filter all
python3 evee_pipeline.py emit   --filter all

# Include non-carriers (full panel view)
python3 evee_pipeline.py fetch --no-carriers-only
python3 evee_pipeline.py emit   --no-carriers-only

# Only Pathogenic / Likely-pathogenic carriers
python3 evee_pipeline.py fetch --filter pathogenic
python3 evee_pipeline.py emit   --filter pathogenic
```

Full CLI: `python3 evee_pipeline.py --help`.

### Paths

All inputs and outputs are configurable via flags:

```
--input PATH         Consumer DNA raw-data file (default: MyHeritage_raw_dna_data.csv)
--clinvar PATH       ClinVar gz     (default: clinvar_variant_summary.txt.gz)
--candidates PATH    Build output / fetch input (default: candidates.json)
--cache PATH         Resumable NDJSON cache (default: evee_cache.ndjson)
--out-json PATH      Final JSON report (default: evee_results.json)
--out-md PATH        Final Markdown report (default: evee_results.md)
```

### Network

```
--api-base URL       EVEE API base (env: EVEE_API_BASE)
--concurrency N      Max concurrent workers (env: EVEE_CONCURRENCY, default 8)
--qps N              Soft req-per-second limit (env: EVEE_QPS, default 15)
```

## Output

`evee_results.json` is a JSON array, one record per qualifying variant,
sorted by (carrier status, pathogenicity rank, gene, rsid). Each record:

```json
{
  "rsid": "rs1800896",
  "variant_id": "chr1:206773551:T:C",
  "resolved_variant_id": "chr1:206773551:T:C",
  "found_via": "direct",
  "genotype": "TC",
  "clinvar": {
    "chrom": "1",
    "pos": 206773552,
    "ref": "T",
    "alt": "C",
    "gene": "IL10",
    "significance": "pathogenic",
    "phenotypes": "Leprosy, susceptibility to, 1|Hepatitis C virus, susceptibility to",
    "review_status": "criteria provided, single submitter",
    "name": "NC_000001.11:g.206773552T>C"
  },
  "evee_status": "ok",
  "evee": {
    "confidence": "high",
    "model": "claude-sonnet-4-6",
    "mechanism": "…",
    "summary": "…",
    "key_evidence": ["…", "…"],
    "variant_id": "chr1:206773551:T:C",
    "generated_at": 1776272912.95
  },
  "zygosity": "het",
  "user_carries_alt": true,
  "palindromic": false
}
```

`evee_results.md` is the same data formatted as a human-readable personal
variant report, grouped by classification and sorted carrier-first.

Both outputs are suitable for pasting into an AI chat (Claude, ChatGPT, …)
for further discussion, though the JSON is typically denser and easier to
cite.

## Caveats

- **Not medical advice.** EVEE is a research tool; its interpretations are
  LLM-generated (via Claude-Sonnet on top of Evo-2 embeddings) and can be
  wrong. ClinVar classifications change over time, and a single ClinVar
  "pathogenic" label — especially for susceptibility / risk-factor /
  drug-response entries — is *not* the same as a high-penetrance disease
  allele.
- **Palindromic SNVs (A/T, C/G) are strand-ambiguous** in principle. The
  pipeline classifies them by direct allele matching (trusting both
  sources' stated forward-strand convention) but flags them in the output
  so you can eyeball suspicious ones.
- **Indels are skipped.** Consumer DNA raw-data files report SNV calls
  only, so the pipeline filters ClinVar to SNV-compatible entries
  accordingly.
- **Coordinate offset.** ClinVar's `PositionVCF` is 1-based; EVEE's
  `variant_id` is 0-based. The pipeline subtracts 1 when constructing the id.
- **EVEE caches server-side.** Variants someone has already looked up come
  back instantly. Cold variants trigger a one-time Claude-Sonnet generation
  on Goodfire's backend — be reasonable about volume (default concurrency
  is 8, which is polite).
- **Low-pass WGS files** commonly contain internal IDs for positions
  without a dbSNP rsid. Those rows are silently dropped — the pipeline only
  keeps rsid-matched entries.

## Fair use of EVEE's API

This tool queries EVEE's public, unauthenticated API. There are no
published rate limits or terms of service, but please be respectful:

- **Cold variants cost real money.** Every variant that hasn't been looked
  up before triggers a one-time Claude-Sonnet generation on Goodfire's
  backend. The result is cached permanently so subsequent lookups are free,
  but the first hit is real compute on their infrastructure.
- **Keep concurrency reasonable.** The default (8 workers, 15 req/s cap) is
  polite. Don't crank `--concurrency` to 100.
- **Use `--carriers-only` (the default).** It reduces your query footprint
  by 5–10x compared to querying every ClinVar-matched variant regardless
  of whether you carry it.
- **The cache is your friend.** `evee_cache.ndjson` means re-runs are free.
  Don't delete it between runs unless you have a reason to.
- **Report issues upstream.** If the API is down, slow, or returns errors,
  that's a Goodfire infrastructure question — reach out to them via
  <https://goodfire.ai>, not via issues on this repo.

## How it handles EVEE's async API

EVEE's `/variants/{id}/analysis` endpoint is asynchronous. The first request
for a cold variant queues a Claude-Sonnet generation and returns
`{"status": "processing", "retry_after": 10}`. The fetcher polls on an
adaptive client-side schedule (3, 4, 5, 7, 10, 12, 15 seconds) rather than
blindly honoring `retry_after=10`, which is a conservative server-side
default — most variants actually complete within ~5–10 seconds. Successful
results are cached persistently server-side, and the local NDJSON cache
means re-runs are cheap.

On 404, the fetcher falls back to `/variants/search?q={rsid}` once, to
handle multi-allelic sites where ClinVar's alt allele happens to differ
from EVEE's canonical allele.

On 429 / 5xx / network errors, the fetcher uses exponential backoff up to 5
retries. A variant that stays queued past `max_poll_seconds` (default 300)
is recorded with `evee_status: "error:poll_timeout"` and can be retried
later — just delete its row from the cache and re-run `fetch`.

## Zygosity + carrier detection

Consumer DNA raw-data files from all major vendors (MyHeritage, 23andMe,
AncestryDNA, FTDNA, Living DNA) report genotypes on the **forward (+) strand
of GRCh37**. ClinVar's `PositionVCF` / `RefVCF` / `AltVCF` columns are on
the **forward strand of GRCh38**. The forward strand at a given physical
locus is the same strand in both assemblies, so the pipeline compares
alleles directly — no complement / strand-flip logic.

Non-matching genotypes are labeled `other:<gt>` and treated as non-carriers,
which is conservative against multi-allelic sites and sequencing errors.

## Tests

```bash
python3 tests/test_parser.py
```

Runs the parser against the synthetic fixtures in `tests/fixtures/` and
prints PASS/FAIL per vendor. Exit code is non-zero if anything fails.

## Acknowledgements

- **EVEE** (the variant-effect interpretation engine this tool wraps) is
  built by [Goodfire AI](https://goodfire.ai). See the
  [Built on top of EVEE](#built-on-top-of-evee) section above for the
  paper, repo, and web-tool links.
- **ClinVar** is maintained by the US National Center for Biotechnology
  Information (NCBI).
- **MyHeritage / 23andMe / AncestryDNA / FTDNA / Living DNA** are the
  consumer DNA testing services whose raw-data export formats this tool
  accepts as input.

## License

MIT — see [LICENSE](LICENSE).
