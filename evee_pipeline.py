#!/usr/bin/env python3
"""
evee-personalized: fetch EVEE variant-effect interpretations for variants
in a consumer DNA raw-data export, intersected against ClinVar.

Supported input formats:
  - MyHeritage (comma-separated, RESULT column, `--` no-call)
  - 23andMe (tab-separated, commented header, single genotype or split alleles)
  - AncestryDNA (tab-separated, split alleles, `0` no-call)
  - FamilyTreeDNA (standard + Family Finder variant)
  - Living DNA (tab-separated despite .csv extension)
  - Anything else with rsid/chromosome/position plus either a single
    genotype column or allele1/allele2 split columns

Pipeline modes:
  build  — parse the raw-data CSV, intersect by rsid against ClinVar
           (GRCh38, classified SNVs), write candidates.json
  fetch  — query https://evee.goodfire.ai /variants/{id}/analysis for each
           candidate, polling EVEE's async generation queue as needed,
           append to a resumable NDJSON cache
  emit   — read the cache, apply filters, write evee_results.json and
           evee_results.md
  all    — run build, fetch, emit in that order

Defaults: --carriers-only is on (only variants you carry at least one alt
allele of) and --filter is non_benign (drops ClinVar benign/LB/B+LB). Widen
with --no-carriers-only and/or --filter all.

This is not medical advice. EVEE is a research tool; ClinVar classifications
can change; LLM interpretations can be wrong. Verify any actionable result
with a qualified clinician.

Standard library only. Python 3.9+.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import os
import re
import sys
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any


# ---------- Defaults ----------

DEFAULT_INPUT = Path("MyHeritage_raw_dna_data.csv")
DEFAULT_CLINVAR = Path("clinvar_variant_summary.txt.gz")
DEFAULT_CANDIDATES = Path("candidates.json")
DEFAULT_CACHE = Path("evee_cache.ndjson")
DEFAULT_OUT_JSON = Path("evee_results.json")
DEFAULT_OUT_MD = Path("evee_results.md")

# EVEE's API Gateway URL. Extracted from https://evee.goodfire.ai/config.js —
# if the tool stops working because Goodfire has moved the backend, re-fetch
# that file and look for API_BASE, or set EVEE_API_BASE / pass --api-base.
DEFAULT_API_BASE = "https://xix0d0o8le.execute-api.us-east-1.amazonaws.com"
DEFAULT_CONCURRENCY = 8
DEFAULT_QPS = 15.0
DEFAULT_MAX_POLL_SECONDS = 300.0

USER_AGENT = "evee-personalized/1.0"


# ---------- Consumer DNA raw-data CSV parser ----------

RS_RE = re.compile(r"^rs\d+$", re.IGNORECASE)


def parse_dna_csv(path: Path) -> dict[str, dict[str, str]]:
    """Parse a consumer DNA raw-data file into {rsid: {chrom, pos, genotype}}.

    Auto-detects:
      - Delimiter (tab vs comma) from the header line — the file extension
        is ignored because Living DNA ships tab-separated content with a
        `.csv` extension.
      - Header location — `#`-commented preamble is skipped, and a
        commented header (as used by 23andMe / Living DNA / FTDNA Family
        Finder) is recognised by looking for a `#` line that contains
        `rsid` (or `name`) and has at least four delimiter-separated tokens.

    Column matching is alias-tolerant:
      rsid  — rsid, rs_id, rs, name
      chrom — chromosome, chrom, chr
      pos   — position, pos, bp
      gt    — result, genotype, gt  (single column)  OR
              allele1 + allele2      (split columns, concatenated)

    Rows whose rsid doesn't match `rs\\d+` are silently dropped — consumer
    files commonly contain internal IDs for positions without a dbSNP rsid.
    """
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()

    header_idx = -1
    header_text = ""
    for i, line in enumerate(lines):
        stripped = line.strip()
        if not stripped:
            continue
        if not stripped.startswith("#"):
            header_idx = i
            header_text = stripped
            break
        candidate = stripped.lstrip("#").strip()
        low = candidate.lower()
        has_rsid_word = "rsid" in low or bool(re.search(r"\bname\b", low))
        token_count = candidate.count("\t") + candidate.count(",") + 1
        if has_rsid_word and token_count >= 4:
            header_idx = i
            header_text = candidate
            break

    if header_idx < 0:
        raise RuntimeError(f"{path}: no header found")

    delimiter = "\t" if header_text.count("\t") > header_text.count(",") else ","

    body: list[str] = [header_text]
    for line in lines[header_idx + 1:]:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        body.append(line)

    reader = csv.DictReader(io.StringIO("\n".join(body)), delimiter=delimiter)
    header_map = {(h or "").lower().strip(): h for h in (reader.fieldnames or [])}

    def col(*aliases: str) -> str | None:
        for alias in aliases:
            if alias in header_map:
                return header_map[alias]
        return None

    c_rs = col("rsid", "rs_id", "rs", "name")
    c_ch = col("chromosome", "chrom", "chr")
    c_po = col("position", "pos", "bp")
    c_gt = col("result", "genotype", "gt")
    c_a1 = col("allele1", "allele 1")
    c_a2 = col("allele2", "allele 2")

    if not (c_rs and c_ch and c_po):
        raise RuntimeError(
            f"{path}: missing required columns. "
            f"Need rsid (or name) / chromosome / position; got {reader.fieldnames}"
        )
    if not c_gt and not (c_a1 and c_a2):
        raise RuntimeError(
            f"{path}: missing genotype column. "
            f"Need result/genotype or allele1+allele2; got {reader.fieldnames}"
        )

    out: dict[str, dict[str, str]] = {}
    for row in reader:
        rs = (row.get(c_rs) or "").strip().strip('"').lower()
        if not RS_RE.match(rs):
            continue
        chrom = (row.get(c_ch) or "").strip().strip('"')
        pos = (row.get(c_po) or "").strip().strip('"')
        if c_gt:
            gt = (row.get(c_gt) or "").strip().strip('"')
        else:
            a1 = (row.get(c_a1) or "").strip().strip('"')
            a2 = (row.get(c_a2) or "").strip().strip('"')
            gt = a1 + a2
        out[rs] = {"chrom": chrom, "pos": pos, "genotype": gt}
    return out


# ---------- ClinVar parser ----------

# Significance strings that mean "no classification" — dropped at parse time.
CLINVAR_UNCLASSIFIED = {
    "",
    "-",
    "not provided",
    "no classification provided",
    "no classifications from unflagged records",
    "no classification for the single variant",
    "no interpretation for the single variant",
}


def _primary_sig(significance: str) -> str:
    """First semicolon-separated token of a ClinVar significance string,
    lowercased and trimmed."""
    return (significance or "").split(";")[0].split(",")[0].strip().lower()


def parse_clinvar(
    path: Path, rsid_whitelist: set[str]
) -> dict[str, list[dict[str, Any]]]:
    """Return {rsid -> [entry, ...]} for rsids in the whitelist.

    Keeps only GRCh38 rows with a classification and with complete VCF-style
    coordinates (PositionVCF, ReferenceAlleleVCF, AlternateAlleleVCF).
    """
    out: dict[str, list[dict[str, Any]]] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as f:
        header_line = f.readline()
        if not header_line.startswith("#"):
            raise RuntimeError(
                f"{path}: unexpected header (expected '#AlleleID...')"
            )
        headers = header_line.lstrip("#").rstrip("\n").split("\t")
        reader = csv.DictReader(f, fieldnames=headers, delimiter="\t")
        for row in reader:
            rs = (row.get("RS# (dbSNP)") or "").strip()
            if not rs or rs == "-1":
                continue
            rsid = "rs" + rs
            if rsid not in rsid_whitelist:
                continue
            if (row.get("Assembly") or "") != "GRCh38":
                continue

            sig_raw = (row.get("ClinicalSignificance") or "").strip()
            if _primary_sig(sig_raw) in CLINVAR_UNCLASSIFIED:
                continue

            pos_vcf = (row.get("PositionVCF") or "").strip()
            ref_vcf = (row.get("ReferenceAlleleVCF") or "").strip()
            alt_vcf = (row.get("AlternateAlleleVCF") or "").strip()
            chrom_raw = (row.get("Chromosome") or "").strip()
            if not (chrom_raw and pos_vcf and ref_vcf and alt_vcf):
                continue
            if ref_vcf in {"na", "NA", "-"} or alt_vcf in {"na", "NA", "-"}:
                continue

            out.setdefault(rsid, []).append({
                "chrom": chrom_raw,
                "pos": int(pos_vcf),
                "ref": ref_vcf,
                "alt": alt_vcf,
                "gene": (row.get("GeneSymbol") or "").strip(),
                "significance": sig_raw,
                "phenotypes": (row.get("PhenotypeList") or "").strip(),
                "review_status": (row.get("ReviewStatus") or "").strip(),
                "name": (row.get("Name") or "").strip(),
            })
    return out


# ---------- Candidate builder ----------

def variant_id_for(entry: dict[str, Any]) -> str:
    """Construct EVEE's variant_id for a ClinVar entry.

    EVEE uses 0-based coordinates; ClinVar's PositionVCF is the standard
    1-based VCF coordinate, so we subtract 1.
    """
    c = entry["chrom"]
    if not c.startswith("chr"):
        c = "chr" + c
    return f"{c}:{int(entry['pos']) - 1}:{entry['ref']}:{entry['alt']}"


def cmd_build(args: argparse.Namespace) -> None:
    print(f"Parsing {args.input} …", flush=True)
    dna = parse_dna_csv(args.input)
    print(f"  {len(dna):,} rsids in input")

    print(f"Parsing {args.clinvar} (filtered to input rsids) …", flush=True)
    cv = parse_clinvar(args.clinvar, set(dna.keys()))
    n_entries = sum(len(v) for v in cv.values())
    print(f"  {len(cv):,} rsids matched, {n_entries:,} ClinVar entries")

    candidates: list[dict[str, Any]] = []
    seen: set[str] = set()
    for rsid, entries in cv.items():
        gt = dna.get(rsid, {}).get("genotype", "")
        for e in entries:
            vid = variant_id_for(e)
            if vid in seen:
                continue
            seen.add(vid)
            candidates.append({
                "rsid": rsid,
                "variant_id": vid,
                "genotype": gt,
                "clinvar": e,
            })
    print(f"  {len(candidates):,} unique candidate variants")

    args.candidates.write_text(json.dumps(candidates, indent=2), encoding="utf-8")
    print(f"Wrote {args.candidates}")


# ---------- Zygosity + carrier detection ----------
#
# Consumer-DNA raw-data files from the major vendors (MyHeritage, 23andMe,
# AncestryDNA, FTDNA, Living DNA) all report genotypes on the forward (+)
# strand of GRCh37. ClinVar's PositionVCF / ReferenceAlleleVCF /
# AlternateAlleleVCF columns are on the forward strand of GRCh38. The
# forward strand at a given physical locus is the same strand in both
# assemblies, so we compare alleles DIRECTLY — no complement / strand-flip
# logic.
#
# Non-matching genotypes become "other:<gt>" and are treated as non-carriers,
# which is conservative against multi-allelic sites and sequencing errors.
# Palindromic SNVs (A/T, C/G) are strand-ambiguous in principle; we still
# classify them normally but flag them in the output.


def is_palindromic(ref: str, alt: str) -> bool:
    pair = {(ref or "").upper(), (alt or "").upper()}
    return pair == {"A", "T"} or pair == {"C", "G"}


def zygosity(ref: str, alt: str, gt: str) -> tuple[str, bool]:
    """Return (label, user_carries_alt).

    Label is one of: hom_ref, het, hom_alt, unknown, other:<gt>.
    """
    if not gt or len(gt) != 2:
        return "unknown", False
    gt = gt.upper()
    ref = (ref or "").upper()
    alt = (alt or "").upper()
    alleles = {gt[0], gt[1]}
    if alleles == {ref}:
        return "hom_ref", False
    if alleles == {alt}:
        return "hom_alt", True
    if alleles == {ref, alt}:
        return "het", True
    return f"other:{gt}", False


def _user_carries(rec: dict) -> bool:
    cv = rec.get("clinvar") or {}
    _, carries = zygosity(cv.get("ref", ""), cv.get("alt", ""), rec.get("genotype", ""))
    return carries


# ---------- Classification filters ----------

BENIGN = {
    "benign",
    "likely benign",
    "benign/likely benign",
}

PATHO_ONLY = {
    "pathogenic",
    "pathogenic/likely pathogenic",
    "likely pathogenic",
}


def _passes_filters(rec: dict, filter_mode: str, carriers_only: bool) -> bool:
    sig = _primary_sig((rec.get("clinvar") or {}).get("significance", ""))
    if filter_mode == "pathogenic" and sig not in PATHO_ONLY:
        return False
    if filter_mode == "non_benign" and sig in BENIGN:
        return False
    if carriers_only:
        carries = rec.get("user_carries_alt")
        if carries is None:
            carries = _user_carries(rec)
        if not carries:
            return False
    return True


# ---------- HTTP + rate limiter ----------


class RateLimiter:
    """Thread-safe per-process rate limiter."""

    def __init__(self, qps: float):
        self.min_interval = 1.0 / qps if qps > 0 else 0.0
        self.lock = threading.Lock()
        self.next_slot = 0.0

    def wait(self) -> None:
        if self.min_interval <= 0:
            return
        with self.lock:
            now = time.monotonic()
            delay = self.next_slot - now
            if delay > 0:
                time.sleep(delay)
                now = time.monotonic()
            self.next_slot = now + self.min_interval


def http_get_json(url: str, timeout: float = 30.0) -> tuple[int, Any]:
    """GET url, return (status_code, parsed_json_or_raw_string).

    Network-level failures are returned as (0, error_string).
    """
    req = urllib.request.Request(
        url,
        headers={"Accept": "application/json", "User-Agent": USER_AGENT},
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            data = resp.read()
            code = resp.getcode()
            try:
                return code, json.loads(data)
            except json.JSONDecodeError:
                return code, data.decode("utf-8", errors="replace")
    except urllib.error.HTTPError as e:
        body = e.read().decode("utf-8", errors="replace") if e.fp else ""
        return e.code, body
    except (urllib.error.URLError, TimeoutError, OSError) as e:
        return 0, str(e)


# ---------- EVEE fetcher ----------
#
# EVEE's /variants/{id}/analysis endpoint is asynchronous. A first request
# for a cold variant kicks off a Claude-Sonnet-4-6 generation and returns a
# 200/202 with {"status": "processing", "retry_after": 10}. Subsequent polls
# return the same shape until the model finishes, at which point the response
# is {"status": "complete", "result": {...}}. Results are cached persistently
# on EVEE's backend so subsequent lookups for the same variant are instant.

READY_STATES = {"ok", "complete"}
QUEUED_STATES = {
    "queued", "pending", "in_progress", "running", "processing", "started",
}

# Client-side poll schedule (seconds). We don't blindly honor the server's
# retry_after (typically 10) because most variants finish within ~3–8s and
# long polls leave workers idle. Still polite — first poll at 3s.
POLL_SCHEDULE = (3.0, 4.0, 5.0, 7.0, 10.0, 12.0, 15.0)


def _is_ready(body: Any) -> bool:
    if not isinstance(body, dict):
        return False
    if (body.get("status") or "").lower() in READY_STATES:
        return True
    return bool(body.get("result"))


def _is_queued(body: Any) -> bool:
    return (
        isinstance(body, dict)
        and (body.get("status") or "").lower() in QUEUED_STATES
    )


def fetch_evee_analysis(
    variant_id: str,
    rsid: str,
    limiter: RateLimiter,
    api_base: str,
    max_poll_s: float = DEFAULT_MAX_POLL_SECONDS,
) -> dict[str, Any]:
    """Fetch the EVEE analysis for a variant, polling as needed.

    Handles:
      * 200 with completed result → status="ok"
      * 202 / 200 with queued/processing status → poll until ready or timeout
      * 404 → fall back to search-by-rsid once (handles multi-allelic sites
        where the alt we constructed doesn't match EVEE's canonical allele)
      * 429 / 5xx / network errors → exponential backoff, up to 5 retries
    """
    current_vid = variant_id
    found_via = "direct"
    tried_search_fallback = False
    deadline = time.monotonic() + max_poll_s
    backoff_tries = 0
    poll_idx = 0

    while True:
        limiter.wait()
        url = f"{api_base}/variants/{urllib.parse.quote(current_vid, safe='')}/analysis"
        code, body = http_get_json(url)

        if code == 200 and _is_ready(body):
            return {
                "status": "ok",
                "code": 200,
                "body": body,
                "resolved_variant_id": current_vid,
                "found_via": found_via,
            }

        if code == 202 or (code == 200 and _is_queued(body)):
            wait_s = POLL_SCHEDULE[min(poll_idx, len(POLL_SCHEDULE) - 1)]
            poll_idx += 1
            if time.monotonic() + wait_s > deadline:
                return {
                    "status": "error:poll_timeout",
                    "code": code,
                    "resolved_variant_id": current_vid,
                    "found_via": found_via,
                }
            time.sleep(wait_s)
            continue

        if code == 404:
            if not tried_search_fallback:
                tried_search_fallback = True
                limiter.wait()
                s_url = f"{api_base}/variants/search?q={urllib.parse.quote(rsid)}"
                s_code, s_body = http_get_json(s_url)
                if s_code == 200 and isinstance(s_body, list) and s_body:
                    v = s_body[0].get("v")
                    if v and v != current_vid:
                        current_vid = v
                        found_via = "rsid_search"
                        continue
            return {"status": "not_in_evee", "code": 404}

        if code in (429, 500, 502, 503, 504, 0):
            if backoff_tries >= 5:
                return {
                    "status": f"error:transient_{code}",
                    "code": code,
                    "body": str(body)[:200],
                }
            time.sleep(min(30.0, 2 ** backoff_tries))
            backoff_tries += 1
            continue

        return {"status": f"error:http_{code}", "code": code, "body": str(body)[:200]}


# ---------- Cache ----------


def read_cache(path: Path) -> dict[str, dict[str, Any]]:
    if not path.exists():
        return {}
    out: dict[str, dict[str, Any]] = {}
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
                out[rec["variant_id"]] = rec
            except Exception:
                continue
    return out


def append_cache(path: Path, rec: dict[str, Any], lock: threading.Lock) -> None:
    line = json.dumps(rec, ensure_ascii=False)
    with lock:
        with path.open("a", encoding="utf-8") as f:
            f.write(line + "\n")


# ---------- Fetch command ----------


def cmd_fetch(args: argparse.Namespace) -> None:
    if not args.candidates.exists():
        print(f"Missing {args.candidates}; run `build` first", file=sys.stderr)
        sys.exit(1)
    candidates = json.loads(args.candidates.read_text(encoding="utf-8"))

    before = len(candidates)
    candidates = [
        c for c in candidates if _passes_filters(c, args.filter, args.carriers_only)
    ]
    carriers_label = "--carriers-only" if args.carriers_only else "--no-carriers-only"
    print(
        f"After filter (--filter {args.filter}, {carriers_label}): "
        f"{len(candidates):,} / {before:,}"
    )

    if args.sample:
        candidates = candidates[: args.sample]
        print(f"SAMPLE MODE: first {len(candidates)} candidates")

    cached = read_cache(args.cache)
    todo = [c for c in candidates if c["variant_id"] not in cached]
    print(
        f"{len(candidates):,} total, {len(cached):,} already cached, "
        f"{len(todo):,} to fetch"
    )
    if not todo:
        return

    limiter = RateLimiter(args.qps)
    write_lock = threading.Lock()
    t0 = time.monotonic()
    ok = miss = err = 0
    done = 0

    def worker(cand: dict[str, Any]) -> dict[str, Any]:
        r = fetch_evee_analysis(cand["variant_id"], cand["rsid"], limiter, args.api_base)
        rec = {
            "rsid": cand["rsid"],
            "variant_id": cand["variant_id"],
            "resolved_variant_id": r.get("resolved_variant_id"),
            "found_via": r.get("found_via"),
            "genotype": cand.get("genotype", ""),
            "clinvar": cand.get("clinvar", {}),
            "evee_status": r.get("status"),
        }
        if r.get("status") == "ok":
            result = r["body"].get("result") or {}
            rec["evee"] = {
                "confidence": result.get("confidence"),
                "model": result.get("model"),
                "mechanism": result.get("mechanism"),
                "summary": result.get("summary"),
                "key_evidence": result.get("key_evidence"),
                "variant_id": result.get("variant_id"),
                "generated_at": result.get("generated_at"),
            }
        elif (r.get("status") or "").startswith("error"):
            rec["error"] = r.get("body") or r.get("status")
        append_cache(args.cache, rec, write_lock)
        return rec

    with ThreadPoolExecutor(max_workers=args.concurrency) as ex:
        futures = [ex.submit(worker, c) for c in todo]
        for fut in as_completed(futures):
            rec = fut.result()
            done += 1
            status = rec.get("evee_status")
            if status == "ok":
                ok += 1
            elif status == "not_in_evee":
                miss += 1
            else:
                err += 1
            if done % 25 == 0 or done == len(todo):
                elapsed = time.monotonic() - t0
                rate = done / elapsed if elapsed else 0.0
                remaining = (len(todo) - done) / rate if rate else 0.0
                print(
                    f"  {done}/{len(todo)}  ok={ok} miss={miss} err={err}  "
                    f"{rate:.1f} req/s  eta={remaining:.0f}s",
                    flush=True,
                )

    print(f"Done in {time.monotonic() - t0:.1f}s; cache at {args.cache}")


# ---------- Emit command ----------

PATHO_RANK = {
    "pathogenic": 0,
    "pathogenic/likely pathogenic": 1,
    "likely pathogenic": 2,
    "conflicting interpretations of pathogenicity": 3,
    "conflicting classifications of pathogenicity": 3,
    "risk factor": 4,
    "drug response": 5,
    "association": 6,
    "uncertain significance": 7,
    "likely benign": 8,
    "benign/likely benign": 9,
    "benign": 10,
}


def _sort_key(r: dict) -> tuple:
    cv = r.get("clinvar") or {}
    sig = _primary_sig(cv.get("significance", ""))
    carrier = 0 if r.get("user_carries_alt") else 1
    return (carrier, PATHO_RANK.get(sig, 99), cv.get("gene", ""), r["rsid"])


def cmd_emit(args: argparse.Namespace) -> None:
    cache = read_cache(args.cache)

    for r in cache.values():
        cv = r.get("clinvar") or {}
        z, carries = zygosity(cv.get("ref", ""), cv.get("alt", ""), r.get("genotype", ""))
        r["zygosity"] = z
        r["user_carries_alt"] = carries
        r["palindromic"] = is_palindromic(cv.get("ref", ""), cv.get("alt", ""))

    ok_recs = [
        r for r in cache.values()
        if r.get("evee_status") == "ok"
        and _passes_filters(r, args.filter, args.carriers_only)
    ]
    not_in_evee = sum(1 for r in cache.values() if r.get("evee_status") == "not_in_evee")
    errored = sum(
        1 for r in cache.values() if (r.get("evee_status") or "").startswith("error")
    )

    ok_recs.sort(key=_sort_key)

    args.out_json.write_text(
        json.dumps(ok_recs, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(f"Wrote {args.out_json} ({len(ok_recs)} records)")

    sig_counts = Counter(
        _primary_sig((r.get("clinvar") or {}).get("significance", "")) for r in ok_recs
    )
    zyg_counts = Counter(r.get("zygosity", "?") for r in ok_recs)
    palindromic_count = sum(1 for r in ok_recs if r.get("palindromic"))
    carriers_label = "yes" if args.carriers_only else "no"

    md: list[str] = [
        "# EVEE personalized variant report",
        "",
        "Generated by `evee_pipeline.py`.",
        "",
        "**This is not medical advice.** ClinVar classifications change; LLM",
        "interpretations can be wrong; susceptibility / risk-factor / drug-response",
        "labels are not the same as high-penetrance disease alleles. Verify any",
        "actionable result with a qualified clinician.",
        "",
        "## Run parameters",
        "",
        f"- `--filter`: **{args.filter}**",
        f"- `--carriers-only`: **{carriers_label}**",
        "",
        "## Summary",
        "",
        f"- Records in report: **{len(ok_recs)}**",
        f"- Queried but not in EVEE: **{not_in_evee}**",
        f"- Errors: **{errored}**",
        f"- Palindromic (A/T or C/G) SNVs flagged as strand-ambiguous: **{palindromic_count}**",
        "",
        "### By ClinVar classification",
        "",
    ]
    for sig, n in sig_counts.most_common():
        md.append(f"- {sig}: **{n}**")
    md += ["", "### By zygosity", ""]
    for z, n in zyg_counts.most_common():
        md.append(f"- {z}: **{n}**")
    md += ["", "---", ""]

    for r in ok_recs:
        _append_md_record(md, r)

    args.out_md.write_text("\n".join(md), encoding="utf-8")
    print(f"Wrote {args.out_md}")


def _append_md_record(lines: list[str], r: dict) -> None:
    cv = r.get("clinvar") or {}
    ev = r.get("evee") or {}
    gene = cv.get("gene") or "?"
    lines.append(f"### {gene} · {r['rsid']} · `{r['variant_id']}`")
    lines.append("")

    tag = r.get("zygosity", "?")
    carried = " — **you carry this**" if r.get("user_carries_alt") else ""
    palin = "  *(palindromic SNV — strand-ambiguous)*" if r.get("palindromic") else ""
    lines.append(f"- **Genotype:** `{r.get('genotype','')}`  ({tag}){carried}{palin}")
    lines.append(f"- **ClinVar:** {cv.get('significance','')}  ({cv.get('review_status','')})")
    if cv.get("name"):
        lines.append(f"- **HGVS:** `{cv['name']}`")
    if cv.get("phenotypes"):
        phen = cv["phenotypes"].replace("|", "; ")
        lines.append(f"- **Phenotypes:** {phen}")
    if ev.get("confidence"):
        lines.append(
            f"- **EVEE:** {ev.get('confidence')} confidence "
            f"(model `{ev.get('model','')}`)"
        )
    lines.append("")

    if ev.get("mechanism"):
        lines.append(f"**Mechanism.** {ev['mechanism']}")
        lines.append("")
    if ev.get("summary"):
        lines.append(f"**Summary.** {ev['summary']}")
        lines.append("")
    if ev.get("key_evidence"):
        lines.append("**Key evidence:**")
        for e in ev["key_evidence"]:
            lines.append(f"- {e}")
        lines.append("")
    lines.append("---")
    lines.append("")


# ---------- all + CLI ----------


def cmd_all(args: argparse.Namespace) -> None:
    cmd_build(args)
    cmd_fetch(args)
    cmd_emit(args)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="evee_pipeline.py",
        description=(
            "Fetch EVEE variant-effect interpretations for the variants in a "
            "consumer DNA raw-data export (MyHeritage / 23andMe / AncestryDNA / "
            "FTDNA / Living DNA), filtered via ClinVar. Not medical advice."
        ),
    )
    p.add_argument(
        "mode",
        choices=["build", "fetch", "emit", "all"],
        help=(
            "build: parse CSV + ClinVar; fetch: query EVEE; "
            "emit: produce JSON/Markdown; all: run all three in order."
        ),
    )

    paths = p.add_argument_group("paths")
    paths.add_argument(
        "--input", type=Path, default=DEFAULT_INPUT,
        help=(
            f"Consumer DNA raw-data file (MyHeritage / 23andMe / AncestryDNA / "
            f"FTDNA / Living DNA). Default: {DEFAULT_INPUT}"
        ),
    )
    paths.add_argument(
        "--clinvar", type=Path, default=DEFAULT_CLINVAR,
        help=f"Gzipped ClinVar variant_summary. Default: {DEFAULT_CLINVAR}",
    )
    paths.add_argument(
        "--candidates", type=Path, default=DEFAULT_CANDIDATES,
        help=f"Candidate list written/read by build/fetch. Default: {DEFAULT_CANDIDATES}",
    )
    paths.add_argument(
        "--cache", type=Path, default=DEFAULT_CACHE,
        help=f"NDJSON resumable cache. Default: {DEFAULT_CACHE}",
    )
    paths.add_argument(
        "--out-json", type=Path, default=DEFAULT_OUT_JSON,
        help=f"Final JSON report. Default: {DEFAULT_OUT_JSON}",
    )
    paths.add_argument(
        "--out-md", type=Path, default=DEFAULT_OUT_MD,
        help=f"Final Markdown report. Default: {DEFAULT_OUT_MD}",
    )

    filters = p.add_argument_group("filters (fetch + emit)")
    filters.add_argument(
        "--filter",
        choices=["non_benign", "all", "pathogenic"],
        default="non_benign",
        help=(
            "Which ClinVar significance classes to include. "
            "non_benign (default) keeps everything except benign/LB/B+LB. "
            "all keeps every class. pathogenic keeps only P/LP."
        ),
    )
    cgroup = filters.add_mutually_exclusive_group()
    cgroup.add_argument(
        "--carriers-only", dest="carriers_only",
        action="store_true", default=True,
        help="Only include variants you carry at least one alt allele of (default: on)",
    )
    cgroup.add_argument(
        "--no-carriers-only", dest="carriers_only", action="store_false",
        help="Include non-carrier variants too",
    )
    filters.add_argument(
        "--sample", type=int, default=0,
        help="Smoke test: only process the first N candidates",
    )

    net = p.add_argument_group("network (fetch)")
    net.add_argument(
        "--api-base",
        default=os.environ.get("EVEE_API_BASE", DEFAULT_API_BASE),
        help=f"EVEE API base URL. Default: {DEFAULT_API_BASE} (env: EVEE_API_BASE)",
    )
    net.add_argument(
        "--concurrency", type=int,
        default=int(os.environ.get("EVEE_CONCURRENCY", str(DEFAULT_CONCURRENCY))),
        help=f"Max concurrent workers. Default: {DEFAULT_CONCURRENCY} (env: EVEE_CONCURRENCY)",
    )
    net.add_argument(
        "--qps", type=float,
        default=float(os.environ.get("EVEE_QPS", str(DEFAULT_QPS))),
        help=f"Soft req-per-second limit. Default: {DEFAULT_QPS} (env: EVEE_QPS)",
    )
    return p


def main() -> None:
    args = build_parser().parse_args()
    {"build": cmd_build, "fetch": cmd_fetch, "emit": cmd_emit, "all": cmd_all}[args.mode](args)


if __name__ == "__main__":
    main()
