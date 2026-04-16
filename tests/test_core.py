#!/usr/bin/env python3
"""Regression tests for the pure-function logic in evee_pipeline.

These cover the decision logic the pipeline actually relies on — coordinate
offset, zygosity, palindromic detection, filter predicates, ClinVar
significance parsing, EVEE status classification, and the sort ordering
used for the final report.

No network, no files, no randomness — every case is an explicit fixture.
Run with `python3 tests/test_core.py`; exit code is non-zero if anything
fails.
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from evee_pipeline import (  # noqa: E402
    PATHO_RANK,
    _is_queued,
    _is_ready,
    _passes_filters,
    _primary_sig,
    _sort_key,
    is_palindromic,
    variant_id_for,
    zygosity,
)


def _check(name: str, got: object, expected: object) -> bool:
    if got == expected:
        print(f"  PASS  {name}")
        return True
    print(f"  FAIL  {name}: expected {expected!r}, got {got!r}")
    return False


# ---------- variant_id_for ----------
#
# ClinVar PositionVCF is 1-based; EVEE's variant_id is 0-based. The coord
# conversion is a single subtraction but it's the most consequential line of
# arithmetic in the whole pipeline, so it gets explicit coverage.

def test_variant_id_for() -> int:
    cases = [
        ({"chrom": "1", "pos": 100, "ref": "A", "alt": "G"}, "chr1:99:A:G"),
        ({"chrom": "chr1", "pos": 100, "ref": "A", "alt": "G"}, "chr1:99:A:G"),
        ({"chrom": "X", "pos": 1, "ref": "C", "alt": "T"}, "chrX:0:C:T"),
        ({"chrom": "Y", "pos": 10000, "ref": "G", "alt": "A"}, "chrY:9999:G:A"),
        ({"chrom": "MT", "pos": 200, "ref": "A", "alt": "G"}, "chrMT:199:A:G"),
        # Matches the README's worked example: rs1800896 IL10 promoter.
        (
            {"chrom": "1", "pos": 206773552, "ref": "T", "alt": "C"},
            "chr1:206773551:T:C",
        ),
        # String pos should also work (int() coerces it).
        ({"chrom": "1", "pos": "50", "ref": "A", "alt": "T"}, "chr1:49:A:T"),
    ]
    failed = 0
    for entry, expected in cases:
        if not _check(f"variant_id_for({entry})", variant_id_for(entry), expected):
            failed += 1
    return failed


# ---------- is_palindromic ----------

def test_is_palindromic() -> int:
    cases = [
        (("A", "T"), True),
        (("T", "A"), True),
        (("C", "G"), True),
        (("G", "C"), True),
        (("a", "t"), True),        # case-insensitive
        (("A", "C"), False),
        (("A", "G"), False),
        (("C", "T"), False),
        (("G", "T"), False),
        (("", "A"), False),        # malformed allele — not palindromic
        (("A", ""), False),
    ]
    failed = 0
    for (ref, alt), expected in cases:
        if not _check(f"is_palindromic({ref!r},{alt!r})", is_palindromic(ref, alt), expected):
            failed += 1
    return failed


# ---------- zygosity ----------
#
# Documents current behavior, including that single-character genotypes
# (e.g. hemizygous calls on male X/Y) fall through to "unknown" rather than
# being classified as hemi_ref / hemi_alt. If that's ever changed,
# regenerate these expectations deliberately.

def test_zygosity() -> int:
    cases = [
        # (ref, alt, gt, expected)
        ("A", "G", "AA", ("hom_ref", False)),
        ("A", "G", "GG", ("hom_alt", True)),
        ("A", "G", "AG", ("het", True)),
        ("A", "G", "GA", ("het", True)),
        ("A", "G", "ag", ("het", True)),     # lowercase genotype
        ("a", "g", "AG", ("het", True)),     # lowercase ref/alt
        ("A", "G", "", ("unknown", False)),
        ("A", "G", "A", ("unknown", False)),    # hemizygous single-char
        ("A", "G", "AAA", ("unknown", False)),  # unexpectedly long
        ("A", "G", "AT", ("other:AT", False)),  # non-matching allele
        ("A", "G", "--", ("other:--", False)),  # no-call from MyHeritage/23andMe
        ("A", "G", "00", ("other:00", False)),  # no-call from AncestryDNA
    ]
    failed = 0
    for ref, alt, gt, expected in cases:
        if not _check(f"zygosity({ref!r},{alt!r},{gt!r})", zygosity(ref, alt, gt), expected):
            failed += 1
    return failed


# ---------- _primary_sig ----------

def test_primary_sig() -> int:
    cases = [
        ("Pathogenic", "pathogenic"),
        ("Pathogenic;drug response", "pathogenic"),          # semicolon
        ("Pathogenic, low penetrance", "pathogenic"),        # comma
        ("Likely pathogenic; association", "likely pathogenic"),
        ("Conflicting interpretations of pathogenicity", "conflicting interpretations of pathogenicity"),
        ("  Benign  ", "benign"),                            # surrounding whitespace
        ("", ""),
        (None, ""),                                          # None input
    ]
    failed = 0
    for inp, expected in cases:
        if not _check(f"_primary_sig({inp!r})", _primary_sig(inp), expected):
            failed += 1
    return failed


# ---------- _is_ready ----------

def test_is_ready() -> int:
    cases = [
        ({"status": "ok"}, True),
        ({"status": "complete"}, True),
        ({"status": "OK"}, True),                # case-insensitive
        ({"result": {"mechanism": "x"}}, True),  # presence of result also counts
        ({"status": "processing"}, False),
        ({"status": "queued"}, False),
        ({}, False),
        ("not a dict", False),
        (None, False),
        ([], False),
    ]
    failed = 0
    for inp, expected in cases:
        if not _check(f"_is_ready({inp!r})", _is_ready(inp), expected):
            failed += 1
    return failed


# ---------- _is_queued ----------

def test_is_queued() -> int:
    cases = [
        ({"status": "queued"}, True),
        ({"status": "pending"}, True),
        ({"status": "in_progress"}, True),
        ({"status": "running"}, True),
        ({"status": "processing"}, True),
        ({"status": "started"}, True),
        ({"status": "PROCESSING"}, True),        # case-insensitive
        ({"status": "ok"}, False),
        ({"status": "complete"}, False),
        ({"status": ""}, False),
        ({}, False),
        ("not a dict", False),
        (None, False),
    ]
    failed = 0
    for inp, expected in cases:
        if not _check(f"_is_queued({inp!r})", _is_queued(inp), expected):
            failed += 1
    return failed


# ---------- _passes_filters ----------
#
# Three filter modes × carriers-only on/off × various ClinVar classes. The
# combinations matter because they drive what ends up in the user's report.

def _rec(sig: str, gt: str, ref: str = "A", alt: str = "G") -> dict:
    return {
        "rsid": "rs1",
        "variant_id": "chr1:99:A:G",
        "genotype": gt,
        "clinvar": {"significance": sig, "ref": ref, "alt": alt},
    }


def test_passes_filters() -> int:
    patho_carrier = _rec("Pathogenic", "AG")
    benign_carrier = _rec("Benign", "AG")
    patho_noncarrier = _rec("Pathogenic", "AA")
    vus_carrier = _rec("Uncertain significance", "AG")
    likely_benign_carrier = _rec("Likely benign", "AG")
    likely_patho_carrier = _rec("Likely pathogenic", "AG")

    # Each row: (label, record, filter_mode, carriers_only, expected_passes)
    cases = [
        # pathogenic + carrier — passes every filter
        ("patho_carrier/all/noco",         patho_carrier,         "all",         False, True),
        ("patho_carrier/all/co",           patho_carrier,         "all",         True,  True),
        ("patho_carrier/non_benign/co",    patho_carrier,         "non_benign",  True,  True),
        ("patho_carrier/pathogenic/co",    patho_carrier,         "pathogenic",  True,  True),

        # benign + carrier — drops under non_benign and pathogenic
        ("benign_carrier/all/co",          benign_carrier,        "all",         True,  True),
        ("benign_carrier/non_benign/co",   benign_carrier,        "non_benign",  True,  False),
        ("benign_carrier/pathogenic/co",   benign_carrier,        "pathogenic",  True,  False),

        # pathogenic + non-carrier — drops under carriers_only
        ("patho_noncarrier/all/co",        patho_noncarrier,      "all",         True,  False),
        ("patho_noncarrier/all/noco",      patho_noncarrier,      "all",         False, True),
        ("patho_noncarrier/pathogenic/co", patho_noncarrier,      "pathogenic",  True,  False),

        # VUS passes non_benign but not pathogenic-only
        ("vus_carrier/non_benign/co",      vus_carrier,           "non_benign",  True,  True),
        ("vus_carrier/pathogenic/co",      vus_carrier,           "pathogenic",  True,  False),

        # likely benign drops under non_benign
        ("likely_benign/non_benign/co",    likely_benign_carrier, "non_benign",  True,  False),
        ("likely_benign/all/co",           likely_benign_carrier, "all",         True,  True),

        # likely pathogenic passes both pathogenic and non_benign
        ("likely_patho/pathogenic/co",     likely_patho_carrier,  "pathogenic",  True,  True),
        ("likely_patho/non_benign/co",     likely_patho_carrier,  "non_benign",  True,  True),
    ]
    failed = 0
    for label, rec, mode, co, expected in cases:
        if not _check(
            f"_passes_filters({label})", _passes_filters(rec, mode, co), expected
        ):
            failed += 1
    return failed


# ---------- _sort_key ----------
#
# Report ordering is what the user actually sees. Carriers first, then by
# pathogenicity rank, then by gene, then by rsid.

def test_sort_key() -> int:
    def rec(rsid: str, sig: str, gene: str, carries: bool) -> dict:
        return {
            "rsid": rsid,
            "user_carries_alt": carries,
            "clinvar": {"significance": sig, "gene": gene},
        }

    records = [
        rec("rs4", "Benign", "GENE_B", True),
        rec("rs1", "Pathogenic", "GENE_A", True),
        rec("rs5", "Pathogenic", "GENE_Z", False),   # non-carrier — should go last
        rec("rs2", "Likely pathogenic", "GENE_A", True),
        rec("rs3", "Uncertain significance", "GENE_C", True),
    ]
    records.sort(key=_sort_key)
    got_order = [r["rsid"] for r in records]
    expected = ["rs1", "rs2", "rs3", "rs4", "rs5"]
    return 0 if _check("_sort_key ordering", got_order, expected) else 1


# ---------- PATHO_RANK ----------

def test_patho_rank() -> int:
    failed = 0
    # Both ClinVar spellings of "conflicting" are covered.
    for k in (
        "conflicting interpretations of pathogenicity",
        "conflicting classifications of pathogenicity",
    ):
        if not _check(f"PATHO_RANK has {k!r}", k in PATHO_RANK, True):
            failed += 1
    # Pathogenic ranks higher (smaller number) than benign.
    if not _check(
        "PATHO_RANK[pathogenic] < PATHO_RANK[benign]",
        PATHO_RANK["pathogenic"] < PATHO_RANK["benign"],
        True,
    ):
        failed += 1
    # Likely pathogenic ranks between P and VUS.
    if not _check(
        "pathogenic < likely_pathogenic < uncertain_significance",
        PATHO_RANK["pathogenic"]
        < PATHO_RANK["likely pathogenic"]
        < PATHO_RANK["uncertain significance"],
        True,
    ):
        failed += 1
    return failed


ALL_TESTS = [
    ("variant_id_for", test_variant_id_for),
    ("is_palindromic", test_is_palindromic),
    ("zygosity", test_zygosity),
    ("_primary_sig", test_primary_sig),
    ("_is_ready", test_is_ready),
    ("_is_queued", test_is_queued),
    ("_passes_filters", test_passes_filters),
    ("_sort_key", test_sort_key),
    ("PATHO_RANK", test_patho_rank),
]


def main() -> int:
    print(f"Running {len(ALL_TESTS)} core-logic test group(s)\n")
    failed_total = 0
    for name, fn in ALL_TESTS:
        print(f"[{name}]")
        failed_total += fn()
    print()
    if failed_total == 0:
        print(f"All {len(ALL_TESTS)} core test groups passed.")
        return 0
    print(f"{failed_total} core test(s) failed.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
