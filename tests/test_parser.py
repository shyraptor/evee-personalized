#!/usr/bin/env python3
"""Regression tests for evee_pipeline.parse_dna_csv().

Runs parse_dna_csv() against each fixture in tests/fixtures/ and verifies the
parsed {rsid: {chrom, pos, genotype}} dict matches a known-good expectation.
Prints PASS/FAIL per fixture and exits non-zero if anything fails.
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from evee_pipeline import parse_dna_csv  # noqa: E402

FIXTURES = HERE / "fixtures"

# Eight synthetic SNVs: the baseline every standard fixture encodes.
STANDARD = {
    "rs1": {"chrom": "1", "pos": "101", "genotype": "AA"},
    "rs2": {"chrom": "1", "pos": "102", "genotype": "CC"},
    "rs3": {"chrom": "1", "pos": "103", "genotype": "GG"},
    "rs4": {"chrom": "1", "pos": "104", "genotype": "TT"},
    "rs5": {"chrom": "1", "pos": "105", "genotype": "--"},
    "rs6": {"chrom": "1", "pos": "106", "genotype": "GC"},
    "rs7": {"chrom": "1", "pos": "107", "genotype": "TC"},
    "rs8": {"chrom": "1", "pos": "108", "genotype": "AT"},
}

# AncestryDNA uses "0" as its no-call character, so rs5 genotype concatenates
# to "00" instead of "--".
ANCESTRY = {
    **STANDARD,
    "rs5": {"chrom": "1", "pos": "105", "genotype": "00"},
}

# Mitochondrial variant of the Ancestry fixture: only three rows, on
# chromosome 26.
ANCESTRY_MT = {
    "rs1": {"chrom": "26", "pos": "101", "genotype": "AA"},
    "rs2": {"chrom": "26", "pos": "102", "genotype": "00"},
    "rs3": {"chrom": "26", "pos": "103", "genotype": "GC"},
}

TESTS = [
    ("myheritage.csv", STANDARD),
    ("ftdna.csv", STANDARD),
    ("ftdna_famfinder.csv", STANDARD),   # '-' + '-' -> '--'
    ("23andme.txt", STANDARD),
    ("23andme_allele.txt", STANDARD),    # '-' + '-' -> '--'
    ("23andme_win.txt", STANDARD),
    ("ancestry.txt", ANCESTRY),
    ("ancestry_mt.txt", ANCESTRY_MT),
    ("livingdna.csv", STANDARD),
]


def run_one(name: str, expected: dict) -> bool:
    path = FIXTURES / name
    try:
        got = parse_dna_csv(path)
    except Exception as e:
        print(f"  FAIL  {name}  {type(e).__name__}: {e}")
        return False
    ok = True
    for rsid, exp in expected.items():
        if rsid not in got:
            print(f"  FAIL  {name}  missing rsid {rsid}")
            ok = False
            continue
        if got[rsid] != exp:
            print(f"  FAIL  {name}  {rsid}: expected {exp}, got {got[rsid]}")
            ok = False
    if ok:
        extra = len(got) - len(expected)
        extra_note = f" (+{extra} extra)" if extra else ""
        print(f"  PASS  {name}  ({len(expected)} rsids parsed{extra_note})")
    return ok


def main() -> int:
    print(f"Running {len(TESTS)} parser tests from {FIXTURES}\n")
    passed = 0
    failed = 0
    for name, expected in TESTS:
        if run_one(name, expected):
            passed += 1
        else:
            failed += 1
    print()
    print(f"{passed} passed, {failed} failed, {len(TESTS)} total")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
