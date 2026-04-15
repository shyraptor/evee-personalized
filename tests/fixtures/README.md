# Test fixtures

Synthetic consumer-DNA raw-data files for regression testing
`parse_dna_csv()`. Each file encodes the same eight synthetic SNVs
(rs1–rs8, chromosome 1, positions 101–108) with a known allele assignment
(or three SNVs on chromosome 26 for `ancestry_mt.txt`). They exist to
exercise format quirks — delimiter choice, commented vs uncommented header,
split vs combined allele columns, no-call conventions — not to represent
real genetic data.

**No personal genetic data is in these files.** The format specifications
themselves (field names, delimiters, header conventions) are public and
come from each vendor's raw-data export documentation; the test content
(`rs1..rs8`, fake positions) is fabricated for this repo.

## Formats covered

| File | Vendor / variant | What it exercises |
|---|---|---|
| `myheritage.csv` | MyHeritage standard | comma-separated, quoted values, `--` no-call |
| `ftdna.csv` | FamilyTreeDNA standard | like MyHeritage, no `#` preamble |
| `ftdna_famfinder.csv` | FTDNA Family Finder | commented header, `name` alias, split alleles |
| `23andme.txt` | 23andMe standard | tab-separated, commented header, `--` no-call |
| `23andme_allele.txt` | 23andMe allele-split variant | tab, split `allele1`/`allele2` |
| `23andme_win.txt` | 23andMe with CRLF line endings | Windows line endings |
| `ancestry.txt` | AncestryDNA standard | tab, uncommented header, split alleles, `0` no-call |
| `ancestry_mt.txt` | AncestryDNA mitochondrial | chromosome 26, 3-row file |
| `livingdna.csv` | Living DNA | **tab-separated inside a `.csv` extension** |

## Running the tests

```bash
python3 tests/test_parser.py
```

Prints `PASS` / `FAIL` per fixture and exits non-zero if anything fails.
