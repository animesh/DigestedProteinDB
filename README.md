# Digested Peptide Database

Search enzymatically digested peptides by precursor mass, fragment ions, and species.
Handles SwissProt to full TrEMBL scale. No Java, no server, no RocksDB.

## How it works

Three Python scripts replace the original Java + RocksDB stack:

| File | What it does |
|---|---|
| `common.py` | Mass tables, trypsin/chymotrypsin digestion, b/y ion computation, fragment scoring |
| `build_db.py` | Streams UniProt XML.gz, digests proteins, writes sorted Parquet chunks, merges into a single DuckDB file |
| `search_mgf.py` | Parses MGF, fires a single batch range join against the DB, scores fragment ions, writes TSV |

The database is a plain DuckDB file with two tables:

- `peptides(mass_int INT, seq TEXT, accid INT, taxid INT)` sorted by `mass_int`
- `accessions(accid INT, acc TEXT, taxid INT)`

`mass_int` is `round(neutral_monoisotopic_mass * 1e4)` so a 20 ppm window at 1500 Da is about 300 integer units wide. DuckDB's row-group statistics on the sorted column make range scans fast without a separate index.

## Install

```bash
pip install duckdb pyarrow lxml
```

Python 3.10+ required. No other dependencies.

## Build the database

```bash
python build_db.py ^
    --fasta L:\promec\FastaDB\uniprot_trembl.fasta ^
    --db C:\Users\animeshs\tmp\combined.duckdb ^
    --tmp-dir C:\Users\animeshs\tmp ^
    --memory 28GB --threads 12 --bands 20 ^
    --log build.log --merge-only
# TrEMBL bacteria only (~50 GB on disk, a few hours)
python build_db.py \
    --uniprot uniprot_trembl.xml.gz \
    --db      trembl_bacteria.duckdb \
    --taxid   2 \
    --nodes   nodes.dmp \
    --mc 1 --min 7 --max 30
# TrEMBL viruses + plants
python build_db.py \
    --uniprot uniprot_trembl.xml.gz \
    --db      trembl_vp.duckdb \
    --taxid   10239 33090 \
    --nodes   nodes.dmp
```

NCBI taxonomy files (`nodes.dmp`, `names.dmp`) are in `taxdump.tar.gz` at
`https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz`.

### Build parameters

| Parameter | Default | Description |
|---|---|---|
| `--uniprot` | required | UniProt XML or XML.gz |
| `--db` | required | Output DuckDB file |
| `--enzyme` | `trypsin` | `trypsin` or `chymotrypsin` |
| `--min` | `7` | Minimum peptide length |
| `--max` | `30` | Maximum peptide length |
| `--mc` | `1` | Max missed cleavages |
| `--taxid` | all | Restrict to descendants of these NCBI taxids |
| `--nodes` | none | `nodes.dmp` (required when `--taxid` is used) |

Memory during build stays constant regardless of input size: proteins are streamed
one at a time, peptide rows are written in Parquet chunks of 5 M rows, and
DuckDB does the final merge-sort in one parallel pass.

## Search an MGF file

```bash
python search_mgf.py \
    --db      trembl_bacteria.duckdb \
    --mgf     sample.mgf \
    --ppm     20 \
    --ptms \
    --taxid   1763 1760 \
    --names   names.dmp \
    --min-ions 3 \
    --out     results.tsv
```

### Search parameters

| Parameter | Default | Description |
|---|---|---|
| `--db` | required | DuckDB file from `build_db.py` |
| `--mgf` | required | Input MGF file |
| `--ppm` | `20` | Precursor mass tolerance in ppm |
| `--ptms` | off | Also search with common PTM mass offsets (see below) |
| `--taxid` | all | Keep hits from these exact NCBI taxids only |
| `--min-ions` | `3` | Min matched fragment ions to report a hit (only applied when MS2 is present) |
| `--names` | none | `names.dmp` for a species name column in the output |
| `--out` | stdout | Output TSV path |

### PTM offsets

When `--ptms` is set, each spectrum is searched against additional mass
windows shifted by the following deltas. A hit from a PTM window is labelled
accordingly in the output, and fragment scoring tries placing the modification
on every applicable residue to localise the best position.

| Label | Delta (Da) | Residues |
|---|---|---|
| `ox_M` | +15.994915 | M |
| `phospho` | +79.966331 | S, T, Y |
| `deamid` | +0.984016 | N, Q |
| `acetyl_K` | +42.010565 | K |
| `carbam_C` | +57.021464 | C |

### Output columns

```
scan_id  title  precursor_mass_da  charge  peptide_seq  accession  taxid  species
ppm_err  ptm  ptm_pos  n_matched_ions  frag_score
```

- `ptm` is `none` for unmodified hits
- `ptm_pos` is the 0-indexed position of the best-localised modification (`-1` if none)
- `frag_score` is the fraction of total observed fragment intensity matched
- `n_matched_ions` counts matched singly-charged b/y ions within the ppm tolerance

## Species filtering

`--taxid` in `search_mgf.py` takes exact NCBI taxids. For a whole clade
(e.g. all Mycobacterium, taxid 1763), the most efficient approach is to filter
at build time:

```bash
python build_db.py --uniprot trembl.xml.gz --db myco.duckdb \
    --taxid 1763 --nodes nodes.dmp
```

If you already have a full-TrEMBL database and want to filter at query time,
pass all descendant taxids to `--taxid`. You can generate them with:

```python
from build_db import load_ncbi_descendants
taxids = load_ncbi_descendants('nodes.dmp', [1763])
# pass as: --taxid 1763 <all descendants>
```

## Fragment ion scoring

Scoring is straightforward: for each candidate peptide compute singly-charged
b and y ions, count how many land within `--ppm` of an observed peak, and sum
the matched intensities. `frag_score` is that sum divided by total spectrum
intensity.

For PTM candidates the scorer tries each possible modification site and reports
the position with the most matched ions. Double modifications and multiply
charged fragment ions are not handled; that level of scoring belongs in a
dedicated engine like MSFragger or Mascot.

## Known limitations

- Fragment scoring covers z=1 b/y ions only.
- Double modifications on the same peptide are not modelled.
- Clade expansion at search time requires pre-computing descendants externally.
- The DuckDB range join for very large MGF files (> 50k spectra) has not been
  benchmarked at full TrEMBL scale yet. If it becomes a bottleneck, splitting
  the MGF into chunks and running in parallel is the straightforward fix.

## File sizes and timing (approximate)

| Dataset | DB size | Build time | Search 1000 spectra |
|---|---|---|---|
| SwissProt (~570k proteins) | ~2 GB | ~10 min | < 5 s |
| TrEMBL bacteria (~50M proteins) | ~50 GB | ~3 h | < 60 s |
| Full TrEMBL | ~500 GB | ~24 h | < 10 min |

Timings on a 16-core workstation with NVMe storage. Build is I/O and
single-thread XML-parse bound. Search is DuckDB multicore.

## TODO

- Benchmark and tune DuckDB range join at full TrEMBL scale
- Add z=2 fragment ions
- Add decoy generation for FDR estimation
- Expose `load_ncbi_descendants` as a CLI helper so clade filtering at search
  time is one extra flag rather than a Python call
- Downloadable pre-built databases for common taxonomy subsets