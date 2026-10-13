"""
build_db.py - Build a DuckDB peptide database from UniProt FASTA.

For uncompressed FASTA: splits file into --workers byte ranges, each
worker reads its own range independently and writes Parquet directly.
No large data passes through IPC; only tiny config args go in,
only worker stats come back.

For gzip FASTA: falls back to single-threaded streaming (byte-range
seeking does not work on gzip streams).

Downloads:
    https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
    https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz

Usage:
    python build_db.py \
        --fasta   uniprot_sprot.fasta uniprot_trembl.fasta \
        --db      combined.duckdb \
        [--enzyme trypsin|chymotrypsin] \
        [--min 7]  [--max 30]  [--mc 1] \
        [--workers 12] \
        [--memory 28GB] \
        [--threads 12] \
        [--log    build.log] \
        [--taxid 2 10239] \
        [--nodes  nodes.dmp]

Dependencies:
    pip install duckdb pyarrow
"""

import argparse
import gzip
import os
import re
import shutil
import sys
import traceback
from collections import defaultdict, deque
from datetime import datetime
from multiprocessing import Pool, cpu_count
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq

from common import digest, neutral_mass

# ── Schemas ───────────────────────────────────────────────────────────────────
# Peptide Parquet files store acc as a string; DuckDB assigns integer accids
# during the merge step so workers need no coordination.

_PEP_SCHEMA = pa.schema([
    ('mass_int', pa.int32()),
    ('seq',      pa.string()),
    ('acc',      pa.string()),
    ('taxid',    pa.int32()),
])
_ACC_SCHEMA = pa.schema([
    ('acc',   pa.string()),
    ('taxid', pa.int32()),
])

WORKER_CHUNK = 2_000_000    # peptide rows per Parquet file per worker

# ── Logging ───────────────────────────────────────────────────────────────────

_log_fh = None

def log(msg: str):
    ts   = datetime.now().strftime('%H:%M:%S')
    line = f"[{ts}] {msg}"
    print(line, file=sys.stderr, flush=True)
    if _log_fh:
        print(line, file=_log_fh, flush=True)

# ── FASTA parsing ─────────────────────────────────────────────────────────────

_OX_RE = re.compile(r'OX=(\d+)')


def _parse_header(line: str):
    """Return (accession, taxid) from a UniProt FASTA header line (without '>')."""
    fields = line.split('|')
    acc    = fields[1] if len(fields) >= 2 else line.split()[0]
    m      = _OX_RE.search(line)
    taxid  = int(m.group(1)) if m else 0
    return acc, taxid


def _iter_fasta_range(path: str, start: int, end: int):
    """
    Yield (acc, taxid, seq) for proteins whose header byte-offset is in
    [start, end).  Uses readline() so line_pos is exact before each read.
    Works only on uncompressed files where seek() is meaningful.
    """
    acc, taxid, parts = None, 0, []
    with open(path, 'r', encoding='ascii', errors='ignore') as fh:
        fh.seek(start)
        while True:
            line_pos = fh.tell()
            line     = fh.readline()
            if not line:                          # EOF
                if acc and parts:
                    yield acc, taxid, ''.join(parts)
                break
            line = line.rstrip()
            if not line:
                continue
            if line.startswith('>'):
                if acc and parts:
                    yield acc, taxid, ''.join(parts)
                if line_pos >= end:               # header belongs to next worker
                    break
                acc, taxid = _parse_header(line[1:])
                parts = []
            elif acc is not None:
                parts.append(line.upper())


def _iter_fasta_gz(path: str):
    """Single-threaded streaming for gzip FASTA (seek not supported)."""
    abs_path = str(Path(path).resolve())
    open_fn  = gzip.open if abs_path.endswith('.gz') else open
    acc, taxid, parts = None, 0, []
    with open_fn(abs_path, 'rt', encoding='ascii', errors='ignore') as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith('>'):
                if acc and parts:
                    yield acc, taxid, ''.join(parts)
                acc, taxid = _parse_header(line[1:])
                parts = []
            elif acc is not None:
                parts.append(line.upper())
    if acc and parts:
        yield acc, taxid, ''.join(parts)

# ── Chunk boundary finder ─────────────────────────────────────────────────────

def _find_chunk_boundaries(path: str, n: int) -> list:
    """
    Return n+1 byte offsets that split path into n chunks, each starting
    at a '>' line.  Reads only a 1 MB window near each split point.
    """
    size = os.path.getsize(path)
    boundaries = [0]
    with open(path, 'rb') as fh:
        for i in range(1, n):
            target = (size * i) // n
            fh.seek(max(0, target - 1))
            buf = fh.read(1 << 20)        # 1 MB
            idx = buf.find(b'\n>')
            boundaries.append(
                (max(0, target - 1) + idx + 1) if idx >= 0 else size)
    boundaries.append(size)
    return boundaries

# ── Worker (module-level for Windows multiprocessing) ─────────────────────────

def _worker_fn(args):
    """
    Reads a byte range from a plain FASTA, digests proteins, writes Parquet.
    Only a small config tuple enters via IPC; only stats exit.
    """
    (path, start, end, allowed_taxids,
     min_len, max_len, mc, enzyme, tmp_dir, worker_id) = args

    ms, seqs, accs, tids = [], [], [], []
    aacc, atids_acc = [], []
    chunk_id  = 0
    n_prot    = 0
    n_pep     = 0

    def _flush():
        nonlocal chunk_id
        # Sort by mass_int within the chunk so DuckDB merge-sort is faster
        rows = sorted(zip(ms, seqs, accs, tids), key=lambda r: r[0])
        if rows:
            ms_s, seqs_s, accs_s, tids_s = zip(*rows)
        else:
            ms_s = seqs_s = accs_s = tids_s = []
        pq.write_table(
            pa.table({
                'mass_int': pa.array(ms_s,   pa.int32()),
                'seq':      pa.array(seqs_s, pa.string()),
                'acc':      pa.array(accs_s, pa.string()),
                'taxid':    pa.array(tids_s, pa.int32()),
            }, schema=_PEP_SCHEMA),
            str(Path(tmp_dir) / f'pep_{worker_id:04d}_{chunk_id:05d}.parquet'),
            compression='zstd', write_statistics=True)
        ms.clear(); seqs.clear(); accs.clear(); tids.clear()
        chunk_id += 1

    for acc, taxid, seq in _iter_fasta_range(path, start, end):
        if allowed_taxids is not None and taxid not in allowed_taxids:
            continue

        aacc.append(acc); atids_acc.append(taxid)
        peps = digest(seq, min_len, max_len, mc, enzyme)

        for pep in peps:
            ms.append(round(neutral_mass(pep) * 1e4))
            seqs.append(pep); accs.append(acc); tids.append(taxid)

        n_prot += 1
        n_pep  += len(peps)

        if len(ms) >= WORKER_CHUNK:
            _flush()

    if ms:
        _flush()

    pq.write_table(
        pa.table({
            'acc':   pa.array(aacc,      pa.string()),
            'taxid': pa.array(atids_acc, pa.int32()),
        }, schema=_ACC_SCHEMA),
        str(Path(tmp_dir) / f'acc_{worker_id:04d}.parquet'),
        compression='zstd', write_statistics=True)

    return worker_id, n_prot, n_pep

# ── NCBI taxonomy ─────────────────────────────────────────────────────────────

def load_ncbi_descendants(nodes_dmp: str, parent_taxids: list) -> set:
    log(f"Loading NCBI taxonomy from {nodes_dmp} ...")
    children = defaultdict(list)
    with open(nodes_dmp) as fh:
        for line in fh:
            parts = line.split('|')
            if len(parts) < 2:
                continue
            try:
                children[int(parts[1].strip())].append(int(parts[0].strip()))
            except ValueError:
                pass
    allowed, queue = set(), deque(parent_taxids)
    while queue:
        tid = queue.popleft()
        if tid in allowed:
            continue
        allowed.add(tid)
        queue.extend(children.get(tid, []))
    log(f"  {len(allowed):,} allowed taxids.")
    return allowed

# ── Main build ────────────────────────────────────────────────────────────────

def _build_parallel(args, abs_paths, allowed_taxids, tmp_dir):
    """Byte-split + multiprocessing path for uncompressed FASTA."""
    n_workers = args.workers or max(1, cpu_count())
    log(f"Workers   : {n_workers}  (parallel byte-range)")

    jobs = []
    for path in abs_paths:
        log(f"  Splitting {path} ...")
        bounds = _find_chunk_boundaries(path, n_workers)
        sizes  = [f"{(e-s)//1024//1024}MB" for s, e in zip(bounds, bounds[1:])]
        log(f"    chunks: {sizes}")
        for s, e in zip(bounds, bounds[1:]):
            wid = len(jobs)
            jobs.append((path, s, e, allowed_taxids,
                         args.min, args.max, args.mc, args.enzyme,
                         str(tmp_dir), wid))

    total_prot, total_pep = 0, 0
    try:
        with Pool(n_workers) as pool:
            for wid, n_prot, n_pep in pool.imap_unordered(
                    _worker_fn, jobs, chunksize=1):
                total_prot += n_prot
                total_pep  += n_pep
                log(f"  Worker {wid:3d} done: {n_prot:,} proteins | {n_pep:,} peptides")
    except KeyboardInterrupt:
        log("Interrupted.")
        sys.exit(1)
    except Exception:
        log("ERROR during digestion:")
        log(traceback.format_exc())
        sys.exit(1)

    return total_prot, total_pep


def _build_serial(args, abs_paths, allowed_taxids, tmp_dir):
    """Single-threaded path for gzip FASTA (or --workers 1)."""
    log("Workers   : 1  (serial streaming)")

    ms, seqs, accs, tids = [], [], [], []
    aacc, atids_acc       = [], []
    chunk_id = 0
    n_prot   = 0
    n_pep    = 0

    def _flush():
        nonlocal chunk_id
        rows = sorted(zip(ms, seqs, accs, tids), key=lambda r: r[0])
        if rows:
            ms_s, seqs_s, accs_s, tids_s = zip(*rows)
        else:
            ms_s = seqs_s = accs_s = tids_s = []
        pq.write_table(
            pa.table({
                'mass_int': pa.array(ms_s,   pa.int32()),
                'seq':      pa.array(seqs_s, pa.string()),
                'acc':      pa.array(accs_s, pa.string()),
                'taxid':    pa.array(tids_s, pa.int32()),
            }, schema=_PEP_SCHEMA),
            str(tmp_dir / f'pep_0000_{chunk_id:05d}.parquet'),
            compression='zstd', write_statistics=True)
        ms.clear(); seqs.clear(); accs.clear(); tids.clear()
        chunk_id += 1

    CHUNK = 5_000_000

    for path in abs_paths:
        log(f"Parsing {path} ...")
        for acc, taxid, seq in _iter_fasta_gz(path):
            if allowed_taxids is not None and taxid not in allowed_taxids:
                continue
            aacc.append(acc); atids_acc.append(taxid)
            peps = digest(seq, args.min, args.max, args.mc, args.enzyme)
            for pep in peps:
                ms.append(round(neutral_mass(pep) * 1e4))
                seqs.append(pep); accs.append(acc); tids.append(taxid)
            n_prot += 1
            n_pep  += len(peps)
            if n_prot % 100_000 == 0:
                log(f"  {n_prot:,} proteins | {n_pep:,} peptides | buffer {len(ms):,}")
            if len(ms) >= CHUNK:
                log(f"  Flushing chunk {chunk_id} ...")
                _flush()

    if ms:
        _flush()

    pq.write_table(
        pa.table({'acc': pa.array(aacc, pa.string()),
                  'taxid': pa.array(atids_acc, pa.int32())}, schema=_ACC_SCHEMA),
        str(tmp_dir / 'acc_0000.parquet'),
        compression='zstd', write_statistics=True)

    return n_prot, n_pep


def _sort_bands(con, pep_glob: str, data_dir: Path, n_bands: int) -> str:
    """
    Partition peptide rows into n_bands non-overlapping mass_int ranges,
    sort each band independently, and write to Parquet in data_dir.

    Each band is ~total/n_bands rows.  With n_bands=20 and 6.65B rows,
    each band is ~332M rows * 33 bytes = ~11 GB: fits in 28 GB RAM,
    so DuckDB never spills to disk.

    Peak disk during this step:
        worker Parquets (read-only) + one band output (~total/n_bands compressed)
    Typically ~92 GB for TrEMBL, well within any local SSD.

    Returns the glob pattern for the written band files.
    """
    log("  Reading mass range from Parquet statistics ...")
    lo, hi = con.execute(
        f"SELECT MIN(mass_int), MAX(mass_int) FROM read_parquet('{pep_glob}')"
    ).fetchone()
    band_width = max(1, (hi - lo + n_bands - 1) // n_bands)
    log(f"  mass_int: {lo} to {hi}, {n_bands} bands of width {band_width}")

    data_dir.mkdir(exist_ok=True)

    for b in range(n_bands):
        b_lo = lo + b * band_width
        b_hi = (lo + (b + 1) * band_width) if b < n_bands - 1 else (hi + 1)
        out  = str(data_dir / f'peptides_{b:05d}.parquet').replace('\\', '/')
        log(f"  Band {b + 1}/{n_bands}: mass_int [{b_lo}, {b_hi}) ...")
        con.execute(f"""
            COPY (
                SELECT * FROM read_parquet('{pep_glob}')
                WHERE mass_int >= {b_lo} AND mass_int < {b_hi}
                ORDER BY mass_int
            ) TO '{out}' (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 500000)
        """)

    return str(data_dir / 'peptides_*.parquet').replace('\\', '/')


def _clean_stale_parquet(tmp_dir: Path):
    """
    Remove any pep_*.parquet files whose schema does not contain an 'acc'
    string column.  These are leftovers from an older build that stored an
    integer 'accid' instead, and would break the DuckDB join.
    """
    removed = 0
    for f in sorted(tmp_dir.glob('pep_*.parquet')):
        try:
            names = pq.read_schema(str(f)).names
            if 'acc' not in names:
                f.unlink()
                removed += 1
        except Exception:
            f.unlink()       # corrupt file - also remove
            removed += 1
    if removed:
        log(f"  Removed {removed} stale pep Parquet files with wrong schema.")


def build(args):
    db_path = Path(args.db)
    tmp_root = Path(args.tmp_dir) if args.tmp_dir else db_path.parent
    tmp_dir  = tmp_root / (db_path.stem + '_build_tmp')
    log(f"Temp dir  : {tmp_dir}")
    log(f"Output DB : {db_path}")

    if args.merge_only:
        # Skip digestion; just clean bad-schema leftovers and go to DuckDB.
        if not tmp_dir.exists():
            log("ERROR: --merge-only requires an existing temp dir with Parquet files.")
            sys.exit(1)
        log("--merge-only: skipping digestion, cleaning stale files ...")
        _clean_stale_parquet(tmp_dir)
        total_prot = total_pep = 0   # not re-counted; DuckDB will report final numbers
    else:
        # Fresh start: wipe any stale temp dir from a previous run.
        if tmp_dir.exists():
            log(f"Removing stale temp dir from previous run ...")
            shutil.rmtree(tmp_dir)
        tmp_dir.mkdir()

        allowed_taxids = None
        if args.taxid:
            if not args.nodes:
                log("ERROR: --nodes nodes.dmp required when --taxid is used")
                sys.exit(1)
            allowed_taxids = load_ncbi_descendants(args.nodes, args.taxid)

        abs_paths    = [str(Path(p).resolve()) for p in args.fasta]
        all_gz       = all(p.endswith('.gz') for p in abs_paths)
        force_serial = (args.workers == 1)

        try:
            if all_gz or force_serial:
                total_prot, total_pep = _build_serial(args, abs_paths,
                                                       allowed_taxids, tmp_dir)
            else:
                total_prot, total_pep = _build_parallel(args, abs_paths,
                                                         allowed_taxids, tmp_dir)
        except Exception:
            log("ERROR during build:")
            log(traceback.format_exc())
            log("Temp files kept: " + str(tmp_dir))
            sys.exit(1)

    log(f"Digestion done. {total_prot:,} proteins | {total_pep:,} peptides")

    try:
        log("Merging, assigning accids, sorting with DuckDB ...")
        if db_path.exists():
            db_path.unlink()

        con = duckdb.connect(str(db_path))
        con.execute(f"PRAGMA threads={args.threads}")
        con.execute(f"PRAGMA memory_limit='{args.memory}'")
        # Halves sort temp-space by not preserving insertion order.
        # Safe here: we ORDER BY mass_int anyway, so insertion order is irrelevant.
        con.execute("SET preserve_insertion_order=false")
        if args.tmp_dir:
            ddb_tmp = str(Path(args.tmp_dir) / 'duckdb_tmp').replace('\\', '/')
            con.execute(f"PRAGMA temp_directory='{ddb_tmp}'")

        # Accessions: deduplicate from the small acc_*.parquet files.
        # Indexed on acc string; no integer accid needed because peptides
        # already carry acc directly, so searches never need a join.
        acc_glob = str(tmp_dir / 'acc_*.parquet').replace('\\', '/')
        log("  Building accessions table ...")
        con.execute(f"""
            CREATE TABLE accessions AS
            SELECT acc, MAX(taxid) AS taxid
            FROM read_parquet('{acc_glob}')
            GROUP BY acc
        """)
        n_acc = con.execute("SELECT COUNT(*) FROM accessions").fetchone()[0]
        con.execute("CREATE INDEX acc_idx ON accessions(acc)")
        log(f"  accessions: {n_acc:,}")

        # Sort peptides into mass-range bands, each small enough to fit in RAM.
        # No global sort, no temp-disk spill.
        data_dir  = db_path.parent / (db_path.stem + '_data')
        pep_glob  = str(tmp_dir / 'pep_*.parquet').replace('\\', '/')
        data_glob = _sort_bands(con, pep_glob, data_dir, args.bands)

        # Worker Parquets are no longer needed; reclaim the space now.
        shutil.rmtree(tmp_dir)
        log("  Worker Parquet files deleted.")

        # VIEW over sorted band files; the band files ARE the data store.
        # Keep data_dir alongside the .duckdb file when moving or copying.
        con.execute(
            f"CREATE VIEW peptides AS SELECT * FROM read_parquet('{data_glob}')")
        n_pep = con.execute("SELECT COUNT(*) FROM peptides").fetchone()[0]
        log(f"  peptides:   {n_pep:,}")
        con.close()

    except Exception:
        log("ERROR during DuckDB merge/sort:")
        log(traceback.format_exc())
        if tmp_dir.exists():
            log("Temp files kept: " + str(tmp_dir))
        sys.exit(1)

    log(f"Done.  {db_path}")
    log(f"  Data dir  : {data_dir}  -- keep next to the .duckdb file")
    log(f"  proteins:   {total_prot:,}")
    log(f"  peptides:   {n_pep:,}")
    log(f"  accessions: {n_acc:,}")

# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    global _log_fh
    p = argparse.ArgumentParser(
        description='Build peptide DuckDB from UniProt FASTA/FASTA.gz')
    p.add_argument('--fasta',    required=True, nargs='+',
                   help='One or more FASTA or FASTA.gz files')
    p.add_argument('--db',       required=True)
    p.add_argument('--enzyme',   default='trypsin',
                   choices=['trypsin', 'chymotrypsin'])
    p.add_argument('--min',      type=int, default=7)
    p.add_argument('--max',      type=int, default=30)
    p.add_argument('--mc',       type=int, default=1)
    p.add_argument('--memory',   default='8GB',
                   help='DuckDB memory limit for sort step (default 8GB)')
    p.add_argument('--threads',  type=int,
                   default=max(1, cpu_count() or 1),
                   help='DuckDB threads for sort step (default: all CPUs)')
    p.add_argument('--workers',  type=int, default=None,
                   help='Parallel workers for digestion (default: cpu_count). '
                        'Set to 1 to force serial mode (required for .gz files).')
    p.add_argument('--log',      default=None,
                   help='Also write progress to this log file')
    p.add_argument('--bands',     type=int, default=20,
                   help='Mass-range bands for in-RAM sort (default 20). '
                        'Each band is sorted without disk spill. '
                        'Increase if you still get OOM during sort.')
    p.add_argument('--tmp-dir',   default=None,
                   help='Directory for temp Parquet files during build. '
                        'Use a fast LOCAL path (e.g. C:\\temp) when the '
                        'FASTA and output DB are on a network drive.')
    p.add_argument('--merge-only', action='store_true',
                   help='Skip digestion; clean stale Parquet files and run '
                        'only the DuckDB merge/sort step on an existing temp dir.')
    p.add_argument('--taxid',    type=int, nargs='+')
    p.add_argument('--nodes',    default=None)
    args = p.parse_args()

    if args.log:
        _log_fh = open(args.log, 'w', buffering=1)

    log("build_db.py starting")
    log(f"  fasta  : {args.fasta}")
    log(f"  db     : {args.db}")
    log(f"  memory : {args.memory}   threads: {args.threads}")

    build(args)

    if _log_fh:
        _log_fh.close()


if __name__ == '__main__':
    main()