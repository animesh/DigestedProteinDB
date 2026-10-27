"""
build_db.py - Build a peptide database from UniProt FASTA.

Design:
  1. Workers digest proteins and write peptides directly to
     mass-band Parquet files (band = mass range bucket).
  2. After digestion, each band is sorted in-memory with PyArrow.
     Each band is ~1 GB so it always fits in RAM. No temp disk.
  3. DuckDB VIEW over sorted band files for fast range queries.
     PSMs and accessions go in meta.duckdb.

No RocksDB. No Rust. No global sort. Dependencies: duckdb pyarrow.

Usage:
  python build_db.py ^
      --fasta uniprot_trembl.fasta ^
      --db    C:\\tmp\\trembl ^
      [--workers 12] [--taxid 2] [--nodes nodes.dmp] [--names names.dmp]
"""

import argparse, gzip, os, re, shutil, sys, traceback
from collections import defaultdict, deque
from datetime import datetime
from multiprocessing import Pool, cpu_count
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
import pyarrow.compute as pc

from common import digest, neutral_mass, pack_entry, MASS_INT_SCALE

# ── Mass band config ──────────────────────────────────────────────────────────
MASS_LO   = 4_000_000      # 400 Da × 1e4
MASS_HI   = 50_000_000     # 5000 Da × 1e4
N_BANDS   = 20
BAND_W    = (MASS_HI - MASS_LO) // N_BANDS
BAND_FLUSH = 100_000        # rows per band buffer; keep per-worker RAM low

def _band(mi: int) -> int:
    return max(0, min(N_BANDS - 1, (mi - MASS_LO) // BAND_W))

# ── Schemas ───────────────────────────────────────────────────────────────────
_PEP_SCHEMA = pa.schema([
    ('mass_int', pa.int32()), ('seq', pa.string()),
    ('acc',      pa.string()), ('taxid', pa.int32()),
])
_ACC_SCHEMA = pa.schema([('acc', pa.string()), ('taxid', pa.int32())])

# ── Logging ───────────────────────────────────────────────────────────────────
_log_fh = None
def log(msg):
    ts = datetime.now().strftime('%H:%M:%S')
    line = f"[{ts}] {msg}"
    print(line, file=sys.stderr, flush=True)
    if _log_fh: print(line, file=_log_fh, flush=True)

# ── FASTA parsing ─────────────────────────────────────────────────────────────
_OX_RE = re.compile(r'OX=(\d+)')

def _parse_header(line):
    fields = line.split('|')
    acc = fields[1] if len(fields) >= 2 else line.split()[0]
    m = _OX_RE.search(line)
    return acc, int(m.group(1)) if m else 0

def _iter_fasta_range(path, start, end):
    acc, taxid, parts = None, 0, []
    with open(path, 'r', encoding='ascii', errors='ignore') as fh:
        fh.seek(start)
        while True:
            lp = fh.tell(); line = fh.readline()
            if not line:
                if acc and parts: yield acc, taxid, ''.join(parts)
                break
            line = line.rstrip()
            if not line: continue
            if line.startswith('>'):
                if acc and parts: yield acc, taxid, ''.join(parts)
                if lp >= end: break
                acc, taxid = _parse_header(line[1:]); parts = []
            elif acc is not None:
                parts.append(line.upper())

def _iter_fasta_gz(path):
    fn = gzip.open if path.endswith('.gz') else open
    acc, taxid, parts = None, 0, []
    with fn(path, 'rt', encoding='ascii', errors='ignore') as fh:
        for line in fh:
            line = line.rstrip()
            if not line: continue
            if line.startswith('>'):
                if acc and parts: yield acc, taxid, ''.join(parts)
                acc, taxid = _parse_header(line[1:]); parts = []
            elif acc is not None: parts.append(line.upper())
    if acc and parts: yield acc, taxid, ''.join(parts)

def _chunk_boundaries(path, n):
    size = os.path.getsize(path); bounds = [0]
    with open(path, 'rb') as fh:
        for i in range(1, n):
            t = (size * i) // n; fh.seek(max(0, t-1))
            buf = fh.read(1 << 20); idx = buf.find(b'\n>')
            bounds.append((max(0,t-1)+idx+1) if idx >= 0 else size)
    bounds.append(size); return bounds

# ── NCBI helpers ──────────────────────────────────────────────────────────────
def _descendants(nodes_dmp, parents):
    log("Loading taxonomy ..."); ch = defaultdict(list)
    with open(nodes_dmp) as fh:
        for line in fh:
            p = line.split('|')
            if len(p) < 2: continue
            try: ch[int(p[1].strip())].append(int(p[0].strip()))
            except ValueError: pass
    allowed, q = set(), deque(parents)
    while q:
        t = q.popleft()
        if t in allowed: continue
        allowed.add(t); q.extend(ch.get(t, []))
    log(f"  {len(allowed):,} allowed taxids."); return allowed

def _names(names_dmp):
    log("Loading species names ..."); names = {}
    with open(names_dmp) as fh:
        for line in fh:
            p = [x.strip() for x in line.split('|')]
            if len(p) >= 4 and p[3] == 'scientific name':
                try: names[int(p[0])] = p[1]
                except ValueError: pass
    log(f"  {len(names):,} names."); return names

# ── Worker (pure Python + pyarrow, no Rust) ───────────────────────────────────
def _flush_band(b, bufs, bchk, tmp_dir, worker_id):
    ms, seqs, accs, tids = bufs[b]
    if not ms: return
    pq.write_table(
        pa.table({'mass_int': pa.array(ms, pa.int32()),
                  'seq':      pa.array(seqs, pa.string()),
                  'acc':      pa.array(accs, pa.string()),
                  'taxid':    pa.array(tids, pa.int32())},
                 schema=_PEP_SCHEMA),
        str(Path(tmp_dir) / f'wband_{b:03d}_{worker_id:04d}_{bchk[b]:05d}.parquet'),
        compression='zstd', write_statistics=True)
    bufs[b] = ([], [], [], []); bchk[b] += 1

def _worker_fn(args):
    path, start, end, allowed, min_l, max_l, mc, enzyme, tmp_dir, wid = args
    bufs = [([], [], [], []) for _ in range(N_BANDS)]
    bchk = [0] * N_BANDS; aacc = []; atid = []; n_prot = n_pep = 0
    for acc, taxid, seq in _iter_fasta_range(path, start, end):
        if allowed is not None and taxid not in allowed: continue
        aacc.append(acc); atid.append(taxid)
        peps = digest(seq, min_l, max_l, mc, enzyme)
        for pep in peps:
            mi = round(neutral_mass(pep) * MASS_INT_SCALE)
            b = _band(mi); ms, seqs, accs, tids = bufs[b]
            ms.append(mi); seqs.append(pep)
            accs.append(acc); tids.append(taxid)
        n_prot += 1; n_pep += len(peps)
        for b in range(N_BANDS):
            if len(bufs[b][0]) >= BAND_FLUSH:
                _flush_band(b, bufs, bchk, tmp_dir, wid)
    for b in range(N_BANDS): _flush_band(b, bufs, bchk, tmp_dir, wid)
    pq.write_table(
        pa.table({'acc': pa.array(aacc), 'taxid': pa.array(atid, pa.int32())},
                 schema=_ACC_SCHEMA),
        str(Path(tmp_dir) / f'acc_{wid:04d}.parquet'), compression='zstd')
    return wid, n_prot, n_pep

# ── In-memory band sort (PyArrow) ─────────────────────────────────────────────
_LARGE_SCHEMA = pa.schema([
    pa.field('mass_int', pa.int32()),
    pa.field('seq',      pa.large_string()),
    pa.field('acc',      pa.large_string()),
    pa.field('taxid',    pa.int32()),
])

# PyArrow sort needs ~3x data size in RAM (input + indices + output).
# At 36 GB machine with 8 GB OS overhead: safe up to ~9 GB = ~270M rows.
_PYARROW_ROW_LIMIT = 250_000_000


def _sort_bands(tmp_dir: Path, data_dir: Path, memory: str = '28GB') -> str:
    """
    Sort each band's files.
    Small bands (≤250M rows): PyArrow in-memory – fastest.
    Large bands (>250M rows): DuckDB COPY ORDER BY – memory-efficient,
      spills to tmp_dir if needed, never crashes.
    """
    data_dir.mkdir(exist_ok=True)
    for b in range(N_BANDS):
        files = sorted(tmp_dir.glob(f'wband_{b:03d}_*.parquet'))
        if not files:
            continue
        total_rows = sum(pq.read_metadata(str(f)).num_rows for f in files)
        out = str(data_dir / f'peptides_{b:05d}.parquet')
        log(f"  Band {b+1}/{N_BANDS}: {len(files)} file(s), {total_rows:,} rows ...")

        if total_rows <= _PYARROW_ROW_LIMIT:
            # In-memory PyArrow sort
            table = pq.read_table([str(f) for f in files], schema=_PEP_SCHEMA)
            table = table.cast(_LARGE_SCHEMA)
            idx   = pc.sort_indices(table, sort_keys=[('mass_int', 'ascending')])
            pq.write_table(table.take(idx), out,
                           compression='zstd', row_group_size=500_000,
                           write_statistics=True)
        else:
            # DuckDB sort: memory-efficient, spills to tmp_dir if needed
            file_list = ', '.join(f"'{f}'" for f in files)
            ddb_tmp   = str(tmp_dir / 'ddb_tmp').replace('\\', '/')
            con = duckdb.connect()
            con.execute(f"SET memory_limit='{memory}'")
            con.execute(f"SET temp_directory='{ddb_tmp}'")
            con.execute("SET preserve_insertion_order=false")
            con.execute(f"""
                COPY (SELECT * FROM read_parquet([{file_list}]) ORDER BY mass_int)
                TO '{out}' (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 500000)
            """)
            con.close()

        for f in files:
            f.unlink()

    return str(data_dir / 'peptides_*.parquet').replace('\\', '/')

# ── Build ─────────────────────────────────────────────────────────────────────
def build(args):
    db_dir  = Path(args.db); db_dir.mkdir(parents=True, exist_ok=True)
    data_dir = db_dir / 'peptides_data'
    meta_p  = db_dir / 'meta.duckdb'
    tmp_dir = (Path(args.tmp) if args.tmp else db_dir.parent) / (db_dir.name + '_tmp')
    log(f"DB  : {db_dir}"); log(f"Tmp : {tmp_dir}")

    if tmp_dir.exists(): shutil.rmtree(tmp_dir)
    tmp_dir.mkdir()

    allowed = _descendants(args.nodes, args.taxid) if args.taxid else None
    if args.taxid and not args.nodes:
        log("ERROR: --nodes required with --taxid"); sys.exit(1)

    abs_paths = [str(Path(p).resolve()) for p in args.fasta]
    all_gz    = all(p.endswith('.gz') for p in abs_paths)
    nw        = args.workers or max(1, cpu_count())

    # ── digestion ─────────────────────────────────────────────────────────────
    if all_gz or nw == 1:
        log("Workers: 1 (serial)")
        bufs = [([], [], [], []) for _ in range(N_BANDS)]
        bchk = [0] * N_BANDS; aacc = []; atid = []; np = npe = 0
        for path in abs_paths:
            log(f"Parsing {path} ...")
            for acc, taxid, seq in _iter_fasta_gz(path):
                if allowed and taxid not in allowed: continue
                aacc.append(acc); atid.append(taxid)
                peps = digest(seq, args.min, args.max, args.mc, args.enzyme)
                for pep in peps:
                    mi = round(neutral_mass(pep) * MASS_INT_SCALE); b = _band(mi)
                    ms, seqs, accs, tids = bufs[b]
                    ms.append(mi); seqs.append(pep); accs.append(acc); tids.append(taxid)
                np += 1; npe += len(peps)
                if np % 100_000 == 0: log(f"  {np:,} proteins | {npe:,} peptides")
                for b in range(N_BANDS):
                    if len(bufs[b][0]) >= BAND_FLUSH:
                        _flush_band(b, bufs, bchk, tmp_dir, 0)
        for b in range(N_BANDS): _flush_band(b, bufs, bchk, tmp_dir, 0)
        pq.write_table(
            pa.table({'acc': pa.array(aacc), 'taxid': pa.array(atid, pa.int32())},
                     schema=_ACC_SCHEMA),
            str(tmp_dir / 'acc_0000.parquet'), compression='zstd')
        total_p, total_pe = np, npe
    else:
        log(f"Workers: {nw} (parallel)")
        jobs = []
        for path in abs_paths:
            log(f"  Splitting {path} ...")
            bounds = _chunk_boundaries(path, nw)
            for s, e in zip(bounds, bounds[1:]):
                jobs.append((path, s, e, allowed,
                             args.min, args.max, args.mc, args.enzyme,
                             str(tmp_dir), len(jobs)))
        total_p = total_pe = 0
        try:
            with Pool(nw) as pool:
                for wid, np, npe in pool.imap_unordered(_worker_fn, jobs, chunksize=1):
                    total_p += np; total_pe += npe
                    log(f"  Worker {wid:3d}: {np:,} proteins | {npe:,} peptides")
        except KeyboardInterrupt: log("Interrupted."); sys.exit(1)
        except Exception:
            log("ERROR:"); log(traceback.format_exc())
            log(f"Temp kept: {tmp_dir}"); sys.exit(1)

    log(f"Digestion done. {total_p:,} proteins | {total_pe:,} peptides")

    # ── sort ──────────────────────────────────────────────────────────────────
    log("Sorting bands in-memory (PyArrow) ...")
    if data_dir.exists(): shutil.rmtree(data_dir)
    try:
        data_glob = _sort_bands(tmp_dir, data_dir)
    except Exception:
        log("ERROR:"); log(traceback.format_exc()); log(f"Temp: {tmp_dir}"); sys.exit(1)

    # ── meta.duckdb ───────────────────────────────────────────────────────────
    log("Building meta.duckdb ...")
    smap = _names(args.names) if args.names else {}
    if meta_p.exists(): meta_p.unlink()
    con = duckdb.connect(str(meta_p))
    acc_glob = str(tmp_dir / 'acc_*.parquet').replace('\\', '/')
    con.execute(f"""
        CREATE TABLE accessions AS
        SELECT acc, MAX(taxid) AS taxid, '' AS species
        FROM read_parquet('{acc_glob}') GROUP BY acc""")
    if smap:
        for tid, name in smap.items():
            con.execute("UPDATE accessions SET species=? WHERE taxid=?", [name, tid])
    con.execute("CREATE INDEX acc_idx ON accessions(acc)")
    n_acc = con.execute("SELECT COUNT(*) FROM accessions").fetchone()[0]
    log(f"  accessions: {n_acc:,}")

    dg = str(data_dir / 'peptides_*.parquet').replace('\\', '/')
    con.execute(f"CREATE VIEW peptides AS SELECT * FROM read_parquet('{dg}')")
    n_pep = con.execute("SELECT COUNT(*) FROM peptides").fetchone()[0]
    log(f"  peptides:   {n_pep:,}")

    con.execute("""CREATE TABLE psms (
        run TEXT, scan_id INTEGER, title TEXT,
        precursor_da DOUBLE, charge INTEGER,
        seq TEXT, acc TEXT, taxid INTEGER,
        ppm_err DOUBLE, ptm TEXT, ptm_pos INTEGER,
        n_ions INTEGER, frag_score DOUBLE,
        created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP)""")
    con.close()

    shutil.rmtree(tmp_dir)
    log(f"\nDone.  {db_dir}/")
    log(f"  peptides    : {n_pep:,}")
    log(f"  accessions  : {n_acc:,}")

# ── CLI ───────────────────────────────────────────────────────────────────────
def main():
    global _log_fh
    p = argparse.ArgumentParser()
    p.add_argument('--fasta',   required=True, nargs='+')
    p.add_argument('--db',      required=True)
    p.add_argument('--tmp',     default=None)
    p.add_argument('--workers', type=int, default=None)
    p.add_argument('--enzyme',  default='trypsin', choices=['trypsin','chymotrypsin'])
    p.add_argument('--min',     type=int, default=7)
    p.add_argument('--max',     type=int, default=30)
    p.add_argument('--mc',      type=int, default=1)
    p.add_argument('--taxid',   type=int, nargs='+')
    p.add_argument('--nodes',   default=None)
    p.add_argument('--names',   default=None)
    p.add_argument('--log',     default=None)
    args = p.parse_args()
    if args.log: _log_fh = open(args.log, 'w', buffering=1)
    log("build_db.py"); log(f"  fasta: {args.fasta}"); log(f"  db: {args.db}")
    build(args)
    if _log_fh: _log_fh.close()

if __name__ == '__main__':
    main()