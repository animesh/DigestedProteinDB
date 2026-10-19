"""
build_db.py - Build a peptide database from UniProt FASTA.

Output: <db>/ directory containing
  peptides.rocks/   RocksDB: 4-byte big-endian mass_int → packed (seq,acc,taxid)
  meta.duckdb       accessions (acc, taxid, species) + empty psms table

Workers use only Python stdlib (no rocksdict, no Rust).  Each worker
writes sorted binary chunks: [4B mass_int BE][4B blob_len BE][blob].
The main process does a k-way merge of all chunks into the final RocksDB.

Usage:
  python build_db.py \\
      --fasta  uniprot_sprot.fasta uniprot_trembl.fasta \\
      --db     C:\\proteomics\\combined \\
      [--tmp   C:\\tmp] [--workers 12] \\
      [--enzyme trypsin] [--min 7] [--max 30] [--mc 1] \\
      [--taxid 2] [--nodes nodes.dmp] [--names names.dmp] \\
      [--log   build.log]

  pip install rocksdict duckdb pyarrow
"""

import argparse, gzip, heapq, os, re, shutil, struct, sys, traceback
from collections import defaultdict, deque
from datetime import datetime
from multiprocessing import Pool, cpu_count
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq

from common import (digest, neutral_mass, pack_entry, mass_key, MASS_INT_SCALE)

# ── Logging ───────────────────────────────────────────────────────────────────

_log_fh = None

def log(msg):
    ts = datetime.now().strftime('%H:%M:%S')
    line = f"[{ts}] {msg}"
    print(line, file=sys.stderr, flush=True)
    if _log_fh:
        print(line, file=_log_fh, flush=True)

# ── Binary chunk I/O ──────────────────────────────────────────────────────────
# Format: [mass_int: 4B BE][blob_len: 4B BE][blob: blob_len B]  (sorted by mass_int)

_HDR = struct.Struct('>II')   # big-endian unsigned int, unsigned int

def _write_chunk(path: str, buf: dict):
    """Write sorted {mass_int_int: bytearray} to a binary chunk file."""
    with open(path, 'wb') as f:
        for mi in sorted(buf):
            blob = bytes(buf[mi])
            f.write(_HDR.pack(mi, len(blob)))
            f.write(blob)

def _iter_chunk(path: str):
    """Yield (mass_key_bytes, blob_bytes) from a sorted binary chunk."""
    hdr_size = _HDR.size
    with open(path, 'rb') as f:
        while True:
            hdr = f.read(hdr_size)
            if len(hdr) < hdr_size:
                break
            mi, blen = _HDR.unpack(hdr)
            blob = f.read(blen)
            yield mass_key(mi), blob

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
            line_pos = fh.tell()
            line = fh.readline()
            if not line:
                if acc and parts: yield acc, taxid, ''.join(parts)
                break
            line = line.rstrip()
            if not line: continue
            if line.startswith('>'):
                if acc and parts: yield acc, taxid, ''.join(parts)
                if line_pos >= end: break
                acc, taxid = _parse_header(line[1:])
                parts = []
            elif acc is not None:
                parts.append(line.upper())

def _iter_fasta_gz(path):
    open_fn = gzip.open if path.endswith('.gz') else open
    acc, taxid, parts = None, 0, []
    with open_fn(path, 'rt', encoding='ascii', errors='ignore') as fh:
        for line in fh:
            line = line.rstrip()
            if not line: continue
            if line.startswith('>'):
                if acc and parts: yield acc, taxid, ''.join(parts)
                acc, taxid = _parse_header(line[1:])
                parts = []
            elif acc is not None:
                parts.append(line.upper())
    if acc and parts: yield acc, taxid, ''.join(parts)

def _find_chunk_boundaries(path, n):
    size = os.path.getsize(path)
    boundaries = [0]
    with open(path, 'rb') as fh:
        for i in range(1, n):
            target = (size * i) // n
            fh.seek(max(0, target - 1))
            buf = fh.read(1 << 20)
            idx = buf.find(b'\n>')
            boundaries.append((max(0, target-1)+idx+1) if idx >= 0 else size)
    boundaries.append(size)
    return boundaries

# ── NCBI taxonomy ─────────────────────────────────────────────────────────────

def load_ncbi_descendants(nodes_dmp, parent_taxids):
    log("Loading taxonomy ...")
    children = defaultdict(list)
    with open(nodes_dmp) as fh:
        for line in fh:
            parts = line.split('|')
            if len(parts) < 2: continue
            try: children[int(parts[1].strip())].append(int(parts[0].strip()))
            except ValueError: pass
    allowed, q = set(), deque(parent_taxids)
    while q:
        t = q.popleft()
        if t in allowed: continue
        allowed.add(t); q.extend(children.get(t, []))
    log(f"  {len(allowed):,} allowed taxids.")
    return allowed

def load_ncbi_names(names_dmp):
    log("Loading species names ...")
    names = {}
    with open(names_dmp) as fh:
        for line in fh:
            parts = [p.strip() for p in line.split('|')]
            if len(parts) >= 4 and parts[3] == 'scientific name':
                try: names[int(parts[0])] = parts[1]
                except ValueError: pass
    log(f"  {len(names):,} names loaded.")
    return names

# ── Worker (no rocksdict - stdlib only) ───────────────────────────────────────

FLUSH_KEYS = 500_000   # unique mass_int keys per binary chunk

def _worker_fn(args):
    """
    Digest proteins and write sorted binary chunk files.
    No rocksdict here: avoids Rust/spawn crash on Windows.
    """
    (path, start, end, allowed_taxids,
     min_len, max_len, mc, enzyme, tmp_dir, worker_id) = args

    buf = {}        # mass_int (int) → bytearray of packed entries
    aacc, atid = [], []
    n_prot = n_pep = chunk_id = 0

    def flush():
        nonlocal chunk_id
        out = str(Path(tmp_dir) / f'w{worker_id:04d}_{chunk_id:05d}.bin')
        _write_chunk(out, buf)
        buf.clear()
        chunk_id += 1

    for acc, taxid, seq in _iter_fasta_range(path, start, end):
        if allowed_taxids is not None and taxid not in allowed_taxids:
            continue
        aacc.append(acc); atid.append(taxid)
        peps = digest(seq, min_len, max_len, mc, enzyme)
        for pep in peps:
            mi = round(neutral_mass(pep) * MASS_INT_SCALE)
            entry = pack_entry(pep, acc, taxid)
            if mi in buf: buf[mi] += entry
            else: buf[mi] = bytearray(entry)
        n_prot += 1; n_pep += len(peps)
        if len(buf) >= FLUSH_KEYS:
            flush()

    if buf:
        flush()

    pq.write_table(
        pa.table({'acc': pa.array(aacc), 'taxid': pa.array(atid, pa.int32())}),
        str(Path(tmp_dir) / f'acc_{worker_id:04d}.parquet'),
        compression='zstd')

    return worker_id, n_prot, n_pep, chunk_id

# ── K-way merge into RocksDB (main process only) ─────────────────────────────

def _kway_merge_to_rocks(chunk_files, rocks_path):
    """
    K-way merge of sorted binary chunk files into a single RocksDB.
    rocksdict is imported here, in the main process only.
    """
    from rocksdict import Rdict, Options, WriteBatch, DBCompressionType

    opts = Options(raw_mode=True)
    opts.create_if_missing(True)
    opts.set_compression_type(DBCompressionType.lz4())
    opts.set_write_buffer_size(256 << 20)
    out = Rdict(str(rocks_path), opts)

    iters = [_iter_chunk(p) for p in chunk_files]
    heap  = []
    for i, it in enumerate(iters):
        try:
            k, v = next(it)
            heapq.heappush(heap, (k, i, v))
        except StopIteration:
            pass

    cur_key, cur_val = None, b''
    batch = WriteBatch(raw_mode=True)
    n_keys = 0

    while heap:
        k, i, v = heapq.heappop(heap)
        if k == cur_key:
            cur_val += v
        else:
            if cur_key is not None:
                batch.put(cur_key, cur_val)
                n_keys += 1
                if n_keys % 1_000_000 == 0:
                    out.write(batch)
                    batch = WriteBatch(raw_mode=True)
                    log(f"  Merged {n_keys:,} unique masses ...")
            cur_key, cur_val = k, v
        try:
            nk, nv = next(iters[i])
            heapq.heappush(heap, (nk, i, nv))
        except StopIteration:
            pass

    if cur_key:
        batch.put(cur_key, cur_val)
        n_keys += 1
    out.write(batch)
    out.close()
    return n_keys

# ── Build ─────────────────────────────────────────────────────────────────────

def build(args):
    db_dir     = Path(args.db)
    db_dir.mkdir(parents=True, exist_ok=True)
    rocks_path = db_dir / 'peptides.rocks'
    meta_path  = db_dir / 'meta.duckdb'
    tmp_root   = Path(args.tmp) if args.tmp else db_dir.parent
    tmp_dir    = tmp_root / (db_dir.name + '_build_tmp')

    log(f"DB        : {db_dir}")
    log(f"Temp      : {tmp_dir}")

    if tmp_dir.exists():
        log("Removing stale temp dir ...")
        shutil.rmtree(tmp_dir)
    tmp_dir.mkdir()

    allowed_taxids = None
    if args.taxid:
        if not args.nodes:
            log("ERROR: --nodes required with --taxid"); sys.exit(1)
        allowed_taxids = load_ncbi_descendants(args.nodes, args.taxid)

    abs_paths = [str(Path(p).resolve()) for p in args.fasta]
    all_gz    = all(p.endswith('.gz') for p in abs_paths)
    n_workers = args.workers or max(1, cpu_count())

    # ── Digestion ─────────────────────────────────────────────────────────────
    if all_gz or n_workers == 1:
        log("Workers: 1 (serial)")
        buf = {}; aacc, atid = [], []; n_prot = n_pep = chunk_id = 0

        def serial_flush():
            nonlocal chunk_id
            _write_chunk(str(tmp_dir / f'w0000_{chunk_id:05d}.bin'), buf)
            buf.clear(); chunk_id += 1  # noqa: E702

        for path in abs_paths:
            log(f"Parsing {path} ...")
            for acc, taxid, seq in _iter_fasta_gz(path):
                if allowed_taxids and taxid not in allowed_taxids: continue
                aacc.append(acc); atid.append(taxid)
                peps = digest(seq, args.min, args.max, args.mc, args.enzyme)
                for pep in peps:
                    mi = round(neutral_mass(pep) * MASS_INT_SCALE)
                    entry = pack_entry(pep, acc, taxid)
                    if mi in buf: buf[mi] += entry
                    else: buf[mi] = bytearray(entry)
                n_prot += 1; n_pep += len(peps)
                if n_prot % 100_000 == 0:
                    log(f"  {n_prot:,} proteins | {n_pep:,} peptides")
                if len(buf) >= FLUSH_KEYS:
                    serial_flush()
        if buf: serial_flush()
        pq.write_table(
            pa.table({'acc': pa.array(aacc), 'taxid': pa.array(atid, pa.int32())}),
            str(tmp_dir / 'acc_0000.parquet'), compression='zstd')
        total_prot, total_pep = n_prot, n_pep
    else:
        log(f"Workers: {n_workers} (parallel byte-range)")
        jobs = []
        for path in abs_paths:
            log(f"  Splitting {path} ...")
            bounds = _find_chunk_boundaries(path, n_workers)
            for s, e in zip(bounds, bounds[1:]):
                wid = len(jobs)
                jobs.append((path, s, e, allowed_taxids,
                             args.min, args.max, args.mc, args.enzyme,
                             str(tmp_dir), wid))
        total_prot = total_pep = 0
        try:
            with Pool(n_workers) as pool:
                for wid, np, npe, nc in pool.imap_unordered(
                        _worker_fn, jobs, chunksize=1):
                    total_prot += np; total_pep += npe
                    log(f"  Worker {wid:3d}: {np:,} proteins | {npe:,} peptides"
                        f" | {nc} chunk(s)")
        except KeyboardInterrupt:
            log("Interrupted."); sys.exit(1)
        except Exception:
            log("ERROR:"); log(traceback.format_exc())
            log("Temp kept: " + str(tmp_dir)); sys.exit(1)

    log(f"Digestion done. {total_prot:,} proteins | {total_pep:,} peptides")

    # ── K-way merge ───────────────────────────────────────────────────────────
    chunk_files = sorted(str(p) for p in tmp_dir.glob('w*.bin'))
    log(f"K-way merging {len(chunk_files)} sorted chunks into RocksDB ...")
    if rocks_path.exists():
        shutil.rmtree(str(rocks_path))
    try:
        n_keys = _kway_merge_to_rocks(chunk_files, rocks_path)
    except Exception:
        log("ERROR:"); log(traceback.format_exc())
        log("Temp kept: " + str(tmp_dir)); sys.exit(1)
    log(f"  {n_keys:,} unique mass keys.")

    # ── meta.duckdb ───────────────────────────────────────────────────────────
    log("Building meta.duckdb ...")
    species_map = load_ncbi_names(args.names) if args.names else {}
    if meta_path.exists():
        meta_path.unlink()
    con = duckdb.connect(str(meta_path))
    acc_glob = str(tmp_dir / 'acc_*.parquet').replace('\\', '/')
    con.execute(f"""
        CREATE TABLE accessions AS
        SELECT acc, MAX(taxid) AS taxid, '' AS species
        FROM read_parquet('{acc_glob}') GROUP BY acc
    """)
    if species_map:
        for taxid, name in species_map.items():
            con.execute("UPDATE accessions SET species=? WHERE taxid=?", [name, taxid])
    con.execute("CREATE INDEX acc_idx ON accessions(acc)")
    n_acc = con.execute("SELECT COUNT(*) FROM accessions").fetchone()[0]
    log(f"  accessions: {n_acc:,}")
    con.execute("""
        CREATE TABLE psms (
            run TEXT, scan_id INTEGER, title TEXT,
            precursor_da DOUBLE, charge INTEGER,
            seq TEXT, acc TEXT, taxid INTEGER,
            ppm_err DOUBLE, ptm TEXT, ptm_pos INTEGER,
            n_ions INTEGER, frag_score DOUBLE,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    """)
    con.close()

    shutil.rmtree(tmp_dir)
    log(f"\nDone.  {db_dir}/")
    log(f"  unique masses : {n_keys:,}")
    log(f"  accessions    : {n_acc:,}")

# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    global _log_fh
    p = argparse.ArgumentParser(description='Build peptide DB from UniProt FASTA')
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
    if args.log:
        _log_fh = open(args.log, 'w', buffering=1)
    log("build_db.py starting")
    log(f"  fasta  : {args.fasta}")
    log(f"  db     : {args.db}")
    build(args)
    if _log_fh: _log_fh.close()

if __name__ == '__main__':
    main()