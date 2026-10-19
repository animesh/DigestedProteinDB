"""
query.py - Query the peptide database interactively.

Examples:
  # All peptides near 1500.6 Da
  python query.py --db combined --mass 1500.6 --ppm 20

  # From m/z
  python query.py --db combined --mz 501.2 --charge 3 --ppm 20

  # Filter to human only
  python query.py --db combined --mass 1500.6 --taxid 9606

  # Include PSMs from prior searches
  python query.py --db combined --mass 1500.6 --show-psms

  # Species distribution across all stored PSMs
  python query.py --db combined --species-dist

  # Species distribution for one run
  python query.py --db combined --species-dist --run sample1

  # Export all PSMs
  python query.py --db combined --export-psms results.tsv
"""

import argparse, sys
from pathlib import Path

import duckdb
from rocksdict import Rdict, Options

from common import (PROTON, mass_key, key_to_mass, unpack_entries,
                    ppm_bounds, ppm_error, MASS_INT_SCALE)

def _open_rocks(db_dir):
    opts = Options(raw_mode=True)
    opts.create_if_missing(False)
    return Rdict(str(Path(db_dir) / 'peptides.rocks'), opts)

def _open_meta(db_dir):
    return duckdb.connect(str(Path(db_dir) / 'meta.duckdb'), read_only=True)

# ── Mass lookup ───────────────────────────────────────────────────────────────

def lookup(args):
    if args.mz and args.charge:
        neutral = args.mz * args.charge - args.charge * PROTON
    elif args.mass:
        neutral = args.mass
    else:
        sys.exit("Provide --mass or --mz + --charge")

    lo, hi = ppm_bounds(neutral, args.ppm)
    taxid_set = set(args.taxid) if args.taxid else None

    rocks = _open_rocks(args.db)
    meta  = _open_meta(args.db)

    species_map = dict(meta.execute(
        "SELECT taxid, species FROM accessions WHERE species != ''").fetchall())

    # Collect hits grouped by sequence
    hits = {}   # seq -> {mass_int, ppm, [(acc, taxid, species)]}
    hi_key = mass_key(hi)
    for k, v in rocks.items(from_key=mass_key(lo)):
        if k > hi_key: break
        mi = key_to_mass(k)
        ppm = ppm_error(neutral, mi / MASS_INT_SCALE)
        for seq, acc, taxid in unpack_entries(v):
            if taxid_set and taxid not in taxid_set:
                continue
            if seq not in hits:
                hits[seq] = {'mass_int': mi, 'ppm': ppm, 'proteins': []}
            hits[seq]['proteins'].append(
                (acc, taxid, species_map.get(taxid, f'taxid:{taxid}')))
    rocks.close()

    if not hits:
        print(f"No hits for mass {neutral:.5f} Da ± {args.ppm} ppm"); return

    print(f"\nQuery: {neutral:.5f} Da ± {args.ppm} ppm  →  {len(hits)} peptide(s)\n")

    for seq, info in sorted(hits.items(), key=lambda x: abs(x[1]['ppm'])):
        print(f"  {seq}  ({info['mass_int']/MASS_INT_SCALE:.5f} Da,"
              f" {info['ppm']:+.2f} ppm)")
        for acc, taxid, sp in sorted(set(info['proteins'])):
            print(f"    {acc}  {sp}  [taxid:{taxid}]")
        if args.show_psms:
            rows = meta.execute("""
                SELECT run, COUNT(*) AS n, MIN(frag_score) AS min_score,
                       MAX(frag_score) AS max_score
                FROM psms WHERE seq = ?
                GROUP BY run ORDER BY run
            """, [seq]).fetchall()
            if rows:
                print(f"    PSMs:")
                for run, n, lo_s, hi_s in rows:
                    print(f"      {run}: {n} PSMs  score [{lo_s:.3f}–{hi_s:.3f}]")
            else:
                print(f"    PSMs: none stored")
        print()

    meta.close()

# ── Species distribution ──────────────────────────────────────────────────────

def species_dist(args):
    meta = _open_meta(args.db)
    where = "WHERE p.run = ?" if args.run else ""
    params = [args.run] if args.run else []
    rows = meta.execute(f"""
        SELECT COALESCE(a.species, 'taxid:' || CAST(p.taxid AS TEXT)) AS species,
               COUNT(*) AS n_psms
        FROM psms p
        LEFT JOIN accessions a ON a.acc = p.acc
        {where}
        GROUP BY 1
        ORDER BY n_psms DESC
    """, params).fetchall()
    meta.close()

    if not rows:
        print("No PSMs stored yet. Run search_mgf.py first."); return

    title = f"run='{args.run}'" if args.run else "all runs"
    total = sum(r[1] for r in rows)
    print(f"\nSpecies distribution ({title}, {total:,} PSMs total)\n")
    for sp, n in rows[:50]:
        bar = '█' * int(30 * n / total)
        print(f"  {n:6,}  {n/total*100:5.1f}%  {bar}  {sp}")
    if len(rows) > 50:
        rest = sum(r[1] for r in rows[50:])
        print(f"  {rest:6,}  {rest/total*100:5.1f}%  ... ({len(rows)-50} more species)")

# ── Export PSMs ───────────────────────────────────────────────────────────────

def export_psms(args):
    meta = duckdb.connect(str(Path(args.db) / 'meta.duckdb'), read_only=True)
    where = "WHERE p.run = ?" if args.run else ""
    params = [args.run] if args.run else []
    rows = meta.execute(f"""
        SELECT p.run, p.scan_id, p.title, p.precursor_da, p.charge,
               p.seq, p.acc, p.taxid,
               COALESCE(a.species,'') AS species,
               p.ppm_err, p.ptm, p.ptm_pos, p.n_ions, p.frag_score
        FROM psms p
        LEFT JOIN accessions a ON a.acc = p.acc
        {where}
        ORDER BY p.run, p.scan_id
    """, params).fetchall()
    meta.close()

    out = args.export_psms
    cols = ['run','scan_id','title','precursor_da','charge','seq','acc',
            'taxid','species','ppm_err','ptm','ptm_pos','n_ions','frag_score']
    with open(out, 'w') as fh:
        print('\t'.join(cols), file=fh)
        for row in rows:
            print('\t'.join(str(v) for v in row), file=fh)
    print(f"Exported {len(rows):,} PSMs → {out}")

# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(description='Query the peptide database')
    p.add_argument('--db',           required=True)
    # Mass lookup
    p.add_argument('--mass',         type=float, default=None,
                   help='Neutral monoisotopic mass in Da')
    p.add_argument('--mz',           type=float, default=None,
                   help='Precursor m/z')
    p.add_argument('--charge',       type=int,   default=None)
    p.add_argument('--ppm',          type=float, default=20.0)
    p.add_argument('--taxid',        type=int,   nargs='+',
                   help='Filter to these NCBI taxids')
    p.add_argument('--show-psms',    action='store_true',
                   help='Show stored PSMs for each matching peptide')
    # Analytics
    p.add_argument('--species-dist', action='store_true',
                   help='Print species distribution of stored PSMs')
    p.add_argument('--run',          default=None,
                   help='Filter analytics to this run label')
    p.add_argument('--export-psms',  default=None, metavar='FILE',
                   help='Export all PSMs to TSV')
    args = p.parse_args()

    if args.species_dist:
        species_dist(args)
    elif args.export_psms:
        export_psms(args)
    elif args.mass or (args.mz and args.charge):
        lookup(args)
    else:
        p.print_help()

if __name__ == '__main__':
    main()