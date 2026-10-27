"""
query.py - Query the peptide database.

Examples:
  python query.py --db C:\\tmp\\trembl --mass 1500.6 --ppm 20
  python query.py --db C:\\tmp\\trembl --mz 501.2 --charge 3
  python query.py --db C:\\tmp\\trembl --mass 1500.6 --taxid 9606
  python query.py --db C:\\tmp\\trembl --mass 1500.6 --show-psms
  python query.py --db C:\\tmp\\trembl --species-dist
  python query.py --db C:\\tmp\\trembl --species-dist --run sample1
  python query.py --db C:\\tmp\\trembl --export-psms all.tsv
"""

import argparse, sys
from pathlib import Path

import duckdb

from common import PROTON, ppm_bounds, ppm_error, precursor_neutral, MASS_INT_SCALE

def _meta(db, read_only=True):
    return duckdb.connect(str(Path(db) / 'meta.duckdb'), read_only=read_only)

# ── Mass lookup ───────────────────────────────────────────────────────────────

def lookup(args):
    if args.mz and args.charge:
        neutral = precursor_neutral(args.mz, args.charge)
    elif args.mass:
        neutral = args.mass
    else:
        sys.exit("Provide --mass or --mz + --charge")

    lo, hi    = ppm_bounds(neutral, args.ppm)
    taxid_set = set(args.taxid) if args.taxid else None
    con       = _meta(args.db, read_only=False)

    species_map = dict(con.execute(
        "SELECT taxid, species FROM accessions WHERE species != ''").fetchall())

    taxid_clause = f"AND taxid IN ({','.join(str(t) for t in taxid_set)})" \
                   if taxid_set else ""

    rows = con.execute(f"""
        SELECT mass_int, seq, acc, taxid
        FROM peptides
        WHERE mass_int >= {lo} AND mass_int <= {hi}
        {taxid_clause}
        ORDER BY mass_int
    """).fetchall()

    if not rows:
        print(f"No hits for {neutral:.5f} Da ± {args.ppm} ppm"); con.close(); return

    # Group by sequence
    hits = {}
    for mass_int, seq, acc, taxid in rows:
        ppm = ppm_error(neutral, mass_int / MASS_INT_SCALE)
        if seq not in hits:
            hits[seq] = {'mass': mass_int / MASS_INT_SCALE, 'ppm': ppm, 'proteins': []}
        hits[seq]['proteins'].append((acc, taxid, species_map.get(taxid, f'taxid:{taxid}')))

    print(f"\nQuery: {neutral:.5f} Da ± {args.ppm} ppm  →  {len(hits)} peptide(s)\n")
    for seq, info in sorted(hits.items(), key=lambda x: abs(x[1]['ppm'])):
        print(f"  {seq}  ({info['mass']:.5f} Da, {info['ppm']:+.2f} ppm)")
        for acc, taxid, sp in sorted(set(info['proteins'])):
            print(f"    {acc}  {sp}  [taxid:{taxid}]")
        if args.show_psms:
            psm_rows = con.execute("""
                SELECT run, COUNT(*) n, MIN(frag_score) lo, MAX(frag_score) hi
                FROM psms WHERE seq = ? GROUP BY run ORDER BY run
            """, [seq]).fetchall()
            if psm_rows:
                print(f"    PSMs:")
                for run, n, lo_s, hi_s in psm_rows:
                    print(f"      {run}: {n} PSM(s)  score [{lo_s:.3f}–{hi_s:.3f}]")
            else:
                print(f"    PSMs: none stored")
        print()
    con.close()

# ── Species distribution ──────────────────────────────────────────────────────

def species_dist(args):
    con   = _meta(args.db)
    where = "WHERE p.run = ?" if args.run else ""
    params = [args.run] if args.run else []
    rows = con.execute(f"""
        SELECT COALESCE(NULLIF(a.species,''), 'taxid:' || CAST(p.taxid AS TEXT)) AS sp,
               COUNT(*) AS n
        FROM psms p
        LEFT JOIN accessions a ON a.acc = p.acc
        {where}
        GROUP BY 1 ORDER BY 2 DESC
    """, params).fetchall()
    con.close()

    if not rows:
        print("No PSMs stored yet. Run search_mgf.py first."); return

    title = f"run='{args.run}'" if args.run else "all runs"
    total = sum(r[1] for r in rows)
    print(f"\nSpecies distribution ({title}, {total:,} PSMs)\n")
    for sp, n in rows[:50]:
        bar = '█' * int(30 * n / total)
        print(f"  {n:7,}  {n/total*100:5.1f}%  {bar}  {sp}")
    if len(rows) > 50:
        rest = sum(r[1] for r in rows[50:])
        print(f"  {rest:7,}  {rest/total*100:5.1f}%  ... ({len(rows)-50} more species)")

# ── Export PSMs ───────────────────────────────────────────────────────────────

def export_psms(args):
    con   = _meta(args.db)
    where = "WHERE p.run = ?" if args.run else ""
    params = [args.run] if args.run else []
    rows = con.execute(f"""
        SELECT p.run, p.scan_id, p.title, p.precursor_da, p.charge,
               p.seq, p.acc, p.taxid,
               COALESCE(a.species,'') AS species,
               p.ppm_err, p.ptm, p.ptm_pos, p.n_ions, p.frag_score
        FROM psms p
        LEFT JOIN accessions a ON a.acc = p.acc
        {where}
        ORDER BY p.run, p.scan_id
    """, params).fetchall()
    con.close()
    cols = ['run','scan_id','title','precursor_da','charge','seq','acc',
            'taxid','species','ppm_err','ptm','ptm_pos','n_ions','frag_score']
    with open(args.export_psms, 'w') as fh:
        print('\t'.join(cols), file=fh)
        for row in rows:
            print('\t'.join(str(v) for v in row), file=fh)
    print(f"Exported {len(rows):,} PSMs → {args.export_psms}")

# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(description='Query the peptide database')
    p.add_argument('--db',           required=True)
    p.add_argument('--mass',         type=float, default=None)
    p.add_argument('--mz',           type=float, default=None)
    p.add_argument('--charge',       type=int,   default=None)
    p.add_argument('--ppm',          type=float, default=20.0)
    p.add_argument('--taxid',        type=int,   nargs='+')
    p.add_argument('--show-psms',    action='store_true')
    p.add_argument('--species-dist', action='store_true')
    p.add_argument('--run',          default=None)
    p.add_argument('--export-psms',  default=None, metavar='FILE')
    args = p.parse_args()

    if args.species_dist:    species_dist(args)
    elif args.export_psms:   export_psms(args)
    elif args.mass or (args.mz and args.charge): lookup(args)
    else: p.print_help()

if __name__ == '__main__':
    main()