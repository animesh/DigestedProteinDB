"""common.py – mass tables, digestion, fragment scoring, pack/unpack."""

PROTON = 1.007276467
H2O    = 18.010564684
MASS_INT_SCALE = 1e4

AA_MASS = {
    'A':  71.03711379,  'C': 103.00918448,  'D': 115.02694303,
    'E': 129.04259310,  'F': 147.06841392,  'G':  57.02146372,
    'H': 137.05891186,  'I': 113.08406398,  'K': 128.09496302,
    'L': 113.08406398,  'M': 131.04048461,  'N': 114.04292745,
    'P':  97.05276385,  'Q': 128.05857751,  'R': 156.10111103,
    'S':  87.03202841,  'T': 101.04767847,  'V':  99.06841392,
    'W': 186.07931295,  'Y': 163.06332854,
}

PTMS = {
    'ox_M':     (15.994915, 'M'),
    'phospho':  (79.966331, 'STY'),
    'deamid':   ( 0.984016, 'NQ'),
    'acetyl_K': (42.010565, 'K'),
    'carbam_C': (57.021464, 'C'),
}

# ── Digestion ─────────────────────────────────────────────────────────────────

def digest(seq: str, min_len: int, max_len: int, mc: int,
           enzyme: str = 'trypsin') -> list:
    """
    Return peptides whose every residue has a known mass in AA_MASS.
    This makes the filter self-consistent with neutral_mass() and
    handles any non-standard residue automatically.
    """
    if enzyme == 'trypsin':
        sites = [0]
        for i in range(len(seq) - 1):
            if seq[i] in 'RK' and seq[i+1] != 'P':
                sites.append(i+1)
    elif enzyme == 'chymotrypsin':
        sites = [0]
        for i in range(len(seq) - 1):
            if seq[i] in 'FWY' and seq[i+1] != 'P':
                sites.append(i+1)
    else:
        raise ValueError(f'Unknown enzyme: {enzyme}')
    sites.append(len(seq))

    seen, result = set(), []
    for n_mc in range(mc + 1):
        for i in range(len(sites) - 1 - n_mc):
            pep = seq[sites[i]: sites[i+1+n_mc]]
            if (min_len <= len(pep) <= max_len
                    and all(aa in AA_MASS for aa in pep)
                    and pep not in seen):
                seen.add(pep)
                result.append(pep)
    return result

# ── Mass ──────────────────────────────────────────────────────────────────────

def neutral_mass(seq: str) -> float:
    return sum(AA_MASS[aa] for aa in seq) + H2O

def precursor_neutral(pepmass: float, charge: int) -> float:
    return pepmass * charge - charge * PROTON

def ppm_bounds(mass: float, ppm: float) -> tuple:
    delta = mass * ppm / 1e6
    return round((mass - delta) * MASS_INT_SCALE), round((mass + delta) * MASS_INT_SCALE)

def ppm_error(observed: float, theoretical: float) -> float:
    return (observed - theoretical) / theoretical * 1e6

# ── RocksDB key helpers ───────────────────────────────────────────────────────

def mass_key(mass_int: int) -> bytes:
    """4-byte big-endian key: integer order == mass order."""
    return mass_int.to_bytes(4, 'big')

def key_to_mass(key: bytes) -> int:
    return int.from_bytes(key, 'big')

# ── Binary value packing ──────────────────────────────────────────────────────

def pack_entry(seq: str, acc: str, taxid: int) -> bytes:
    """
    Pack one (seq, acc, taxid) into compact binary.
    Format: [seq_len:1][seq][acc_len:1][acc][taxid:4be]
    Max: 1+30+1+15+4 = 51 bytes per entry.
    """
    sb = seq.encode('ascii')
    ab = acc.encode('ascii')
    return bytes([len(sb)]) + sb + bytes([len(ab)]) + ab + taxid.to_bytes(4, 'big')

def unpack_entries(blob: bytes):
    """Yield (seq, acc, taxid) from a packed binary blob."""
    i = 0
    while i < len(blob):
        sl = blob[i]; i += 1
        seq = blob[i:i+sl].decode('ascii'); i += sl
        al = blob[i]; i += 1
        acc = blob[i:i+al].decode('ascii'); i += al
        taxid = int.from_bytes(blob[i:i+4], 'big'); i += 4
        yield seq, acc, taxid

# ── Fragment ions ─────────────────────────────────────────────────────────────

def b_ions(seq: str, ptm_mass: float = 0.0, ptm_pos: int = -1) -> list:
    ions, mass = [], 0.0
    for i, aa in enumerate(seq[:-1]):
        mass += AA_MASS[aa]
        if i == ptm_pos:
            mass += ptm_mass
        ions.append(mass + PROTON)
    return ions

def y_ions(seq: str, ptm_mass: float = 0.0, ptm_pos: int = -1) -> list:
    ions, mass = [], H2O
    for i in range(len(seq)-1, 0, -1):
        mass += AA_MASS[seq[i]]
        if i == ptm_pos:
            mass += ptm_mass
        ions.append(mass + PROTON)
    return ions

def score_fragments(seq: str, obs: list, tol_ppm: float = 20.0,
                    ptm_mass: float = 0.0, ptm_sites: str = '') -> tuple:
    """Returns (n_matched, frac_intensity, best_ptm_pos)."""
    if not obs:
        return 0, 0.0, -1
    if ptm_mass == 0.0 or not ptm_sites:
        theory = b_ions(seq) + y_ions(seq)
        total = sum(i for _, i in obs)
        matched_int, n = 0.0, 0
        for omz, oint in obs:
            for tmz in theory:
                if abs(omz - tmz) / tmz * 1e6 <= tol_ppm:
                    n += 1; matched_int += oint; break
        return n, matched_int/total if total else 0.0, -1
    best_n, best_frac, best_pos = 0, 0.0, -1
    total = sum(i for _, i in obs)
    for pos, aa in enumerate(seq):
        if aa not in ptm_sites:
            continue
        theory = b_ions(seq, ptm_mass, pos) + y_ions(seq, ptm_mass, pos)
        matched_int, n = 0.0, 0
        for omz, oint in obs:
            for tmz in theory:
                if abs(omz - tmz) / tmz * 1e6 <= tol_ppm:
                    n += 1; matched_int += oint; break
        if n > best_n or (n == best_n and matched_int > best_frac * total):
            best_n, best_frac, best_pos = n, matched_int/total if total else 0.0, pos
    return best_n, best_frac, best_pos