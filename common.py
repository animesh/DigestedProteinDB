"""common.py – shared mass tables, digestion, and fragment scoring."""

PROTON = 1.007276467
H2O    = 18.010564684

AA_MASS = {
    'A':  71.03711379,  'C': 103.00918448,  'D': 115.02694303,
    'E': 129.04259310,  'F': 147.06841392,  'G':  57.02146372,
    'H': 137.05891186,  'I': 113.08406398,  'K': 128.09496302,
    'L': 113.08406398,  'M': 131.04048461,  'N': 114.04292745,
    'P':  97.05276385,  'Q': 128.05857751,  'R': 156.10111103,
    'S':  87.03202841,  'T': 101.04767847,  'V':  99.06841392,
    'W': 186.07931295,  'Y': 163.06332854,
}

# name -> (delta_mass_da, residues_where_applicable)
# None for residues means N-terminal modification
PTMS = {
    'ox_M':     (15.994915, 'M'),
    'phospho':  (79.966331, 'STY'),
    'deamid':   ( 0.984016, 'NQ'),
    'acetyl_K': (42.010565, 'K'),
    'carbam_C': (57.021464, 'C'),
}

# ── Digestion ────────────────────────────────────────────────────────────────

def _cleavage_sites(seq: str, enzyme: str) -> list:
    sites = [0]
    n = len(seq)
    if enzyme == 'trypsin':
        for i in range(n - 1):
            if seq[i] in 'RK' and seq[i + 1] != 'P':
                sites.append(i + 1)
    elif enzyme == 'chymotrypsin':
        for i in range(n - 1):
            if seq[i] in 'FWY' and seq[i + 1] != 'P':
                sites.append(i + 1)
    else:
        raise ValueError(f'Unknown enzyme: {enzyme}')
    sites.append(n)
    return sites


def digest(seq: str, min_len: int, max_len: int, mc: int,
           enzyme: str = 'trypsin') -> list:
    """
    Return tryptic (or chymotryptic) peptides within [min_len, max_len].
    A peptide is kept only if every residue has a known mass in AA_MASS.
    This makes the filter self-consistent with neutral_mass() and handles
    any non-standard residue (U, O, B, Z, J, X, ...) without a hardcoded
    exclusion list that can silently miss newly encountered characters.
    """
    sites = _cleavage_sites(seq, enzyme)
    seen, result = set(), []
    for n_mc in range(mc + 1):
        for i in range(len(sites) - 1 - n_mc):
            pep = seq[sites[i]: sites[i + 1 + n_mc]]
            if (min_len <= len(pep) <= max_len
                    and all(aa in AA_MASS for aa in pep)
                    and pep not in seen):
                seen.add(pep)
                result.append(pep)
    return result


# ── Mass ─────────────────────────────────────────────────────────────────────

def neutral_mass(seq: str) -> float:
    return sum(AA_MASS[aa] for aa in seq) + H2O


def precursor_neutral(pepmass: float, charge: int) -> float:
    """Neutral monoisotopic mass from MGF PEPMASS and CHARGE."""
    return pepmass * charge - charge * PROTON


def ppm_bounds(mass: float, ppm: float) -> tuple:
    delta = mass * ppm / 1e6
    return round((mass - delta) * 1e4), round((mass + delta) * 1e4)


def ppm_error(observed: float, theoretical: float) -> float:
    return (observed - theoretical) / theoretical * 1e6


# ── Fragment ions ─────────────────────────────────────────────────────────────

def b_ions(seq: str, ptm_mass: float = 0.0, ptm_pos: int = -1) -> list:
    """
    Singly-charged b-ions for seq.
    ptm_pos: 0-indexed position of modification on the sequence.
    ptm_mass added to residue at ptm_pos; propagates to all b-ions that
    span that position.
    """
    ions, mass = [], 0.0
    for i, aa in enumerate(seq[:-1]):
        mass += AA_MASS[aa]
        if i == ptm_pos:
            mass += ptm_mass
        ions.append(mass + PROTON)
    return ions


def y_ions(seq: str, ptm_mass: float = 0.0, ptm_pos: int = -1) -> list:
    """
    Singly-charged y-ions for seq.
    ptm_pos: 0-indexed from N-terminus.
    """
    ions, mass = [], H2O
    n = len(seq)
    for i in range(n - 1, 0, -1):          # C-terminal residues first
        mass += AA_MASS[seq[i]]
        if i == ptm_pos:
            mass += ptm_mass
        ions.append(mass + PROTON)
    return ions


def _count_matches(theory: list, obs: list, tol_ppm: float) -> tuple:
    """(n_matched, frac_intensity_matched) against observed [(mz, intensity)]."""
    if not obs or not theory:
        return 0, 0.0
    total_int = sum(i for _, i in obs)
    matched_int, n = 0.0, 0
    for omz, oint in obs:
        for tmz in theory:
            if abs(omz - tmz) / tmz * 1e6 <= tol_ppm:
                n += 1
                matched_int += oint
                break
    return n, matched_int / total_int if total_int else 0.0


def score_fragments(seq: str, obs: list, tol_ppm: float = 20.0,
                    ptm_mass: float = 0.0, ptm_sites: str = '') -> tuple:
    """
    Score sequence against observed MS2 fragments.

    For unmodified candidates (ptm_mass == 0): straightforward b/y matching.

    For PTM candidates: try placing the modification on every applicable
    residue and return the score for the best-localised position.

    Returns (n_matched_ions, frac_intensity_matched, best_ptm_pos).
    best_ptm_pos is -1 when no PTM.
    """
    if not obs:
        return 0, 0.0, -1

    if ptm_mass == 0.0 or not ptm_sites:
        b = b_ions(seq)
        y = y_ions(seq)
        n, frac = _count_matches(b + y, obs, tol_ppm)
        return n, frac, -1

    # Try each applicable site, keep the best
    best_n, best_frac, best_pos = 0, 0.0, -1
    for pos, aa in enumerate(seq):
        if aa not in ptm_sites:
            continue
        b = b_ions(seq, ptm_mass, pos)
        y = y_ions(seq, ptm_mass, pos)
        n, frac = _count_matches(b + y, obs, tol_ppm)
        if n > best_n or (n == best_n and frac > best_frac):
            best_n, best_frac, best_pos = n, frac, pos

    return best_n, best_frac, best_pos