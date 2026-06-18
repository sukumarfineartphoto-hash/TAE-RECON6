"""PDB reader - translation of rpdb.f + pdbread.f."""
from ..element import ele
from ..periodic import atomic_number
from ..sparse import SparseMatrix


def _element_from_pdb_line(line, atm_name):
    """Determine the element for a PDB ATOM/HETATM line.

    The original Fortran (rpdb.f) only ever looked at the first 2
    characters of the atom-name field (columns 13-14) and ran them
    through ele.f's name-based heuristic. That heuristic silently fails
    for atoms whose *name* doesn't start with their element symbol -
    notably deuterium, which PDB commonly names like "DNAD" (deuterium
    bonded to atom NAD) rather than starting with "D" as a standalone
    symbol in the way ele.f expects for other elements.

    This port additionally consults the dedicated element column
    (PDB columns 77-78, 0-indexed [76:78]) when present, which is the
    authoritative source for this. Deuterium ('D') is folded into
    hydrogen, since it's chemically just an isotope of H and the TAE
    database has no separate deuterium atom-type entries.
    """
    elem_col = line[76:78].strip() if len(line) >= 78 else ''
    if elem_col:
        if elem_col.upper() == 'D':
            return 'H '
        sym = (elem_col + '  ')[:2]
        # Element column is authoritative when it parses to a known
        # element; otherwise fall through to the name-based heuristic.
        if atomic_number(sym) is not None or atomic_number(sym.upper()) is not None:
            return sym
    return ele(atm_name)


def read_pdb(filepath, bondl=None, iovr=0, auto_add_h=False, bond_length=1.09):
    """Read a PDB file. Returns same dict shape as read_sdf.

    If bondl is provided (54x54 bond-length table) and iovr > 0,
    connectivity is determined by distance criterion (pdbread.f logic).
    Otherwise CONECT records are used.

    If auto_add_h is True, missing hydrogens are added using the same
    standard-valence logic as SDF/MOL2 input (see hydrogenate.py) when
    the structure has more than 10 atoms and zero explicit hydrogens
    (deuterium counts as hydrogen for this check, so a partially
    deuterium-labeled structure is left alone rather than risking a
    double-saturation - matching the same conservative behavior used
    for SDF/MOL2 input).
    """
    MAX = 1000
    atm_raw = [None] * (MAX + 1)
    atom = [None] * (MAX + 1)
    nuc = [0] * (MAX + 1)
    coords = [[0.0, 0.0, 0.0] for _ in range(MAX + 1)]
    idcon = [[0] * 5 for _ in range(MAX + 1)]
    ival = [0] * (MAX + 1)
    nbo = SparseMatrix()
    icon = [[0] * 5 for _ in range(MAX + 1)]
    isum = [0] * (MAX + 1)
    natom = 0

    with open(filepath) as fh:
        lines = fh.readlines()

    conect_lines = []
    raw_lines_by_atom = {}
    for line in lines:
        rec = line[:6]
        if rec in ('ATOM  ', 'HETATM'):
            if line[12:14].strip() == 'LP':
                continue
            natom += 1
            atm_raw[natom] = line[12:14].strip()
            raw_lines_by_atom[natom] = line
            coords[natom][0] = float(line[30:38])
            coords[natom][1] = float(line[38:46])
            coords[natom][2] = float(line[46:54])
        elif rec in ('CONECT', 'BOND  '):
            conect_lines.append(line)

    # Normalise element symbols (element column takes priority over
    # the name-based heuristic when present and recognized)
    for i in range(1, natom + 1):
        atom[i] = _element_from_pdb_line(raw_lines_by_atom[i], atm_raw[i])

    # Assign atomic numbers
    for i in range(1, natom + 1):
        an = atomic_number(atom[i])
        if an is not None:
            nuc[i] = an

    numbond = 0
    use_distance = (bondl is not None and iovr > 0)
    if not use_distance and not conect_lines and bondl is not None:
        # No CONECT records in the file at all, and the caller didn't
        # explicitly request CONECT-only behavior (iovr <= 0 just means
        # "no override was given" in that case) - rather than silently
        # returning a fully disconnected molecule, fall back to
        # distance-based connectivity automatically. This matches the
        # common case (PDB ligand extracts with no CONECT records).
        use_distance = True

    if use_distance:
        from ..connectivity import determine_connectivity
        numbond = determine_connectivity(natom, nuc, coords, bondl, idcon, ival, nbo)
    else:
        # CONECT records
        idbond = [[False] * (natom + 1) for _ in range(natom + 1)]
        for line in conect_lines:
            try:
                n = int(line[6:11])
            except ValueError:
                continue
            for col in (11, 16, 21, 26):
                s = line[col:col+5].strip()
                if not s:
                    break
                try:
                    j = int(s)
                except ValueError:
                    break
                if j == 0:
                    break
                ival[n] = min(ival[n] + 1, 4)
                idcon[n][ival[n]] = j
                if not idbond[n][j] and not idbond[j][n]:
                    numbond += 1
                    idbond[n][j] = idbond[j][n] = True

    # Build icon / isum
    for i in range(1, natom + 1):
        s = 0
        for j in range(1, ival[i] + 1):
            nb = nuc[idcon[i][j]]
            icon[i][j] = nb
            s += nb
        isum[i] = s

    mol = dict(natom=natom, numbond=numbond, atom=atom, nuc=nuc,
               coords=coords, idcon=idcon, ival=ival, nbo=nbo,
               icon=icon, isum=isum, charge={})

    if auto_add_h:
        from ..hydrogenate import needs_hydrogens, add_missing_hydrogens
        if needs_hydrogens(mol):
            n_before = mol['natom']
            mol = add_missing_hydrogens(mol, bond_length=bond_length)
            mol['_pdb_hydrogens_added'] = mol['hydrogens_added']
            mol['_pdb_natom_before_h'] = n_before

    return mol
