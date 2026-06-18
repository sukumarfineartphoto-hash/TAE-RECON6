"""
Gaussian .com/.gjf input file reader - Cartesian coordinate format only.

A Gaussian input file's molecule specification section has this shape:

    %chk=...                  <- link0 commands (any number, optional)
    %mem=...
    # route line / method     <- route section (one or more lines)
                               <- blank line (ends route section)
    Title line                <- free-text title (any number of lines)
                               <- blank line (ends title section)
    0 1                       <- charge and multiplicity
    C   0.000  0.000  0.000   <- atom block: element symbol OR atomic
    H   0.629  0.629  0.629      number, followed by x, y, z
    ...
                               <- blank line (ends atom block)

This reader skips everything up through the blank-line-then-charge/
multiplicity pattern, then reads atoms (by symbol or integer atomic
number) until the terminating blank line. Z-matrix (internal
coordinate) input is intentionally not supported - only Cartesian
coordinates, which is what RECON's actual use case needs.

Since a Gaussian input has no bond list at all (just atoms and
coordinates), connectivity is always determined the same way as a
CONECT-less PDB file: distance-based bonding plus bond-order inference
(see connectivity.py) using the bond-length reference table.
"""
from ..periodic import TABLE, TABLECAP, atomic_number
from ..sparse import SparseMatrix
from ..connectivity import determine_connectivity


def _parse_atom_symbol_or_number(token):
    """Return a 2-char element symbol for a Gaussian atom-spec token,
    which may be an element symbol (any case, e.g. 'C', 'Cl', 'CL') or
    an integer atomic number (e.g. '6' for carbon). Returns None if it
    can't be resolved to a known element (1-54).

    The returned symbol is always normalized to the canonical
    upper-case form (e.g. 'CL', 'BR'), matching the convention used
    throughout the rest of the codebase (see element.py).
    """
    try:
        n = int(token)
        if 1 <= n <= 54:
            return TABLECAP[n]
        return None
    except ValueError:
        pass
    sym = (token + '  ')[:2]
    if atomic_number(sym) is not None:
        return sym.upper()
    if atomic_number(sym.upper()) is not None:
        return sym.upper()
    # Try matching case-insensitively against the canonical table
    for i in range(1, 55):
        if TABLE[i].strip().upper() == token.strip().upper():
            return TABLECAP[i]
    return None


def read_gaussian_com(filepath, bondl=None, bond_length=1.09):
    """
    Read a Gaussian .com/.gjf file (Cartesian coordinates only).
    Returns the same dict shape as read_sdf.

    Connectivity is always distance-based (Gaussian input has no bond
    list), using `bondl` (the standard RECON bond-length reference
    table) exactly as for a CONECT-less PDB file.

    Note: unlike SDF/MOL2/PDB, this reader never adds hydrogens.
    Gaussian input conventionally specifies every atom including
    hydrogens explicitly, so automatic H-addition is not offered here
    at all (a heavy-atom-only Gaussian file is not a supported input -
    it would silently misrepresent the actual quantum chemistry input
    if "completed" by guesswork).
    """
    MAX = 1000
    atom = [None] * (MAX + 1)
    nuc = [0] * (MAX + 1)
    coords = [[0.0, 0.0, 0.0] for _ in range(MAX + 1)]
    idcon = [[0] * 5 for _ in range(MAX + 1)]
    ival = [0] * (MAX + 1)
    nbo = SparseMatrix()
    icon = [[0] * 5 for _ in range(MAX + 1)]
    isum = [0] * (MAX + 1)

    with open(filepath) as fh:
        lines = fh.readlines()

    pos = 0
    n = len(lines)

    def next_line():
        nonlocal pos
        if pos >= n:
            return None
        line = lines[pos]
        pos += 1
        return line

    def is_blank(line):
        return line is not None and line.strip() == ''

    # Skip link0/route lines through the blank line that ends them
    line = next_line()
    while line is not None and not is_blank(line):
        line = next_line()
    if line is None:
        raise ValueError("Gaussian input ended before route section blank line")

    # Skip the title section through its terminating blank line
    line = next_line()
    while line is not None and not is_blank(line):
        line = next_line()
    if line is None:
        raise ValueError("Gaussian input ended before title section blank line")

    # Charge/multiplicity line: two integers
    charge_mult_line = next_line()
    if charge_mult_line is None:
        raise ValueError("Missing charge/multiplicity line")
    parts = charge_mult_line.split()
    if len(parts) < 2:
        raise ValueError("Malformed charge/multiplicity line: %r" % charge_mult_line)
    try:
        mol_charge = int(parts[0])
        mol_mult = int(parts[1])
    except ValueError:
        raise ValueError("Malformed charge/multiplicity line: %r" % charge_mult_line)

    # Atom block: element symbol or atomic number, then x, y, z;
    # terminated by a blank line (or EOF).
    natom = 0
    while True:
        line = next_line()
        if line is None or is_blank(line):
            break
        parts = line.split()
        if len(parts) < 4:
            raise ValueError("Malformed atom line: %r" % line)
        sym = _parse_atom_symbol_or_number(parts[0])
        if sym is None:
            raise ValueError("Unrecognized element/atomic number: %r" % parts[0])
        natom += 1
        if natom > MAX:
            raise ValueError("Too many atoms (> %d)" % MAX)
        atom[natom] = sym
        nuc[natom] = atomic_number(sym)
        coords[natom][0] = float(parts[1])
        coords[natom][1] = float(parts[2])
        coords[natom][2] = float(parts[3])

    if natom == 0:
        raise ValueError("No atoms found in Gaussian input atom block")

    numbond = 0
    if bondl is not None:
        numbond = determine_connectivity(natom, nuc, coords, bondl, idcon, ival, nbo)

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
               icon=icon, isum=isum, charge={}, mol_charge=mol_charge,
               mol_multiplicity=mol_mult)

    return mol
