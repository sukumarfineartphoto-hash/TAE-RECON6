"""
ORCA input file reader - Cartesian ("* xyz") coordinate block format.

An ORCA input file's molecule specification has this shape:

    ! method basis keywords    <- simple keyword line(s), any number
    %scf                       <- optional %block ... end sections
      MaxIter 200
    end
    * xyz 0 1                  <- "* xyz <charge> <multiplicity>" marker
    C   0.000  0.000  0.000    <- atom block: element symbol OR atomic
    H   0.629  0.629  0.629       number, followed by x, y, z
    ...
    *                          <- bare "*" terminates the atom block

Unlike Gaussian's fixed blank-line-delimited sections, ORCA's
Cartesian block is identified by scanning for a line starting with
"*" whose first token (after stripping leading "*"s and whitespace)
is "xyz" (case-insensitive; both "* xyz" and "*xyz" are accepted).
Everything before that marker (route lines, %block...end sections,
comments) is skipped without needing to be parsed. The charge and
multiplicity are the two tokens following "xyz" on that same line.
The atom block then reads element-symbol-or-atomic-number + x + y + z
lines until a line that is just "*" (optionally with surrounding
whitespace).

As with the Gaussian reader, only Cartesian input is supported - ORCA
also supports internal coordinates and a few other input styles, which
are out of scope here. Connectivity has no bond list in this format
either, so it is always determined the same way as a CONECT-less PDB
file: distance-based bonding plus bond-order inference (connectivity.py).
"""
from ..sparse import SparseMatrix
from ..connectivity import determine_connectivity
from .gaussian import _parse_atom_symbol_or_number
from ..periodic import atomic_number


def read_orca_xyz(filepath, bondl=None, bond_length=1.09):
    """
    Read an ORCA input file's Cartesian ("* xyz") coordinate block.
    Returns the same dict shape as read_sdf.

    Connectivity is always distance-based (ORCA xyz input has no bond
    list), using `bondl` exactly as for a CONECT-less PDB file or
    Gaussian Cartesian input.

    Note: unlike SDF/MOL2/PDB, this reader never adds hydrogens. ORCA
    input conventionally specifies every atom including hydrogens
    explicitly, so automatic H-addition is not offered here at all (a
    heavy-atom-only ORCA file is not a supported input - it would
    silently misrepresent the actual quantum chemistry input if
    "completed" by guesswork).
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

    # Scan forward for the "* xyz <charge> <mult>" marker line,
    # skipping any preceding keyword lines or %block...end sections.
    mol_charge = None
    mol_mult = None
    while True:
        line = next_line()
        if line is None:
            raise ValueError("No '* xyz <charge> <multiplicity>' block found")
        stripped = line.strip()
        if not stripped.startswith('*'):
            continue
        # Strip leading '*' characters and whitespace, then check for
        # an "xyz" token (handles both "* xyz" and "*xyz").
        after_star = stripped.lstrip('*').strip()
        parts = after_star.split()
        if len(parts) >= 3 and parts[0].lower() == 'xyz':
            try:
                mol_charge = int(parts[1])
                mol_mult = int(parts[2])
            except ValueError:
                raise ValueError("Malformed '* xyz' line: %r" % line)
            break
        # A bare '*' or a '*'-prefixed line that isn't an xyz marker
        # (e.g. closing a different block) is just skipped.

    # Atom block: element symbol or atomic number, then x, y, z;
    # terminated by a bare '*' line (or EOF).
    natom = 0
    while True:
        line = next_line()
        if line is None:
            break
        stripped = line.strip()
        if stripped == '*':
            break
        if stripped == '':
            continue
        parts = stripped.split()
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
        raise ValueError("No atoms found in ORCA '* xyz' coordinate block")

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
