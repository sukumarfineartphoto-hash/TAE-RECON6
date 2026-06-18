"""SMILES reader - translation of rsmiles.f (readsmile + bondorder subroutines)."""
from ..periodic import TABLE, atomic_number
from ..sparse import SparseMatrix

# Normal valences indexed by atomic number 1-54 (mirroring nvalence in rsmiles.f)
_NVALENCE = [0,
    1, 0, 1, 2, 3, 4, 3, 2, 1, 0,  # H He Li Be B C N O F Ne
    1, 2, 3, 4, 3, 2, 1, 0,         # Na Mg Al Si P S Cl Ar
    1, 2, 3, 4, 5, 4, 4, 3, 2, 2, 2, 2, 3, 4, 3, 2, 1, 0,  # K-Kr
    1, 2, 3, 4, 5, 4, 4, 3, 2, 2, 2, 2, 3, 4, 3, 2, 1, 0,  # Rb-Xe
]


def parse_smiles(smiles_str):
    """Parse a SMILES string, returning a molecule dict.

    Translates readsmile + bondorder + H-saturation block from rsmiles.f.
    Returned dict has same schema as read_sdf (1-based arrays).

    Raises ValueError for syntax this (and the original Fortran rsmiles.f)
    does not support: stereochemistry brackets ([...]), bond-direction
    markers (/ \\), or a string that does not open with an atom.
    """
    line = smiles_str.strip()
    if not line:
        raise ValueError("empty SMILES string")
    if line[0] == '(':
        raise ValueError("SMILES cannot begin with '(' (unsupported by rsmiles.f)")
    if '[' in line or ']' in line:
        raise ValueError("stereochemistry brackets [...] not supported")
    if '/' in line or '\\' in line:
        raise ValueError("bond-direction markers / \\ not supported")

    MAX = 1000
    atom = [None] * (MAX + 1)
    nuc = [0] * (MAX + 1)
    coords = [[0.0, 0.0, 0.0] for _ in range(MAX + 1)]
    ival = [0] * (MAX + 1)
    nbo = SparseMatrix()
    idcon = [[0] * 5 for _ in range(MAX + 1)]
    icon = [[0] * 5 for _ in range(MAX + 1)]
    isum = [0] * (MAX + 1)
    ifree = [0] * (MAX + 1)
    nv = [0] * (MAX + 1)  # heavy-atom bond count per atom
    kbeg = [0] * (MAX + 1)
    kend = [0] * (MAX + 1)

    lstring = len(line)

    # ---- readsmile ----
    n = 1
    k = 0
    while k < lstring:
        s = line[k]
        if s in ('C', 'c'):
            if line[k:k+2] == 'Cl':
                nuc[n] = 17
                kbeg[n] = k + 2
            else:
                nuc[n] = 6
                kbeg[n] = k + 1
        elif s in ('N', 'n'):
            nuc[n] = 7
            kbeg[n] = k + 1
        elif s in ('O', 'o'):
            nuc[n] = 8
            kbeg[n] = k + 1
        elif s == 'F':
            nuc[n] = 9
            kbeg[n] = k + 1
        elif s in ('P', 'p'):
            nuc[n] = 15
            kbeg[n] = k + 1
        elif s == 'I':
            nuc[n] = 53
            kbeg[n] = k + 1
        elif s in ('S', 's'):
            if line[k:k+2] == 'Si':
                nuc[n] = 14
                kbeg[n] = k + 2
            else:
                nuc[n] = 16
                kbeg[n] = k + 1
        elif s in ('B', 'b'):
            if line[k:k+2] == 'Br':
                nuc[n] = 35
                kbeg[n] = k + 2
            else:
                nuc[n] = 5
                kbeg[n] = k + 1
        else:
            k += 1
            continue

        atom[n] = TABLE[nuc[n]]
        # aromatic atoms get valence reduced by 1
        if s in ('c', 'n', 'o', 'p'):
            ival[n] = _NVALENCE[nuc[n]] - 1
        else:
            ival[n] = _NVALENCE[nuc[n]]

        if n != 1:
            kend[n - 1] = k  # Fortran: kend(n-1) = k (1-indexed, current k)
        n += 1
        k += 1
        continue

    kend[n - 1] = k
    natom_heavy = n - 1

    # ---- bondorder ----
    lpar = 0
    npar = [0] * 10  # max ring depth 9
    jring = [0] * 10
    nring = [0] * 10
    mbond = [1] * 10

    for i in range(1, natom_heavy + 1):
        nbond = 1
        jpar = 0
        kpar = 0
        lo = kbeg[i]
        hi = kend[i]
        lo = max(lo, hi)  # mirroring: l2=max(kbeg(i),kend(i))
        # iterate from kbeg[i] to lo (inclusive, 1-based -> 0-based)
        for j in range(kbeg[i] - 1, lo):
            if j < 0 or j >= lstring:
                continue
            s = line[j]
            p = line[j - 1] if j > 0 else ''
            if s == '=':
                nbond = 2
            elif s == '#':
                nbond = 3
            elif s == '(':
                lpar += 1
                if p != ')':
                    npar[lpar] = i
            elif s == ')':
                jpar = 1
                kpar = npar[lpar]
                lpar -= 1
            elif s.isdigit():
                is_ = int(s)
                if 1 <= is_ <= 9:
                    if jring[is_] == 0:
                        mbond[is_] = 2 if p == '=' else (3 if p == '#' else 1)
                        if p in ('=', '#'):
                            nbond = 1
                        nring[is_] = i
                        jring[is_] = 1
                    else:
                        jring[is_] = 0
                        nv[i] += 1
                        nv[nring[is_]] += 1
                        idcon[i][nv[i]] = nring[is_]
                        idcon[nring[is_]][nv[nring[is_]]] = i
                        nbo[i][nring[is_]] = mbond[is_]
                        nbo[nring[is_]][i] = mbond[is_]

        if i == natom_heavy:
            continue
        if jpar == 0:
            nv[i] += 1
            nv[i + 1] += 1
            idcon[i][nv[i]] = i + 1
            idcon[i + 1][nv[i + 1]] = i
            nbo[i][i + 1] = nbond
            nbo[i + 1][i] = nbond
        else:
            nv[i + 1] += 1
            nv[kpar] += 1
            idcon[i + 1][nv[i + 1]] = kpar
            idcon[kpar][nv[kpar]] = i + 1
            nbo[i + 1][kpar] = nbond
            nbo[kpar][i + 1] = nbond

    # Symmetrize nbo
    for i in range(1, natom_heavy + 1):
        for j in range(i + 1, natom_heavy + 1):
            nbo[j][i] = nbo[i][j]

    # Fix valences and count free valencies
    jfree = 0
    for i in range(1, natom_heavy + 1):
        for j in range(1, natom_heavy + 1):
            if nbo[i][j] > 1:
                ival[i] -= (nbo[i][j] - 1)
        if ival[i] < nv[i]:
            ival[i] = nv[i]
        ifree[i] = ival[i] - nv[i]
        if ifree[i] < 0:
            ifree[i] = 0
        if ival[i] > 4:
            ival[i] = 4
        jfree += ifree[i]

    # Saturate with Hs
    hi = natom_heavy
    for j in range(1, natom_heavy + 1):
        while ifree[j] > 0:
            hi += 1
            nuc[hi] = 1
            ival[hi] = 1
            atom[hi] = 'H '
            idcon[hi][1] = j
            nv[j] += 1
            idcon[j][nv[j]] = hi
            ifree[j] -= 1
            nbo[hi][j] = 1
            nbo[j][hi] = 1

    natom = hi
    numbond = 0
    for l in range(1, natom + 1):
        for j in range(l + 1, natom + 1):
            if nbo[l][j] != 0:
                numbond += 1

    # Rebuild ival from nv for heavy atoms + H
    for i in range(1, natom + 1):
        ival[i] = nv[i] if i <= natom_heavy else 1

    # Build icon / isum
    for i in range(1, natom + 1):
        s = 0
        for j in range(1, ival[i] + 1):
            nb = nuc[idcon[i][j]]
            icon[i][j] = nb
            s += nb
        isum[i] = s

    return dict(natom=natom, numbond=numbond, atom=atom, nuc=nuc,
                coords=coords, idcon=idcon, ival=ival, nbo=nbo,
                icon=icon, isum=isum)
