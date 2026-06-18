"""MOL2 (Tripos SYBYL) reader - translation of rmol2.f."""
from ..element import ele
from ..periodic import atomic_number
from ..sparse import SparseMatrix


def read_mol2(fileobj):
    """Read next molecule from an open .mol2 fileobj.
    Returns same dict shape as read_sdf, or raises StopIteration at EOF."""
    MAX = 1000
    atom = [None] * (MAX + 1)
    nuc = [0] * (MAX + 1)
    coords = [[0.0, 0.0, 0.0] for _ in range(MAX + 1)]
    idcon = [[0] * 5 for _ in range(MAX + 1)]
    ival = [0] * (MAX + 1)
    nbo = SparseMatrix()
    icon = [[0] * 5 for _ in range(MAX + 1)]
    isum = [0] * (MAX + 1)

    # Advance to @<TRIPOS>MOLECULE
    found = False
    while True:
        line = fileobj.readline()
        if not line:
            raise StopIteration
        if line.strip() == '@<TRIPOS>MOLECULE':
            found = True
            break
    if not found:
        raise StopIteration

    fileobj.readline()  # mol name
    counts_line = fileobj.readline().split()
    natom = int(counts_line[0])
    numbond = int(counts_line[1]) if len(counts_line) > 1 else 0

    # Find @<TRIPOS>ATOM
    while True:
        line = fileobj.readline()
        if not line:
            raise StopIteration
        if line.strip() == '@<TRIPOS>ATOM':
            break

    for i in range(1, natom + 1):
        line = fileobj.readline()
        parts = line.split()
        # fields: id name x y z sybyl_type [subst_id subst_name charge]
        coords[i][0] = float(parts[2])
        coords[i][1] = float(parts[3])
        coords[i][2] = float(parts[4])
        sybyl = parts[5] if len(parts) > 5 else parts[1]
        # element is first part of sybyl type before '.'
        elem_raw = sybyl.split('.')[0]
        atom[i] = ele(elem_raw)

    # Assign atomic numbers
    for i in range(1, natom + 1):
        an = atomic_number(atom[i])
        if an is not None:
            nuc[i] = an

    # Find @<TRIPOS>BOND
    while True:
        line = fileobj.readline()
        if not line:
            break
        if line.strip() == '@<TRIPOS>BOND':
            break

    for _ in range(numbond):
        line = fileobj.readline()
        if not line:
            break
        parts = line.split()
        iat, jat = int(parts[1]), int(parts[2])
        bo_str = parts[3] if len(parts) > 3 else '1'
        try:
            bo = int(bo_str)
        except ValueError:
            bo = 1  # 'ar', 'am', etc.
        nbo[iat][jat] = bo
        nbo[jat][iat] = bo
        if ival[iat] < 4 and ival[jat] < 4:
            ival[iat] += 1
            ival[jat] += 1
            idcon[iat][ival[iat]] = jat
            idcon[jat][ival[jat]] = iat

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
