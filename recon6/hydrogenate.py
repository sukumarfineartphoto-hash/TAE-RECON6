"""
Hydrogen-addition utility for molecules read without explicit
hydrogens (e.g. many SDF files only list heavy atoms).

Uses standard valence rules (the same default valences RDKit/OpenBabel
apply for implicit-H calculation) rather than replicating the original
Fortran's behavior, which simply logged a warning and left such
molecules unsaturated. New hydrogen positions are placed using local
geometry heuristics (linear / trigonal planar / tetrahedral) based on
each atom's existing bonded-neighbor count and geometry.

This is a standalone post-processing step: call `add_missing_hydrogens`
on a molecule dict (as returned by readers.sdf.read_sdf or
readers.mol2.read_mol2) before it's fed into the rest of the RECON5
pipeline. PDB/SMILES inputs are not run through this — PDB structures
normally already specify hydrogens or rely on distance-based
connectivity, and SMILES already gets H-saturated by readers/smiles.py.
"""
import math

from .sparse import SparseMatrix

# Standard (default) valence per element, used to compute how many
# implicit hydrogens an atom needs. Matches common RDKit/OpenBabel
# default valences for neutral atoms in their most common bonding state.
_STANDARD_VALENCE = {
    1: 1,    # H
    5: 3,    # B
    6: 4,    # C
    7: 3,    # N
    8: 2,    # O
    9: 1,    # F
    14: 4,   # Si
    15: 3,   # P  (5 also common, e.g. phosphate - handled via charge/explicit bonds)
    16: 2,   # S  (often 2; higher valences in sulfoxides/sulfones are
             #     usually explicit in the input bond table already)
    17: 1,   # Cl
    35: 1,   # Br
    53: 1,   # I
}

# Elements that commonly appear with a higher standard valence when
# they already have 3+ heavy-atom connections (e.g. P in phosphates,
# S in sulfones) - used as a secondary lookup so we don't try to stuff
# extra H's onto an already fully-substituted center.
_EXTENDED_VALENCE = {
    15: (3, 5),
    16: (2, 4, 6),
}


def _bond_order_sum(mol, i):
    """Total bond order (sum of all bond multiplicities) for atom i."""
    total = 0
    for j in range(1, mol['ival'][i] + 1):
        nbr = mol['idcon'][i][j]
        bo = mol['nbo'][i][nbr]
        total += bo if bo else 1
    return total


def _best_standard_valence(elem, used):
    """Pick the smallest standard/extended valence for elem that is
    >= used (so we don't ask for a negative number of hydrogens)."""
    if elem in _EXTENDED_VALENCE:
        options = _EXTENDED_VALENCE[elem]
    elif elem in _STANDARD_VALENCE:
        options = (_STANDARD_VALENCE[elem],)
    else:
        return None
    for v in sorted(options):
        if v >= used:
            return v
    return max(options)


def needs_hydrogens(mol, threshold_atoms=10):
    """
    Decide whether a molecule looks like it's missing hydrogens.

    Mirrors the trigger condition the original Fortran used for its
    warning (more than `threshold_atoms` atoms and zero hydrogens
    present), but this is just the *detection* heuristic - remediation
    differs (see module docstring).
    """
    natom = mol['natom']
    num_h = sum(1 for i in range(1, natom + 1) if mol['nuc'][i] == 1)
    if natom > threshold_atoms and num_h == 0:
        return True
    return False


def _unit(v):
    n = math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
    if n < 1e-9:
        return [0.0, 0.0, 1.0]
    return [v[0]/n, v[1]/n, v[2]/n]


def _cross(a, b):
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]


def _add(a, b, scale=1.0):
    return [a[0]+b[0]*scale, a[1]+b[1]*scale, a[2]+b[2]*scale]


def _arbitrary_perp(v):
    """Return a unit vector perpendicular to v."""
    ref = [1.0, 0.0, 0.0] if abs(v[0]) < 0.9 else [0.0, 1.0, 0.0]
    p = _cross(v, ref)
    return _unit(p)


def _new_h_directions(center, neighbor_positions, n_new, bond_length=1.09):
    """
    Compute n_new unit-vector directions for placing new hydrogens on
    `center`, given the existing bonded-neighbor positions (list of
    [x,y,z]). Uses simple geometric heuristics based on total
    coordination number (existing neighbors + new H's):
      - 1 total: arbitrary direction
      - 2 total (linear, e.g. on a 1-coordinate center): opposite the
        single existing neighbor
      - 3 total (trigonal planar): in-plane, 120 degrees apart
      - 4 total (tetrahedral): standard tetrahedral geometry
    """
    n_exist = len(neighbor_positions)
    n_total = n_exist + n_new
    existing_dirs = [_unit([p[0]-center[0], p[1]-center[1], p[2]-center[2]])
                     for p in neighbor_positions]

    if n_total <= 1:
        return [[0.0, 0.0, 1.0]][:n_new]

    if n_exist == 0:
        # No existing neighbors at all (isolated atom needing all its
        # H's) - distribute evenly around a tetrahedron / arbitrary axis.
        base_dirs = {
            1: [[0, 0, 1]],
            2: [[1, 0, 0], [-1, 0, 0]],
            3: [[1, 0, 0], [-0.5, 0.866, 0], [-0.5, -0.866, 0]],
            4: [[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]],
        }
        dirs = base_dirs.get(n_total, [[1, 0, 0]] * n_total)
        return [_unit(d) for d in dirs][:n_new]

    if n_exist == 1:
        d0 = existing_dirs[0]
        if n_total == 2:
            # Linear: new H opposite the existing bond
            return [[-x for x in d0]][:n_new]
        if n_total == 3:
            # Trigonal planar: two new H's at 120 deg from d0 and from
            # each other, in an arbitrary plane containing d0.
            perp = _arbitrary_perp(d0)
            ang = math.radians(120)
            d1 = _unit(_add([c*math.cos(ang) for c in d0],
                             perp, math.sin(ang)))
            d2 = _unit(_add([c*math.cos(-ang) for c in d0],
                             perp, math.sin(-ang)))
            return [d1, d2][:n_new]
        if n_total == 4:
            # Tetrahedral: three new H's arranged symmetrically around d0
            perp = _arbitrary_perp(d0)
            perp2 = _unit(_cross(d0, perp))
            tet_ang = math.radians(109.5)
            dirs = []
            for k in range(3):
                phi = math.radians(120 * k)
                radial = _add([c*math.cos(phi) for c in perp],
                               perp2, math.sin(phi))
                d = _unit(_add([c*math.cos(tet_ang) for c in d0],
                                radial, math.sin(tet_ang)))
                dirs.append(d)
            return dirs[:n_new]

    if n_exist == 2:
        d0, d1 = existing_dirs[0], existing_dirs[1]
        bisector = _unit(_add(d0, d1))
        normal = _unit(_cross(d0, d1))
        if n_total == 3:
            # Trigonal planar: remaining direction is in-plane, opposite
            # the bisector of the two existing bonds.
            new_dir = [-x for x in bisector]
            return [new_dir][:n_new]
        if n_total == 4:
            # Tetrahedral: two new H's placed symmetrically out of the
            # plane defined by the two existing bonds.
            opp = [-x for x in bisector]
            tet_half = math.radians(54.75)
            d_a = _unit(_add([c*math.cos(tet_half) for c in opp],
                              normal, math.sin(tet_half)))
            d_b = _unit(_add([c*math.cos(tet_half) for c in opp],
                              normal, -math.sin(tet_half)))
            return [d_a, d_b][:n_new]

    if n_exist == 3 and n_total == 4:
        # Tetrahedral: remaining direction is opposite the sum of the
        # three existing bond directions.
        s = [0.0, 0.0, 0.0]
        for d in existing_dirs:
            s = _add(s, d)
        new_dir = _unit([-x for x in s])
        return [new_dir][:n_new]

    # Fallback for unusual coordination numbers: spread remaining
    # directions roughly away from the average of existing bonds.
    s = [0.0, 0.0, 0.0]
    for d in existing_dirs:
        s = _add(s, d)
    avg = _unit(s) if any(s) else [0.0, 0.0, 1.0]
    base = [-x for x in avg]
    perp = _arbitrary_perp(base)
    dirs = []
    for k in range(n_new):
        ang = math.radians(360 * k / max(n_new, 1))
        d = _unit(_add([c*math.cos(ang*0+1) for c in base], perp,
                        0.3 * math.sin(ang)))
        dirs.append(d)
    return dirs


def add_missing_hydrogens(mol, bond_length=1.09):
    """
    Return a new molecule dict with implicit hydrogens added based on
    standard valence rules, with 3D positions placed using simple
    geometric heuristics from each atom's existing neighbors.

    Atoms that already satisfy their standard valence (including
    correctly-placed existing hydrogens) are left untouched. Charged
    atoms recorded via `mol.get('charge')` (1-based dict, atom index ->
    formal charge) are accounted for when computing target valence
    (e.g. a +1 nitrogen gets one more bond than neutral nitrogen).

    Returns a new dict; does not mutate the input.
    """
    natom = mol['natom']
    atom = list(mol['atom'])
    nuc = list(mol['nuc'])
    coords = [list(c) for c in mol['coords']]
    idcon = [list(row) for row in mol['idcon']]
    ival = list(mol['ival'])
    nbo = SparseMatrix()
    charge = mol.get('charge', {})

    # Copy existing bond orders into the new sparse matrix
    for i in range(1, natom + 1):
        for j in range(1, ival[i] + 1):
            nbr = idcon[i][j]
            bo = mol['nbo'][i][nbr]
            nbo[i][nbr] = bo if bo else 1

    next_idx = natom

    for i in range(1, natom + 1):
        elem = nuc[i]
        if elem == 1 or elem == 0:
            continue
        used = _bond_order_sum(mol, i)
        target = _best_standard_valence(elem, used)
        if target is None:
            continue
        q = charge.get(i, 0)
        # A positive formal charge on N/O-type centers typically means
        # one extra bond is "normal" (e.g. ammonium); negative charge
        # means one fewer.
        target_adj = target + q
        n_missing = target_adj - used
        if n_missing <= 0:
            continue
        if ival[i] + n_missing > 4:
            # Respect the same 4-neighbor cap used elsewhere in the
            # connectivity tables (idcon is sized for max 4 neighbors).
            n_missing = 4 - ival[i]
            if n_missing <= 0:
                continue

        neighbor_positions = [coords[idcon[i][j]] for j in range(1, ival[i] + 1)]
        dirs = _new_h_directions(coords[i], neighbor_positions, n_missing,
                                  bond_length=bond_length)

        for d in dirs:
            next_idx += 1
            if next_idx >= len(atom):
                # Grow arrays if we exceed the original MAX allocation
                grow = next_idx + 100
                atom.extend([None] * (grow - len(atom) + 1))
                nuc.extend([0] * (grow - len(nuc) + 1))
                coords.extend([[0.0, 0.0, 0.0] for _ in range(grow - len(coords) + 1)])
                idcon.extend([[0]*5 for _ in range(grow - len(idcon) + 1)])
                ival.extend([0] * (grow - len(ival) + 1))

            atom[next_idx] = 'H '
            nuc[next_idx] = 1
            new_pos = _add(coords[i], d, bond_length)
            coords[next_idx][0], coords[next_idx][1], coords[next_idx][2] = new_pos
            ival[next_idx] = 1
            idcon[next_idx][1] = i

            ival[i] += 1
            idcon[i][ival[i]] = next_idx

            nbo[i][next_idx] = 1
            nbo[next_idx][i] = 1

    new_natom = next_idx

    # Rebuild icon / isum over the (possibly larger) atom set
    icon = [[0]*5 for _ in range(new_natom + 1)]
    isum = [0] * (new_natom + 1)
    for i in range(1, new_natom + 1):
        s = 0
        for j in range(1, ival[i] + 1):
            nb = nuc[idcon[i][j]]
            icon[i][j] = nb
            s += nb
        isum[i] = s

    numbond = sum(ival[i] for i in range(1, new_natom + 1)) // 2

    result = dict(mol)
    result.update(
        natom=new_natom,
        numbond=numbond,
        atom=atom,
        nuc=nuc,
        coords=coords,
        idcon=idcon,
        ival=ival,
        nbo=nbo,
        icon=icon,
        isum=isum,
        hydrogens_added=new_natom - natom,
    )
    return result
