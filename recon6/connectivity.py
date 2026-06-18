"""
Distance-based connectivity determination and bond-order inference,
shared by any reader that has 3D atomic coordinates but no explicit
bond list (PDB without CONECT records, Gaussian .com Cartesian input,
etc.).

Given only atom identities and positions, this:
  1. Determines which atoms are bonded, using a per-element-pair
     reference bond length table and a tolerance factor (mirrors the
     distance criterion in the original pdbread.f).
  2. Infers likely double-bond character for terminal (single-
     heavy-neighbor) N/O atoms and sulfonyl/sulfate sulfur centers,
     using bond-length shortening as the signal - since pure distance/
     connectivity data carries no explicit bond order, but a
     measurably shortened bond is a reliable, well-established
     indicator of unsaturation (e.g. C=O is consistently ~0.85-0.89x
     the single-bond reference length; S=O is consistently shortened
     but by a smaller margin, requiring its own threshold).

This bond-order inference only matters for standard-valence hydrogen
counting (see hydrogenate.py) - a terminal O or N that's actually
double-bonded would otherwise look under-valent and incorrectly
receive a spurious extra hydrogen.
"""
import math

# A general bond-length-ratio threshold below which a terminal N or O
# bond is inferred to be a double bond (C=O ~0.85-0.89x single-bond
# reference; C=N similar). Single bonds, including slightly-short ones
# like some C-OH centers, measure ~0.95-1.0x and are left alone.
_DOUBLE_BOND_RATIO_THRESHOLD = 0.93

# S=O bonds (sulfonyl/sulfate) are shortened by a smaller margin than
# C=O, so they need their own, looser threshold plus the additional
# structural requirement (>=2 terminal O's on the same S) to avoid
# false positives from the general threshold above.
_SULFONYL_RATIO_THRESHOLD = 0.97


def determine_connectivity(natom, nuc, coords, bondl, idcon, ival, nbo,
                            tolerance=1.1, infer_bond_orders=True):
    """
    Populate idcon/ival/nbo in place from atom positions and a
    reference bond-length table, then (optionally) infer double-bond
    character for terminal N/O atoms and sulfonyl sulfur centers.

    Parameters
    ----------
    natom : int
    nuc : list[int], 1-based, atomic number per atom
    coords : list[[x,y,z]], 1-based
    bondl : 55x55 reference single-bond-length table (1-based)
    idcon, ival, nbo : output structures to populate (same shapes as
        used throughout the readers package - idcon is list-of-5-lists,
        ival is a flat list, nbo is a SparseMatrix or compatible object)
    tolerance : bonded if measured distance <= tolerance * reference length
    infer_bond_orders : whether to run the double-bond inference step

    Returns
    -------
    numbond : int, total number of bonds found
    """
    numbond = 0
    seen = [[False] * (natom + 1) for _ in range(natom + 1)]
    measured_dist = {}

    for i in range(1, natom + 1):
        cnt = 0
        for j in range(1, natom + 1):
            if i == j:
                continue
            ni, nj = nuc[i], nuc[j]
            if ni > 54 or nj > 54 or ni == 0 or nj == 0:
                continue
            bl = bondl[ni][nj]
            if bl == 0.0:
                continue
            dx = coords[i][0] - coords[j][0]
            dy = coords[i][1] - coords[j][1]
            dz = coords[i][2] - coords[j][2]
            dist = math.sqrt(dx*dx + dy*dy + dz*dz)
            if dist > tolerance * bl:
                continue
            cnt += 1
            if cnt <= 4:
                idcon[i][cnt] = j
            if not seen[i][j] and not seen[j][i]:
                numbond += 1
                seen[i][j] = seen[j][i] = True
                measured_dist[(i, j)] = (dist, bl)
        ival[i] = cnt

    if infer_bond_orders:
        for i in range(1, natom + 1):
            if ival[i] != 1:
                continue
            if nuc[i] not in (7, 8):  # only N/O, where this matters for H-counting
                continue
            j = idcon[i][1]
            pair = (i, j) if (i, j) in measured_dist else (j, i)
            if pair not in measured_dist:
                continue
            dist, bl = measured_dist[pair]
            if bl > 0 and (dist / bl) < _DOUBLE_BOND_RATIO_THRESHOLD:
                nbo[i][j] = 2
                nbo[j][i] = 2

        for i in range(1, natom + 1):
            if nuc[i] != 16:
                continue
            terminal_o = []
            for k in range(1, ival[i] + 1):
                j = idcon[i][k]
                if nuc[j] == 8 and ival[j] == 1:
                    pair = (i, j) if (i, j) in measured_dist else (j, i)
                    if pair in measured_dist:
                        dist, bl = measured_dist[pair]
                        if bl > 0 and (dist / bl) < _SULFONYL_RATIO_THRESHOLD:
                            terminal_o.append(j)
            if len(terminal_o) >= 2:
                for j in terminal_o:
                    nbo[i][j] = 2
                    nbo[j][i] = 2

    return numbond
