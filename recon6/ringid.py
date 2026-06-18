"""Translation of ringid.f.

Determines, for each atom, whether it is part of a ring system and
records the ring size (or fused-ring code) in ``isize``.

All arrays here use 1-based atom indexing to mirror the Fortran source
as closely as possible; index 0 is unused.
"""


def ringid(numatom, numval, idcon):
    """
    Parameters
    ----------
    numatom : int
    numval : list[int]   1-based, length >= numatom+1, valence per atom
    idcon : list[list[int]]  1-based both dims; idcon[i][j] = atom index
        connected to atom i via its j-th bond (j=1..numval[i])

    Returns
    -------
    isize : list[int]  1-based, length numatom+1; ring-size code per atom
    """
    isize = [0] * (numatom + 1)

    for i in range(1, numatom + 1):
        ihit3 = ihit4 = ihit5 = ihit6 = 0
        isize[i] = 0

        if numval[i] < 2:
            continue

        # NODE 1
        node1 = []
        for j in range(1, numval[i] + 1):
            cand = idcon[i][j]
            if numval[cand] >= 2:
                node1.append(cand)

        node2 = []
        ip2 = []
        if node1:
            for a in node1:
                for k in range(1, numval[a] + 1):
                    cand = idcon[a][k]
                    if numval[cand] >= 2 and cand != i:
                        node2.append(cand)
                        ip2.append(a)

        node3 = []
        ip3 = []
        if node2:
            for j, a in enumerate(node2):
                parent = ip2[j]
                for k in range(1, numval[a] + 1):
                    cand = idcon[a][k]
                    if numval[cand] >= 2 and cand != parent:
                        node3.append(cand)
                        ip3.append(a)
                        if cand == i:
                            if ihit3 == 0:
                                ihit3 = 1
                            elif ihit3 == 1:
                                ihit3 = 2
                            elif ihit3 == 2:
                                ihit3 = 3
                            elif ihit3 == 3:
                                ihit3 = 4

        node4 = []
        ip4 = []
        if node3:
            for j, a in enumerate(node3):
                parent = ip3[j]
                for k in range(1, numval[a] + 1):
                    cand = idcon[a][k]
                    if numval[cand] >= 2 and cand != parent:
                        node4.append(cand)
                        if cand == i:
                            if ihit4 == 0:
                                ihit4 = 1
                            elif ihit4 == 1:
                                ihit4 = 2
                            elif ihit4 == 2:
                                ihit4 = 3
                            elif ihit4 == 3:
                                ihit4 = 4
                        ip4.append(a)

        node5 = []
        ip5 = []
        if node4:
            for j, a in enumerate(node4):
                parent = ip4[j]
                for k in range(1, numval[a] + 1):
                    cand = idcon[a][k]
                    if numval[cand] >= 2 and cand != parent:
                        node5.append(cand)
                        if cand == i:
                            if ihit5 == 0:
                                ihit5 = 1
                            elif ihit5 == 1:
                                ihit5 = 2
                            elif ihit5 == 2:
                                ihit5 = 3
                            elif ihit5 == 3:
                                ihit5 = 4
                        ip5.append(a)

        if node5:
            for j, a in enumerate(node5):
                parent = ip5[j]
                for k in range(1, numval[a] + 1):
                    cand = idcon[a][k]
                    if numval[cand] >= 2 and cand != parent:
                        if cand == i:
                            if ihit6 == 0:
                                ihit6 = 1
                            elif ihit6 == 1:
                                ihit6 = 2
                            elif ihit6 == 2:
                                ihit6 = 3
                            elif ihit6 == 3:
                                ihit6 = 4

        if ihit3 == 2 and ihit4 == 0 and ihit5 == 0 and ihit6 == 0:
            isize[i] = 3
        elif ihit3 == 4:
            isize[i] = 33
        elif ihit3 == 0 and ihit4 == 2 and ihit5 == 0 and ihit6 == 0:
            isize[i] = 4
        elif ihit4 == 4:
            isize[i] = 44
        elif ihit3 == 0 and ihit4 == 0 and ihit5 == 2 and ihit6 == 0:
            isize[i] = 5
        elif ihit5 == 4:
            isize[i] = 55
        elif ihit3 == 0 and ihit4 == 0 and ihit5 == 0 and ihit6 == 2:
            isize[i] = 6
        elif ihit6 == 4:
            isize[i] = 66
        elif ihit3 == 2 and ihit4 == 2 and ihit5 == 0 and ihit6 == 0:
            isize[i] = 34
        elif ihit3 == 0 and ihit4 == 2 and ihit5 == 2 and ihit6 == 0:
            isize[i] = 45
        elif ihit3 == 0 and ihit4 == 0 and ihit5 == 2 and ihit6 == 2:
            isize[i] = 56

    return isize
