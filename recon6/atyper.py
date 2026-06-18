"""Translation of the ``atyper`` and ``atomsort2`` subroutines from
recon5-5.f.  Produces the 49-character "atom type" code string used to
look up TAE descriptor data.
"""


def atomsort2(numval, iatomnum, ivalence, imax, nnay):
    """In-place sort of 4-element lists, mirroring atomsort2.

    All four lists are length 4; ``numval`` is the number of valid
    leading entries (<=4).
    """
    # Bubble sort by atomic number (descending)
    for i in range(numval - 1):
        for j in range(i + 1, numval):
            if iatomnum[i] < iatomnum[j]:
                iatomnum[i], iatomnum[j] = iatomnum[j], iatomnum[i]
                ivalence[i], ivalence[j] = ivalence[j], ivalence[i]
                imax[i], imax[j] = imax[j], imax[i]
                nnay[i], nnay[j] = nnay[j], nnay[i]

    # Then sort by valence (descending) among entries with equal atomic number
    for i in range(3):
        for j in range(i + 1, 4):
            if iatomnum[i] == iatomnum[j] and ivalence[i] < ivalence[j]:
                ivalence[i], ivalence[j] = ivalence[j], ivalence[i]
                imax[i], imax[j] = imax[j], imax[i]
                nnay[i], nnay[j] = nnay[j], nnay[i]


def atyper(numval, idcon, icon, isize, nuc, nselect):
    """
    Determine the atom-type code string for atom ``nselect``.

    Parameters (all 1-based dicts/lists, index 0 unused):
        numval : valence per atom
        idcon  : idcon[i][j] -> connected atom index (j=1..numval[i])
        icon   : icon[i][j]  -> atomic number of connected atom j
        isize  : ring-size code per atom (from ringid)
        nuc    : atomic number per atom
        nselect : index of the atom to type

    Returns
    -------
    temp : str, 49-character atom type code (digits 0-9, '0'-padded)
    """
    if numval[nselect] == 1:
        n = nuc[nselect]
        if n in (1, 8, 9, 17, 16, 35, 53):
            i = idcon[nselect][1]
        else:
            i = nselect
    else:
        i = nselect

    natom = [0, 0, 0, 0]
    nval = [0, 0, 0, 0]
    isum = [0, 0, 0, 0]
    idummy = [0, 0, 0, 0]

    nnNUC1 = [0, 0, 0, 0]
    nnNUC2 = [0, 0, 0, 0]
    nnNUC3 = [0, 0, 0, 0]
    nnNUC4 = [0, 0, 0, 0]

    for j in range(1, numval[i] + 1):
        naynum = idcon[i][j]
        natom[j - 1] = icon[i][j]
        nval[j - 1] = numval[naynum]
        idummy[j - 1] = j

        nnNUC = [0, 0, 0, 0]
        for m in range(1, nval[j - 1] + 1):
            nnNUC[m - 1] = icon[naynum][m]

        # Bubble sort nnNUC (descending) over its first nval[j-1] entries
        nv = nval[j - 1]
        for m in range(nv - 1):
            for n in range(m + 1, nv):
                if nnNUC[m] < nnNUC[n]:
                    nnNUC[m], nnNUC[n] = nnNUC[n], nnNUC[m]

        for m in range(1, nv + 1):
            if m == 1:
                nnNUC1[j - 1] = nnNUC[m - 1]
            elif m == 2:
                nnNUC2[j - 1] = nnNUC[m - 1]
            elif m == 3:
                nnNUC3[j - 1] = nnNUC[m - 1]
            elif m == 4:
                nnNUC4[j - 1] = nnNUC[m - 1]

    m = numval[i]
    atomsort2(m, natom, nval, isum, idummy)

    id1, id2, id3, id4 = (idummy[0] - 1, idummy[1] - 1, idummy[2] - 1, idummy[3] - 1)

    # Fortran format 7777: i2,i1,4i2,4i1,i2,16i2
    def fmt_i(val, width):
        s = "%d" % val
        if len(s) > width:
            return s[-width:]
        return s.rjust(width)

    parts = []
    parts.append(fmt_i(nuc[i], 2))
    parts.append(fmt_i(numval[i], 1))
    for k in range(4):
        parts.append(fmt_i(natom[k], 2))
    for k in range(4):
        parts.append(fmt_i(nval[k], 1))
    parts.append(fmt_i(isize[i], 2))

    def nn_block(idx):
        return [nnNUC1[idx], nnNUC2[idx], nnNUC3[idx], nnNUC4[idx]]

    for idx in (id1, id2, id3, id4):
        for v in nn_block(idx):
            parts.append(fmt_i(v, 2))

    temp = "".join(parts)
    # temp should be exactly 49 chars: 2+1+8+4+2+32 = 49
    temp = (temp + " " * 49)[:49]
    # Replace blanks with '0' (Fortran: any blank char in field -> '0')
    temp = temp.replace(" ", "0")
    return temp
