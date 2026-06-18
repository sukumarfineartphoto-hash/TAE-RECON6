"""TAE database index (TAE.LIST equivalent) and the ``gettae`` matching
algorithm from recon5-5.f, plus ``build_modtype``.
"""
import os


class TaeIndex:
    """Indexes .dat files from the DATA directory, mirroring TAE.LIST."""
    def __init__(self, data_dir):
        self.data_dir = data_dir
        self.entries = sorted(
            fn for fn in os.listdir(data_dir) if fn.lower().endswith(".dat")
        )
        self.num = len(self.entries)

    def path(self, basename):
        return os.path.join(self.data_dir, basename)


def _match_level(temp, temp2, k):
    """Compute match level (-3..3) between modtype string and a TAE entry name,
    mirroring the fall-through goto-9 structure of gettae in recon5-5.f."""
    level = -3
    if k != 0:
        if temp[0:k] != temp2[0:k]:
            return level
        level = -2
    else:
        level = -2

    if temp[k:k+2] != temp2[k:k+2]:
        return level
    level = -1

    if temp[k+2:k+3] != temp2[k+2:k+3]:
        return level
    level = 0

    if temp[k+3:k+11] != temp2[k+3:k+11]:
        return level
    level = 1

    if temp[k+11:k+15] != temp2[k+11:k+15]:
        return level
    level = 2

    if temp[k+15:k+17] != temp2[k+15:k+17]:
        return level
    level = 3

    return level


def build_modtype(numatom, ival, idcon, icon, isize, nuc):
    """Build 51-char modtype strings for all atoms (main program lines ~395-432)."""
    from .atyper import atyper
    modtype = [None] * (numatom + 1)
    for i in range(1, numatom + 1):
        temp = atyper(ival, idcon, icon, isize, nuc, i)  # 49 chars
        n = nuc[i]
        if ival[i] == 1:
            if n == 1:
                temp2 = 'H' + temp + ' '
            elif n == 8:
                temp2 = 'O' + temp + ' '
            elif n == 9:
                temp2 = 'F' + temp + ' '
            elif n == 16:
                temp2 = 'S' + temp + ' '
            elif n == 17:
                temp2 = 'CL' + temp
            elif n == 35:
                temp2 = 'BR' + temp
            elif n == 53:
                temp2 = 'I' + temp + ' '
            else:
                temp2 = temp + '  '
        else:
            temp2 = temp + '  '
        modtype[i] = (temp2 + ' ' * 51)[:51]
    return modtype


def gettae(numatom, modtype, tae_index, warnings=None):
    """Match each atom's modtype to a TAE .dat file (translation of gettae).

    Returns (atomtype_list, lev_list), both 1-based.
    atomtype[i] is the .dat filename (basename) for atom i.
    lev[i] is the match level (-3..3).
    """
    entries = tae_index.entries
    num = tae_index.num
    atomtype = [None] * (numatom + 1)
    lev = [0] * (numatom + 1)

    for i in range(1, numatom + 1):
        temp = modtype[i]
        if temp.strip().rstrip('0').strip() == '' or len(temp.strip()) == 0:
            # Degenerate/blank atom-type code, typically caused by an
            # atom with zero bonds (disconnected atom in the input
            # structure). The TAE match below will be essentially
            # arbitrary; this is a data-quality issue in the input,
            # not a bug, but it's worth surfacing.
            if warnings is not None:
                warnings.append(
                    "Atom %d has no bonds (disconnected); TAE match will "
                    "be unreliable for this atom." % i
                )
        c0 = temp[0]
        if c0 in ('H', 'F', 'O', 'S', 'I'):
            k = 1
        elif c0 in ('C', 'B'):
            k = 2
        else:
            k = 0

        last = -3
        best_entry = entries[num - 1] if num else ''
        level = -3

        for j in range(num):
            temp2 = entries[j]
            lv = _match_level(temp, temp2, k)
            if lv == 3:
                best_entry = temp2
                level = lv
                break
            if lv > last:
                last = lv
                best_entry = temp2
                level = lv

        ndot = best_entry.find('.')
        if ndot < 0:
            ndot = len(best_entry)
        atomtype[i] = best_entry[:ndot] + '.dat'
        lev[i] = level

    return atomtype, lev
