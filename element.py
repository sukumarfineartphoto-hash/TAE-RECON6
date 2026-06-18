"""Translation of ele.f - normalizes raw atom-name strings (as found in
PDB/MOL2 files) into 2-character element symbols matching the periodic
table convention used elsewhere in RECON5."""


def ele(atm):
    """Translate a raw atom name (e.g. from a PDB/MOL2 file) into a
    2-character element symbol, matching ele.f's logic.

    ``atm`` should be a string of at least 2 characters (will be
    space-padded if shorter).
    """
    temp = (atm + "  ")[:2]
    c1 = temp[0]
    c2 = temp[1]

    if c1 == 'H':
        return 'H '
    elif c1 == 'N':
        if c2 in ('a', 'A'):
            return 'NA'
        elif c2 in ('b', 'B'):
            return 'NB'
        elif c2 in ('i', 'I'):
            return 'NI'
        else:
            return 'N '
    elif c1 == 'S':
        if c2 in ('i', 'I'):
            return 'SI'
        elif c2 in ('b', 'B'):
            return 'SB'
        elif c2 in ('c', 'C'):
            return 'SC'
        elif c2 in ('e', 'E'):
            return 'SE'
        elif c2 in ('n', 'N'):
            return 'SN'
        elif c2 in ('r', 'R'):
            return 'SR'
        else:
            return 'S '
    elif c1 == 'O':
        return 'O '
    elif c1 == 'F':
        if c2 in ('e', 'E'):
            return 'FE'
        else:
            return 'F '
    elif c1 == 'P':
        if c2 in ('d', 'D'):
            return 'PD'
        else:
            return 'P '
    elif c1 == 'C':
        if c2 in ('l', 'L'):
            return 'CL'
        elif c2 in ('a', 'A'):
            return 'CA'
        elif c2 in ('d', 'D'):
            return 'CD'
        elif c2 in ('o', 'O'):
            return 'CO'
        elif c2 in ('r', 'R'):
            return 'CR'
        elif c2 in ('u', 'U'):
            return 'CU'
        else:
            return 'C '
    elif c1 == 'B':
        if c2 in ('r', 'R'):
            return 'BR'
        elif c2 in ('e', 'E'):
            return 'BE'
        else:
            return 'B '
    elif c1 == 'I':
        if c2 in ('n', 'N'):
            return 'IN'
        else:
            return 'I '
    elif c1 == 'V':
        return 'V '
    elif c1 == 'Y':
        return 'Y '
    else:
        return temp
