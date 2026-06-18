"""Periodic table data translated from the DATA TABLE/TABLECAP arrays
in recon5-5.f (elements 1-54, H through Xe)."""

# Element symbols indexed 1..54 (index 0 unused), mixed case as in source
TABLE = [None,
    'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne',
    'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca', 'Sc', 'Ti',
    'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As',
    'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru',
    'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I ', 'Xe',
]

# Upper-case variants
TABLECAP = [None,
    'H ', 'HE', 'LI', 'BE', 'B ', 'C ', 'N ', 'O ', 'F ', 'NE',
    'NA', 'MG', 'AL', 'SI', 'P ', 'S ', 'CL', 'AR', 'K ', 'CA', 'SC', 'TI',
    'V ', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN', 'GA', 'GE', 'AS',
    'SE', 'BR', 'KR', 'RB', 'SR', 'Y ', 'ZR', 'NB', 'MO', 'TC', 'RU',
    'RH', 'PD', 'AG', 'CD', 'IN', 'SN', 'SB', 'TE', 'I ', 'XE',
]


def atomic_number(symbol):
    """Return atomic number (1-54) for a 2-char element symbol, matching
    either case variant, or None if not found."""
    for i in range(1, 55):
        if symbol == TABLE[i] or symbol == TABLECAP[i]:
            return i
    return None
