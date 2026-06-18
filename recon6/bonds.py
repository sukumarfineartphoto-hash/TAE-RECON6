"""Bond length table loader."""


def load_bond_table(bond_file):
    """Read a RECON bond-length file and return a 55x55 table (1-based).
    Format per line: int3 int3 6x float9.6  (handles Windows line endings)
    """
    bondl = [[0.0] * 55 for _ in range(55)]
    with open(bond_file, 'r', errors='replace') as fh:
        for line in fh:
            line = line.rstrip('\r\n')
            if len(line) < 12:
                continue
            try:
                n1 = int(line[0:3])
                n2 = int(line[3:6])
                # field starts at col 6+6=12 (format: 2I3,6X,F9.6)
                val = float(line[12:].strip().split()[0])
            except (ValueError, IndexError):
                continue
            if 1 <= n1 <= 54 and 1 <= n2 <= 54:
                bondl[n1][n2] = val
                if n1 != n2:
                    bondl[n2][n1] = val
    return bondl
