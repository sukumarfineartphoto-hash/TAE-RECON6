"""SDF (MDL SD file) reader - translation of rsdf.f."""
from ..periodic import atomic_number
from ..sparse import SparseMatrix


def _read_molecule_block(fileobj):
    """Read raw lines for one molecule record, from just after the
    previous '$$$$' (or start of file) through and including the next
    '$$$$'. Returns the list of lines (without the '$$$$' line itself),
    or raises StopIteration if no more molecule blocks remain.

    Buffering the whole block before parsing means a malformed record
    can't desync the file position for subsequent molecules - the
    caller always lands cleanly on the next record boundary.
    """
    lines = []
    saw_any_content = False
    while True:
        line = fileobj.readline()
        if not line:
            if saw_any_content:
                # Trailing block with no '$$$$' terminator (malformed
                # end of file) - return what we have.
                return lines
            raise StopIteration
        if line.startswith('$$$$'):
            if not saw_any_content:
                # Stray/duplicate '$$$$' with nothing before it - skip
                # and keep looking for the next real block.
                continue
            return lines
        if line.strip():
            saw_any_content = True
        lines.append(line)


def read_sdf(fileobj, molnum=1):
    """Read next molecule from an open SDF fileobj.
    Returns a dict with keys:
        natom, numbond, atom, nuc, coords, idcon, ival, nbo, icon, isum
    or raises StopIteration at end of file.
    All arrays are 1-based (index 0 unused/zero).
    """
    lines = _read_molecule_block(fileobj)
    return _parse_molecule_block(lines)


def _parse_molecule_block(lines):
    MAX = 1000
    atom = [None] * (MAX + 1)
    nuc = [0] * (MAX + 1)
    coords = [[0.0, 0.0, 0.0] for _ in range(MAX + 1)]
    idcon = [[0] * 5 for _ in range(MAX + 1)]
    ival = [0] * (MAX + 1)
    nbo = SparseMatrix()
    icon = [[0] * 5 for _ in range(MAX + 1)]
    isum = [0] * (MAX + 1)
    jch = [0] * (MAX + 1)

    pos = 0
    n = len(lines)

    def next_line():
        nonlocal pos
        if pos >= n:
            return ''
        line = lines[pos]
        pos += 1
        return line

    title_line = next_line()
    mol_name = title_line.strip()
    next_line()  # program/timestamp line
    next_line()  # comment line (often blank)

    counts = next_line()
    if not counts or not counts[0:3].strip():
        raise ValueError("Missing or malformed counts line")
    natom = int(counts[0:3])
    numbond = int(counts[3:6])

    # Atom block
    for i in range(1, natom + 1):
        line = next_line()
        coords[i][0] = float(line[0:10])
        coords[i][1] = float(line[10:20])
        coords[i][2] = float(line[20:30])
        sym = line[31:33].strip()
        atom[i] = (sym + '  ')[:2]

    # Bond block
    for _ in range(numbond):
        line = next_line()
        iat = int(line[0:3])
        jat = int(line[3:6])
        ibo = int(line[6:9])
        nbo[iat][jat] = ibo
        nbo[jat][iat] = ibo
        if ival[iat] < 4 and ival[jat] < 4:
            ival[iat] += 1
            ival[jat] += 1
            idcon[iat][ival[iat]] = jat
            idcon[jat][ival[jat]] = iat

    # Assign atomic numbers
    for i in range(1, natom + 1):
        an = atomic_number(atom[i])
        if an is not None:
            nuc[i] = an

    # Remaining lines (M  CHG, M  END, data items). Also capture
    # "> <TAG>" data-block fields (the SDF "properties" block),
    # preserving field order, for downstream merging into the combined
    # output CSV alongside computed descriptors.
    data_fields = {}
    current_tag = None
    while pos < n:
        line = next_line()
        if line.startswith('>'):
            # Tag line, e.g. "> <LIC50>" or ">  <EXTREG>  (some name)"
            tag = line.strip().lstrip('>').strip()
            tag = tag.replace('<', '').replace('>', '').strip()
            # If the tag line also has trailing text after the closing
            # '<...>' (as in ">  <EXTREG>  (name)"), keep only the
            # bracketed part as the tag name; the parenthetical text is
            # SDF-author commentary, not part of the field's identity.
            if '(' in tag:
                tag = tag.split('(')[0].strip()
            current_tag = tag
            if current_tag not in data_fields:
                data_fields[current_tag] = []
            continue
        stripped = line.rstrip('\r\n')
        if stripped.strip().startswith('M  CHG'):
            parts = stripped.split()
            nch = int(parts[2])
            for k in range(nch):
                idx = int(parts[3 + 2 * k])
                chg = int(parts[4 + 2 * k])
                jch[idx] = chg
            current_tag = None
            continue
        if stripped.strip() == 'M  END':
            current_tag = None
            continue
        if current_tag is not None:
            if stripped.strip() == '':
                # Blank line ends the current tag's value block
                current_tag = None
            else:
                data_fields[current_tag].append(stripped)

    # Collapse each tag's value lines into a single string (joined by
    # space if a tag had multiple lines, matching the common case of
    # one value per tag while still tolerating multi-line values).
    data_fields = {k: ' '.join(v) for k, v in data_fields.items()}

    # Build icon / isum
    for i in range(1, natom + 1):
        s = 0
        for j in range(1, ival[i] + 1):
            nb = nuc[idcon[i][j]]
            icon[i][j] = nb
            s += nb
        isum[i] = s

    charge = {i: jch[i] for i in range(1, natom + 1) if jch[i] != 0}

    return dict(natom=natom, numbond=numbond, atom=atom, nuc=nuc,
                coords=coords, idcon=idcon, ival=ival, nbo=nbo,
                icon=icon, isum=isum, mol_name=mol_name, charge=charge,
                data_fields=data_fields)


def read_sdf_data_fields_only(lines):
    """Best-effort extraction of just the title line and '> <TAG>'
    data-block fields from a raw molecule block, used when full
    structural parsing failed but we still want to preserve the
    record's response/property data for positional alignment."""
    if not lines:
        return None, {}
    mol_name = lines[0].strip()
    data_fields = {}
    current_tag = None
    in_data_block = False
    for line in lines:
        if line.strip() == 'M  END':
            in_data_block = True
            continue
        if not in_data_block:
            continue
        if line.startswith('>'):
            tag = line.strip().lstrip('>').strip()
            tag = tag.replace('<', '').replace('>', '').strip()
            if '(' in tag:
                tag = tag.split('(')[0].strip()
            current_tag = tag
            if current_tag not in data_fields:
                data_fields[current_tag] = []
            continue
        stripped = line.rstrip('\r\n')
        if current_tag is not None:
            if stripped.strip() == '':
                current_tag = None
            else:
                data_fields[current_tag].append(stripped)
    data_fields = {k: ' '.join(v) for k, v in data_fields.items()}
    return mol_name, data_fields
