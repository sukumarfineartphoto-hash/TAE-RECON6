"""
Descriptor accumulation (translation of the qmf subroutine from recon5-5.f).

Reads per-atom TAE .dat records, accumulates molecular totals, and returns
a flat dict of all descriptor values matching the recon.ff column names.
"""
from .taedat import read_tae_dat

# Scalar groups used in the accumulation loop
_HIST_SUFFIXES = ['1','2','3','4','5','6','7','8','9','10']
_HIST_GROUPS = [
    # (total_prefix, dat_prefix, min_key, max_key, ia_key, count)
    ('TSIDRN',  'SIDRN',  'TDRNMn',  'TDRNMx',  None,     'drna',  'TDRNA'),
    ('TSIDKN',  'SIDKN',  'TDKMn',   'TDKMx',   None,     'dkna',  'TDKNA'),
    ('TSIK',    'SIK',    'TSIKMn',  'TSIKMx',  None,     'sika',  'TSIKA'),
    ('TSIDGN',  'SIDGN',  'TDGNMn',  'TDGNMx',  None,     'dgna',  'TDGNA'),
    ('TSIG',    'SIG',    'TSIGMn',  'TSIGMx',  None,     'siga',  'TSIGA'),
    ('TSIEP',   'SIEP',   'TSIEPMn', 'TSIEPMx', None,     'siepa', 'TSIEPA'),
]


def compute_descriptors(mol, atomtype_list, tae_index, lflag=True):
    """
    Parameters
    ----------
    mol : dict from a reader (natom, ival, idcon, nuc, atom, nbo, coords, ...)
    atomtype_list : list[str], 1-based, .dat basenames from gettae
    tae_index : TaeIndex
    lflag : bool; True means skip geom (multipole) loading (icheck != 2)

    Returns dict with all accumulated molecular descriptors.
    """
    natom = mol['natom']
    ival = mol['ival']
    idcon = mol['idcon']
    nuc = mol['nuc']
    nbo = mol['nbo']

    # Molecule-level accumulators
    AE = 0.0; Epop = 0.0; Volume = 0.0; TSArea = 0.0
    TSIDRN = 0.0; TDRNMn = 20.0; TDRNMx = -20.0; TDRNA = [0.0]*10
    TSIDKN = 0.0; TDKMn  = 20.0; TDKMx  = -20.0; TDKNA = [0.0]*10
    TSIK   = 0.0; TSIKMn = 20.0; TSIKMx = -20.0; TSIKA = [0.0]*10
    TSIDGN = 0.0; TDGNMn = 20.0; TDGNMx = -20.0; TDGNA = [0.0]*10
    TSIG   = 0.0; TSIGMn = 20.0; TSIGMx = -20.0; TSIGA = [0.0]*10
    TSIEP  = 0.0; TSIEPMn= 20.0; TSIEPMx= -20.0; TSIEPA= [0.0]*10
    TEP = [0.0]*10
    TPIPMin = 20.0; TPIPMax = -20.0; SUMPIP = 0.0; TPIP = [0.0]*20
    TLapl = 0.0; TLaplMin = 20.0; TLaplMax = -20.0; TLaplBins = [0.0]*10
    TFuk  = 0.0; TFukMin  = 20.0; TFukMax  = -20.0; TFukBins  = [0.0]*10

    atom_records = []

    for i in range(1, natom + 1):
        dat_path = tae_index.path(atomtype_list[i])
        rec = read_tae_dat(dat_path)
        f = rec.fields
        atom_records.append(f)

        AE      += f['energy']
        Epop    += f['pop']
        Volume  += f['vol']
        TSArea  += f['sarea']

        TSIDRN += f['sidrn']
        if f['drnmn'] < TDRNMn: TDRNMn = f['drnmn']
        if f['drnmx'] > TDRNMx: TDRNMx = f['drnmx']
        for k in range(10): TDRNA[k] += f['drna%d'%(k+1)]

        TSIDKN += f['sidkn']
        if f['dkmn'] < TDKMn: TDKMn = f['dkmn']
        if f['dkmx'] > TDKMx: TDKMx = f['dkmx']
        for k in range(10): TDKNA[k] += f['dkna%d'%(k+1)]

        TSIK   += f['sik']
        if f['sikmn'] < TSIKMn: TSIKMn = f['sikmn']
        if f['sikmx'] > TSIKMx: TSIKMx = f['sikmx']
        for k in range(10): TSIKA[k] += f['sika%d'%(k+1)]

        TSIDGN += f['sidgn']
        if f['dgnmn'] < TDGNMn: TDGNMn = f['dgnmn']
        if f['dgnmx'] > TDGNMx: TDGNMx = f['dgnmx']
        for k in range(10): TDGNA[k] += f['dgna%d'%(k+1)]

        TSIG   += f['sig']
        if f['sigmn'] < TSIGMn: TSIGMn = f['sigmn']
        if f['sigmx'] > TSIGMx: TSIGMx = f['sigmx']
        for k in range(10): TSIGA[k] += f['siga%d'%(k+1)]

        TSIEP  += f['siep']
        if f['siepmn'] < TSIEPMn: TSIEPMn = f['siepmn']
        if f['siepmx'] > TSIEPMx: TSIEPMx = f['siepmx']
        for k in range(10): TSIEPA[k] += f['siepa%d'%(k+1)]

        for k in range(10): TEP[k] += f['ep%d'%(k+1)]

        if f['pipmin'] < TPIPMin: TPIPMin = f['pipmin']
        if f['pipmax'] > TPIPMax: TPIPMax = f['pipmax']
        SUMPIP += f['pipavg']
        for k in range(20): TPIP[k] += f['p%d'%(k+1)]

        TLapl += f['lapl']
        if f['laplmin'] < TLaplMin: TLaplMin = f['laplmin']
        if f['laplmax'] > TLaplMax: TLaplMax = f['laplmax']
        for k in range(10): TLaplBins[k] += f['lapl%d'%(k+1)]

        TFuk  += f['fuk']
        if f['fukmin'] < TFukMin: TFukMin = f['fukmin']
        if f['fukmax'] > TFukMax: TFukMax = f['fukmax']
        for k in range(10): TFukBins[k]  += f['fuk%d'%(k+1)]

    # Derived / normalised quantities
    if TSArea == 0.0:
        TSArea = 1.0  # guard

    TDRNIA  = TSIDRN / TSArea
    TDKIA   = TSIDKN / TSArea
    TSIKIA  = TSIK   / TSArea
    TDGNIA  = TSIDGN / TSArea
    TSIGIA  = TSIG   / TSArea
    TSIEPIA = TSIEP  / TSArea
    TLaplAVG = TLapl / TSArea
    TFukAVG  = TFuk  / TSArea
    TPIPAvg  = SUMPIP / natom

    def frac(arr): return [v / TSArea for v in arr]

    FDRNA  = frac(TDRNA)
    FDKNA  = frac(TDKNA)
    FSIKA  = frac(TSIKA)
    FDGNA  = frac(TDGNA)
    FSIGA  = frac(TSIGA)
    FEP    = frac(TEP)
    FPIP   = frac(TPIP)
    FFuk   = frac(TFukBins)
    FLapl  = frac(TLaplBins)

    # Randic chi index (H-suppressed graph)
    idegcon = [0] * (natom + 1)
    idb = [[0] * (natom + 1) for _ in range(natom + 1)]
    for i in range(1, natom + 1):
        if nuc[i] == 1: continue
        for j in range(1, ival[i] + 1):
            if nuc[idcon[i][j]] != 1:
                idegcon[i] += 1
                idb[i][idcon[i][j]] = 1
    import math
    chi = 0.0
    for i in range(2, natom + 1):
        if nuc[i] == 1: continue
        for j in range(1, i):
            if nuc[j] == 1: continue
            if idb[i][j] == 1:
                tt = idegcon[i] * idegcon[j]
                if tt > 0:
                    chi += 1.0 / math.sqrt(tt)

    # Build output dict using column names from recon.ff header
    out = dict(
        Energy=AE, Population=Epop, VOLTAE=Volume, SurfArea=TSArea,
        SIDel_RhoN=TSIDRN, Del_RhoNMin=TDRNMn, Del_RhoNMax=TDRNMx,
        Del_RhoNIA=TDRNIA,
        SIDel_KN=TSIDKN, Del_KMin=TDKMn, Del_KMax=TDKMx, Del_KIA=TDKIA,
        SIK=TSIK, SIKMin=TSIKMn, SIKMax=TSIKMx, SIKIA=TSIKIA,
        SIDel_GN=TSIDGN, Del_GNMin=TDGNMn, Del_GNMax=TDGNMx, Del_GNIA=TDGNIA,
        SIG=TSIG, SIGMin=TSIGMn, SIGMax=TSIGMx, SIGIA=TSIGIA,
        SIEP=TSIEP, SIEPMin=TSIEPMn, SIEPMax=TSIEPMx, SIEPIA=TSIEPIA,
        PIPMin=TPIPMin, PIPMax=TPIPMax, PIPAvg=TPIPAvg,
        Fuk=TFuk, FukMin=TFukMin, FukMax=TFukMax, FukAvg=TFukAVG,
        Lapl=TLapl, LaplMin=TLaplMin, LaplMax=TLaplMax, LaplAvg=TLaplAVG,
        chi=chi,
    )

    for k in range(10):
        out['Del_RhoNA%d'%(k+1)] = TDRNA[k]
        out['Del_KNA%d'%(k+1)]   = TDKNA[k]
        out['SIKA%d'%(k+1)]      = TSIKA[k]
        out['Del_GNA%d'%(k+1)]   = TDGNA[k]
        out['SIGA%d'%(k+1)]      = TSIGA[k]
        out['SIEPA%d'%(k+1)]     = TSIEPA[k]
        out['EP%d'%(k+1)]        = TEP[k]
        out['Fuk%d'%(k+1)]       = TFukBins[k]
        out['Lapl%d'%(k+1)]      = TLaplBins[k]
        out['FDRNA%d'%(k+1)]     = FDRNA[k]
        out['FDKNA%d'%(k+1)]     = FDKNA[k]
        out['FSIKA%d'%(k+1)]     = FSIKA[k]
        out['FDGNA%d'%(k+1)]     = FDGNA[k]
        out['FSIGA%d'%(k+1)]     = FSIGA[k]
        out['FEP%d'%(k+1)]       = FEP[k]
        out['FFuk%d'%(k+1)]      = FFuk[k]
        out['FLapl%d'%(k+1)]     = FLapl[k]

    for k in range(20):
        out['PIP%d'%(k+1)]  = TPIP[k]
        out['FPIP%d'%(k+1)] = FPIP[k]

    out['atom_records'] = atom_records
    return out
