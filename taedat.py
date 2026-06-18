"""
Reader for TAE atom-type ``.dat`` files (binary, gfortran unformatted
sequential format).  The record layout mirrors the ``write(4) ...``
statements in the original ``binarydat4.f`` converter program.
"""
from .fortranio import FortranSequentialReader


class TaeRecord(object):
    """Holds all fields read from one TAE atom-type .dat file."""
    __slots__ = ("fields",)

    def __init__(self, fields):
        self.fields = fields

    def __getattr__(self, name):
        try:
            return self.fields[name]
        except KeyError:
            raise AttributeError(name)


def _read_group(r, names):
    vals = r.read_doubles(len(names))
    return dict(zip(names, vals))


def read_tae_dat(path):
    """
    Read a TAE .dat file and return a TaeRecord with all named fields.

    This mirrors the read sequence in binarydat4.f.
    """
    fields = {}
    with FortranSequentialReader(path) as r:
        # First two records are read with bare `read(11)` in qmf.f (no
        # variables), so their actual sizes/contents don't matter and are
        # simply skipped.
        fields["_rec1"] = r.skip_record()
        fields["_rec2"] = r.skip_record()

        numval = r.read_int()
        fields["numval"] = numval

        mvaln = None
        for _ in range(numval):
            # itype (int4), blen (real*8)
            r.skip_record()
            # numvaln (int4), ii(1:numvaln-1) (int4 array)
            data2 = r.skip_record()
            import struct
            numvaln = struct.unpack("<i", data2[:4])[0]
            mvaln = numvaln
            # rad,theta1,phi1,rhoval
            r.read_doubles(4)
            # xeigen,yeigen,zeigen,eival
            r.read_doubles(4)

        fields["mvaln"] = mvaln

        if numval == 1 and mvaln is not None:
            for _ in range(mvaln - 1):
                r.skip_record()  # iat, ivalnn, xvec1, yvec1, zvec1

        fields.update(_read_group(r, ["atdx", "atdy", "atdz"]))
        fields.update(_read_group(r, ["xatom", "yatom", "zatom", "energy", "pop"]))
        fields.update(_read_group(r, ["vol", "dipm", "dxm", "dym", "dzm", "qaa"]))
        fields.update(_read_group(r, ["qbb", "qcc", "qxxm", "qxym", "qxzm", "qyym"]))
        fields.update(_read_group(r, ["qyzm", "qzzm", "soxxx", "soyyy", "sozzz"]))
        fields.update(_read_group(r, ["soxxy", "soxxz", "soyyx", "soyyz", "sozzx", "sozzy"]))
        fields.update(_read_group(r, ["soxyz", "shxxxx", "shyyyy", "shzzzz", "shxxxy", "shxxxz"]))
        fields.update(_read_group(r, ["shyyyx", "shyyyz", "shzzzx", "shzzzy", "shxxyy", "shxxzz"]))
        fields.update(_read_group(r, ["shyyzz", "shxxyz", "shyyxz", "shzzxy"]))

        fields.update(_read_group(r, ["sarea", "sidrn", "drnmn", "drnmx", "drnia"]))
        fields.update(_read_group(r, ["drna1", "drna2", "drna3", "drna4", "drna5"]))
        fields.update(_read_group(r, ["drna6", "drna7", "drna8", "drna9", "drna10"]))

        fields.update(_read_group(r, ["sidkn", "dkmn", "dkmx", "dkia"]))
        fields.update(_read_group(r, ["dkna1", "dkna2", "dkna3", "dkna4", "dkna5"]))
        fields.update(_read_group(r, ["dkna6", "dkna7", "dkna8", "dkna9", "dkna10"]))

        fields.update(_read_group(r, ["sik", "sikmn", "sikmx", "sikia"]))
        fields.update(_read_group(r, ["sika1", "sika2", "sika3", "sika4", "sika5"]))
        fields.update(_read_group(r, ["sika6", "sika7", "sika8", "sika9", "sika10"]))

        fields.update(_read_group(r, ["sidgn", "dgnmn", "dgnmx", "dgnia"]))
        fields.update(_read_group(r, ["dgna1", "dgna2", "dgna3", "dgna4", "dgna5"]))
        fields.update(_read_group(r, ["dgna6", "dgna7", "dgna8", "dgna9", "dgna10"]))

        fields.update(_read_group(r, ["sig", "sigmn", "sigmx", "sigia"]))
        fields.update(_read_group(r, ["siga1", "siga2", "siga3", "siga4", "siga5"]))
        fields.update(_read_group(r, ["siga6", "siga7", "siga8", "siga9", "siga10"]))

        fields.update(_read_group(r, ["siep", "siepmn", "siepmx", "siepia"]))
        fields.update(_read_group(r, ["siepa1", "siepa2", "siepa3", "siepa4", "siepa5"]))
        fields.update(_read_group(r, ["siepa6", "siepa7", "siepa8", "siepa9", "siepa10"]))

        fields.update(_read_group(r, ["piv", "sigmapv", "sigmanv", "sumsigma", "sigmanew"]))

        fields.update(_read_group(r, ["ep1", "ep2", "ep3", "ep4", "ep5"]))
        fields.update(_read_group(r, ["ep6", "ep7", "ep8", "ep9", "ep10"]))

        fields.update(_read_group(r, ["pipmin", "pipmax", "pipavg"]))
        fields.update(_read_group(r, ["p1", "p2", "p3", "p4"]))
        fields.update(_read_group(r, ["p5", "p6", "p7", "p8"]))
        fields.update(_read_group(r, ["p9", "p10", "p11", "p12"]))
        fields.update(_read_group(r, ["p13", "p14", "p15", "p16"]))
        fields.update(_read_group(r, ["p17", "p18", "p19", "p20"]))

        # NOTE: field names follow qmf.f's read-side variable names
        # (positional), which differ in min/max naming from the writer
        # (binarydat4.f) - qmf is the consumer we must match.
        fields.update(_read_group(r, ["lapl", "laplmin", "laplmax", "laplavg"]))
        fields.update(_read_group(r, ["lapl1", "lapl2", "lapl3", "lapl4"]))
        fields.update(_read_group(r, ["lapl5", "lapl6", "lapl7", "lapl8"]))
        fields.update(_read_group(r, ["lapl9", "lapl10"]))

        fields.update(_read_group(r, ["fuk", "fukmin", "fukmax", "fukavg"]))
        fields.update(_read_group(r, ["fuk1", "fuk2", "fuk3", "fuk4"]))
        fields.update(_read_group(r, ["fuk5", "fuk6", "fuk7", "fuk8"]))
        fields.update(_read_group(r, ["fuk9", "fuk10"]))

    return TaeRecord(fields)
