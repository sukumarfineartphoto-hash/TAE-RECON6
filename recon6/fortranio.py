"""
Minimal reader for gfortran unformatted sequential files.

Each Fortran ``write(unit) ...`` statement produces one "record" on disk,
framed by a 4-byte little-endian record-length marker before and after
the record payload (the classic gfortran/g77 convention).
"""
import struct


class FortranSequentialReader:
    """Read records from a gfortran unformatted sequential file."""

    def __init__(self, path):
        self._f = open(path, "rb")

    def close(self):
        self._f.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def read_record(self):
        """Return the raw bytes of the next record, or None at EOF."""
        head = self._f.read(4)
        if not head:
            return None
        if len(head) != 4:
            raise IOError("Truncated record header")
        (length,) = struct.unpack("<i", head)
        data = self._f.read(length)
        if len(data) != length:
            raise IOError("Truncated record body")
        tail = self._f.read(4)
        if len(tail) != 4:
            raise IOError("Truncated record trailer")
        (length2,) = struct.unpack("<i", tail)
        if length2 != length:
            raise IOError("Record length mismatch")
        return data

    def skip_record(self):
        data = self.read_record()
        if data is None:
            raise EOFError("Unexpected end of file")
        return data

    def read_doubles(self, n):
        data = self.skip_record()
        if len(data) != 8 * n:
            raise IOError(
                "Expected %d doubles (%d bytes), got %d bytes" % (n, 8 * n, len(data))
            )
        return struct.unpack("<%dd" % n, data)

    def read_ints(self, n):
        data = self.skip_record()
        if len(data) != 4 * n:
            raise IOError(
                "Expected %d ints (%d bytes), got %d bytes" % (n, 4 * n, len(data))
            )
        return struct.unpack("<%di" % n, data)

    def read_int(self):
        return self.read_ints(1)[0]

    def read_char(self, n):
        data = self.skip_record()
        if len(data) != n:
            raise IOError("Expected %d-byte char record, got %d bytes" % (n, len(data)))
        return data.decode("ascii", errors="replace")
