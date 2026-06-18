"""Lightweight sparse 2D integer matrix used for bond-order tables.

Behaves like a dense ``[[0]*N for _ in range(N)]`` matrix for the
``nbo[i][j]`` read/write pattern used throughout the readers, but only
stores entries that have actually been set, avoiding the O(MAX^2)
memory cost of a fully dense matrix sized for the Fortran
``parameter(max=1000)`` worst case.
"""
from collections import defaultdict


class SparseRow(object):
    __slots__ = ("_data", "_row")

    def __init__(self, data, row):
        self._data = data
        self._row = row

    def __getitem__(self, col):
        return self._data.get((self._row, col), 0)

    def __setitem__(self, col, value):
        if value:
            self._data[(self._row, col)] = value
        else:
            self._data.pop((self._row, col), None)


class SparseMatrix(object):
    """Drop-in replacement for a dense NxN int matrix, used as ``nbo``."""
    __slots__ = ("_data",)

    def __init__(self):
        self._data = {}

    def __getitem__(self, row):
        return SparseRow(self._data, row)
