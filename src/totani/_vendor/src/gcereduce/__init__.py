"""Fermi-LAT data reduction (Track B, Task 6).

DELIBERATELY OUTSIDE `gcefisher`, which enables jax x64 and imports gcepy in its
`__init__`. This package runs under the separate `fermi` conda environment and
must never import jax -- keeping it a distinct package is what enforces that.
The two halves of the project communicate through files on disk.
"""
