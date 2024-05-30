"""Seawater: EOS-80 equation of state."""

import warnings

try:
    from ._version import __version__
except ImportError:
    __version__ = "unknown"

from .eos80 import (
    adtg,
    alpha,
    aonb,
    beta,
    cp,
    dens,
    dens0,
    dpth,
    fp,
    g,
    pden,
    pres,
    ptmp,
    salt,
    svel,
    temp,
)
from .extras import dist, f, satAr, satN2, satO2, swvel
from .geostrophic import bfrq, gpan, gvel, svan
from .library import cndr, salds, salrp, salrt, sals, seck, smow

__all__ = [
    "adtg",
    "alpha",
    "aonb",
    "beta",
    "bfrq",
    "cndr",
    "cp",
    "dens",
    "dens0",
    "dist",
    "dpth",
    "f",
    "fp",
    "g",
    "gpan",
    "gvel",
    "pden",
    "pres",
    "ptmp",
    "salds",
    "salrp",
    "salrt",
    "sals",
    "salt",
    "satAr",
    "satN2",
    "satO2",
    "seck",
    "smow",
    "svan",
    "svel",
    "swvel",
    "temp",
]

warnings.warn(
    "The seawater library is deprecated! Please use gsw instead.",
    stacklevel=2,
)
