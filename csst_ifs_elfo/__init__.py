"""
Identifier:     csst_ifs_elfo/__init__.py
Name:           __init__.py
Description:    IFS elfo package
Author:         Hui Guo
Created:        2025-02-09
Modified-History:
    2025-02-09, Hui Guo, created

"""

import os
from .para import create_qsopar
from .elfo import process_i_refit
from .elfo import process_j_refit
from .elfo import downward
from .elfo import upward
from .elfo import leftward
from .elfo import rightward
from .elfo import process_down
from .elfo import process_spectrum_down
from .elfo import process_up
from .elfo import process_spectrum_up
from .elfo import process_left
from .elfo import process_spectrum_left
from .elfo import process_right
from .elfo import process_spectrum_right
from .elfo import process_i
from .elfo import process_j
from .elfo import process_spectrum
from .prepare_data import combine_pixels
__version__ = "0.1.0"

PACKAGE_PATH = os.path.dirname(__file__)
__all__ = [
    'create_qsopar',
    'process_i_refit',
    'combine_pixels',
    'process_j_refit',
    'downward',
    'upward',
    'leftward',
    'rightward',
    'process_down',
    'process_spectrum_down',
    'process_up',
    'process_spectrum_up',
    'process_left',
    'process_spectrum_left',
    'process_right',
    'process_spectrum_right',
    'process_i',
    'process_j',
    'process_spectrum',
]
