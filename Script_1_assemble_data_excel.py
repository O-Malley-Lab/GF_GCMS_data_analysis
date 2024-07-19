"""
GF GCMS Data Analysis Script 1
Lazarina Butkovich 7/19/24

This script compiles outputs from different GC-MS data analysis tools to create a combined data excel to reference.

Tool outputs:
MS-DIAL
GNPS

"""

import pandas as pd
import numpy as np
from os.path import join as pjoin

"""""""""""""""""""""""""""""""""""""""""""""
Values
"""""""""""""""""""""""""""""""""""""""""""""
INPUT_FOLDER = r'input' 
TEMP_OVERALL_FOLDER = r'temp'