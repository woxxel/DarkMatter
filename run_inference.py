import os, sys

root_dir = os.path.dirname(os.path.abspath(''))
if not root_dir in sys.path: sys.path.append(root_dir)

from empirical import *
from inference import *

mP = ModelParams('empirical',filePath='../data/BuscheLab/2P_data.xlsx',population_keys=['*mouse_type','animal','val1'])
mP.regularize_rates()
