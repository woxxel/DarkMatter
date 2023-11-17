import os,sys

print(os.path.dirname(__file__))
sys.path.append(os.path.dirname(__file__))

from darkMatter import *

import demo_scripts as demo