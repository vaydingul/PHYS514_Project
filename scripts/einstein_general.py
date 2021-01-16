# Ready-to-use module import

from scipy.interpolate import interp1d
import numpy as np
import sys
sys.path.insert(1, "./")
sys.path.insert(2, "./../")

# Custom module import
from src import Utils as u
from src import Solvers as s
from src import Plotter as p
from src import Constants as c


p_c_sample = np.linspace(1e31)