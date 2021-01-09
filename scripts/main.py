import sys
sys.path.insert(1, "./")
sys.path.insert(2, "./../")

# Custom module import
from src import Utils as u
from src import Solvers as s
from src import Plotter as p
# Ready-to-use module import
import numpy as np
import matplotlib.pyplot as plt

DATA_PATH = "data/white_dwarf_data.csv"

# Data reading and partitioning
data = u.read_dw_data(DATA_PATH)
logg = data[:,0]; mass = data[:,1]
# Variable conversions
mass_kg, R_m, R_aer = u.logg_2_aer(mass, logg)

p.scatter_([R_aer], [mass])
p.show_()
