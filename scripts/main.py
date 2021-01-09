import sys
sys.path.insert(1, "./")
sys.path.insert(2, "./../")

# Custom module import
from src import Utils as u
from src import Solvers as s
from src import Plotter as p
from src import Constants as c
# Ready-to-use module import
import numpy as np

DATA_PATH = "data/white_dwarf_data.csv"

# Data reading and partitioning
data = u.read_dw_data(DATA_PATH)
logg = data[:,0]; mass = data[:,1]; mean_mass = np.mean(mass)

low_mass_ix = mass < mean_mass
mass_low = mass[low_mass_ix]
# Variable conversions
mass_kg, R_m, R_aer = u.logg_2_aer(mass, logg)


#t,y = s.solve_lane_emden(1.5)

mass_kg_low = mass_kg[low_mass_ix]
R_m_low = R_m[low_mass_ix]

q_0 = 3.0
K_0 = 5e13
K, q = s.white_dwarf_fit_Kq(mass_kg_low, R_m_low, K_0, q_0)
K = s.white_dwarf_fit_Kq(mass_kg_low, R_m_low, K_0)

fit_Kq = s.mass_radius_relation_Kq(R_m, K, q) / c.SOLAR_MASS
fit_K = s.mass_radius_relation_K(R_m, K) / c.SOLAR_MASS

p.figure_()
p.scatter_([R_m], [mass, fit_Kq, fit_K], labels = ["data","fit_Kq", "fit_K"],
        xlabel = "Radius [m]", ylabel = "Mass [$M_\odot\$]")
p.figure_()
p.scatter_([R_m[low_mass_ix]], [mass[low_mass_ix], fit_Kq[low_mass_ix], fit_K[low_mass_ix]],
        labels = ["data","fit_Kq", "fit_K"], xlabel = "Radius [m]", ylabel = "Mass [$M_\odot$]", title = "Zoomed Version")

p.show_()