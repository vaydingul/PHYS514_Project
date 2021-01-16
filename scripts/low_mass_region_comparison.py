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

# Data path
DATA_PATH = "data/white_dwarf_data.csv"

# Data reading and partitioning
data = u.read_dw_data(DATA_PATH)
logg = data[:, 0]
mass = data[:, 1]

# Variable conversions
mass_kg, R_m, R_aer = u.general_conversion(mass, logg)

# Low mass part
low_mass_ix = mass < 0.34 


################ Low Mass ==> mass < 0.34 #################################
################ Curve Fit to find K and q ################################
# Initial guesses for K and q
q_0 = 3.0
A_0 = 9e50
A1, q1 = s.white_dwarf_fit(
    M=mass_kg[low_mass_ix], R=R_m[low_mass_ix], A=A_0, q=q_0, type="Aq")
K1 = u.Aq2K(A1, np.round(q1))
print("Estimated:\nK = {0}\nq = {1}\nA = {2}".format(K1, q1, A1))

M_Aq_fit = s.mass_radius_relation(R=R_m, K=K1, q=q1)
rho_c_Aq = s.calculate_rho_c(M=M_Aq_fit, K=K1, n=u.q2n(q1))
###########################################################################

################ Curve Fit to find K, where q is known ####################
# Initial guesses for K and q
q_0 = np.round(q1)
A_0 = A1
A2 = s.white_dwarf_fit(
    M=mass_kg[low_mass_ix], R=R_m[low_mass_ix], A=A_0, q=q_0, type="A")

K2 = u.Aq2K(A2, q_0)
print("Estimated:\nK = {0}\nq = {1}\nA = {2}".format(K2, q_0, A2))

M_A_fit = s.mass_radius_relation(R=R_m, K=K1, q=q1)
rho_c_A = s.calculate_rho_c(M=M_A_fit, K=K2, n=u.q2n(q_0))
###########################################################################
###########################################################################


# Low mass part
low_mass_ix = mass < 0.5 


################ Low Mass ==> mass < 0.5 #################################
################ Curve Fit to find K and q ################################
# Initial guesses for K and q
q_0 = 3.0
A_0 = 9e50
A1, q1 = s.white_dwarf_fit(
    M=mass_kg[low_mass_ix], R=R_m[low_mass_ix], A=A_0, q=q_0, type="Aq")
K1 = u.Aq2K(A1, np.round(q1))
print("Estimated:\nK = {0}\nq = {1}\nA = {2}".format(K1, q1, A1))

M_Aq_fit_ = s.mass_radius_relation(R=R_m, K=K1, q=q1)
rho_c_Aq_ = s.calculate_rho_c(M=M_Aq_fit_, K=K1, n=u.q2n(q1))
###########################################################################

################ Curve Fit to find K, where q is known ####################
# Initial guesses for K and q
q_0 = np.round(q1)
A_0 = A1
A2 = s.white_dwarf_fit(
    M=mass_kg[low_mass_ix], R=R_m[low_mass_ix], A=A_0, q=q_0, type="A")

K2 = u.Aq2K(A2, q_0)
print("Estimated:\nK = {0}\nq = {1}\nA = {2}".format(K2, q_0, A2))

M_A_fit_ = s.mass_radius_relation(R=R_m, K=K1, q=q1)
rho_c_A_ = s.calculate_rho_c(M=M_A_fit_, K=K2, n=u.q2n(q_0))
###########################################################################
###########################################################################

################### PLOTTING ##############################################

'''
p.figure_()
p.draw(handle="loglog", x_list=[mass], y_list=[rho_c_Aq, rho_c_A], labels=["Curve Fit for $K$ and $q$", "Curve Fit for $K$"],
       xlabel="Mass [$M_\odot$]", ylabel="Central Density [$\dfrac{kg}{m^3}$]", title="Central Density Distributions", ls=" ", savefn="report/figures/13_n_ll_rho_m_034.png")

p.figure_()
p.draw(handle="loglog", x_list=[mass[low_mass_ix]], y_list=[rho_c_Aq[low_mass_ix], rho_c_A[low_mass_ix]], labels=["Curve Fit for $K$ and $q$", "Curve Fit for $K$"],
       xlabel="Mass [$M_\odot$]", ylabel="Central Density [$\dfrac{kg}{m^3}$]", title="Central Density Distributions\n[$M \leq 0.34 M_\odot$]",
       ls=" ",  savefn="report/figures/14_n_ll_rho_m__034.png")

p.figure_()
p.draw(handle="loglog", x_list=[mass], y_list=[rho_c_Aq_, rho_c_A_], labels=["Curve Fit for $K$ and $q$", "Curve Fit for $K$"],
       xlabel="Mass [$M_\odot$]", ylabel="Central Density [$\dfrac{kg}{m^3}$]", title="Central Density Distributions", ls=" ", savefn="report/figures/15_n_ll_rho_m_050.png")

p.figure_()
p.draw(handle="loglog", x_list=[mass[low_mass_ix]], y_list=[rho_c_Aq_[low_mass_ix], rho_c_A_[low_mass_ix]], labels=["Curve Fit for $K$ and $q$", "Curve Fit for $K$"],
       xlabel="Mass [$M_\odot$]", ylabel="Central Density [$\dfrac{kg}{m^3}$]", title="Central Density Distributions\n[$M \leq 0.5 M_\odot$]",
       ls=" ",  savefn="report/figures/16_n_ll_rho_m__050.png")

p.show_()
'''

############################################################################