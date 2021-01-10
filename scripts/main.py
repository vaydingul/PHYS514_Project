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
mean_mass = np.mean(mass)

# Low mass part
low_mass_ix = mass < mean_mass

# Variable conversions
mass_kg, R_m, R_aer = u.general_conversion(mass, logg)

################ Curve Fit to find K and q ################################
# Initial guesses for K and q
q_0 = 3.0
K_0 = 5e13
K1, q1 = s.white_dwarf_fit(
    M=mass_kg[low_mass_ix], R=R_m[low_mass_ix], K=K_0, q=q_0, type="Kq")
print("Estimated:\nK = {0}\nq = {1}".format(K1, q1))

M_Kq_fit = s.mass_radius_relation(R = R_m,K = K1, q = q1)
rho_c_Kq = s.calculate_rho_c(M_Kq_fit, K1, u.q2n(q1))
###########################################################################

################ Curve Fit to find K, where q is known ####################
# Initial guesses for K and q
q_0 = np.round(q1)
K_0 = K1
K2 = s.white_dwarf_fit(
    M=mass_kg[low_mass_ix], R=R_m[low_mass_ix], K=K_0, q=q_0, type="K")
print("Estimated:\nK = {0}\nq = {1}".format(K2, q_0))

M_K_fit = s.mass_radius_relation(R = R_m,K = K2, q = q_0)
rho_c_K = s.calculate_rho_c(M_K_fit, K1, u.q2n(q_0))
###########################################################################
# ! TODO: Refactor calculate_rho_c into a scheme where it is calculated via R, rather than M

################ Curve Fit to find D where K and q are known ##############
rho_c_0 = 1e9
# Sampling points for rho_c
rho_c_list = np.linspace(rho_c_0 / 10, rho_c_0 * 10, 10)
# Initial guess for D
K_0 = K2
q_0 = np.round(q1)
# ! TODO: Play with the initial guess of the D_0
D_0 = 1e7
D1 = s.white_dwarf_fit(M=mass_kg[low_mass_ix], R=R_m[low_mass_ix],
                      K=K_0, q=q_0, rho_c_list=rho_c_list, D=D_0, type="D")

print("Estimated:\nK = {0}\nq = {1}\nD = {2}".format(K_0, q_0, D1))

MR = [s.mass_radius_relation(K = K_0, q = q_0, D = D1, rho_c = rho_c) for rho_c in rho_c_K]
MR = np.array([*MR])
M_D_fit = MR[:, 0]
R_D_fit = MR[:, 1]
###########################################################################


################# PLOTTING ################################################
p.figure_()
p.scatter_([mass], [rho_c_Kq, rho_c_K], labels = ["Curve Fit for $K$ and $q$", "Curve Fit for $K$"],
        xlabel = "Mass [Solar Mass]", ylabel = "Central Density [$\dfrac{kg}{m^3}$]", title = "Central Density Distributions")

p.figure_()
p.scatter_([mass[low_mass_ix]], [rho_c_Kq[low_mass_ix], rho_c_K[low_mass_ix]], labels = ["Curve Fit for $K$ and $q$", "Curve Fit for $K$"],
        xlabel = "Mass [Solar Mass]", ylabel = "Central Density [$\dfrac{kg}{m^3}$]", title = "Central Density Distributions\n[Only Low Mass Starts]")


p.figure_()
p.scatter_([R_aer], [mass_kg, M_Kq_fit, M_K_fit, M_D_fit], labels = ["Actual Data", "Curve Fit for $K$ and $q$", "Curve Fit for only $K$", "Curve Fit for only $D$"],
        xlabel = "Radius [Average Earth Radius]", ylabel = "Mass [kg]", title = "Mass Distributions")

p.figure_()
p.scatter_([R_aer[low_mass_ix]], [mass_kg[low_mass_ix], M_Kq_fit[low_mass_ix], M_K_fit[low_mass_ix], M_D_fit[low_mass_ix]], labels = ["Actual Data", "Curve Fit for $K$ and $q$", "Curve Fit for only $K$", "Curve Fit for only $D$"],
        xlabel = "Radius [Average Earth Radius]", ylabel = "Mass [kg]", title = "Mass Distributions\n[Only Low Mass Starts]")


p.show_()

'''
rho_c_ = np.linspace(rho_c/10, rho_c * 10, 5)
MR = [s.mass_radius_relation_D(rho_c, D_0) for rho_c in np.random.choice(rho_K, 20)]
MR = np.array([*MR])
M_ = MR[:, 0]/c.SOLAR_MASS; R_ = MR[:, 1]
sp_f = interp1d(R_, M_, kind = "cubic", fill_value = "extrapolate")
M_spline = sp_f(R_m[low_mass_ix])
p.figure_()
p.scatter_([R_m[low_mass_ix]], [mass[low_mass_ix], M_spline], labels = ["orgi", "splinişko"])
p.show_()
'''
#res = s.white_dwarf_fit_D(mass_kg, R_m, 1e8, 1e5)
# print(res)
# p.figure_()
#p.scatter_([R_m[low_mass_ix]], [rho_Kq[low_mass_ix], rho_K[low_mass_ix]])
#p.scatter_([R_m, R_], [mass_kg/c.SOLAR_MASS, M_/c.SOLAR_MASS])
# p.show_()

'''
p.figure_()
p.scatter_([R_m], [mass, M_Kq, M_K], labels = ["data","fit_Kq", "fit_K"],
        xlabel = "Radius [m]", ylabel = "Mass [$M_\odot$]")
p.figure_()
p.scatter_([R_m[low_mass_ix]], [mass[low_mass_ix], M_Kq[low_mass_ix], M_K[low_mass_ix]],
        labels = ["data","fit_Kq", "fit_K"], xlabel = "Radius [m]", ylabel = "Mass [$M_\odot$]", title = "Zoomed Version")

p.show_()
'''
