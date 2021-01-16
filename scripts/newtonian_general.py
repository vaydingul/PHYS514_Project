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
low_mass_ix = mass < 0.34  # mean_mass

# Variable conversions
mass_kg, R_m, R_aer = u.general_conversion(mass, logg)

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

M_A_fit = s.mass_radius_relation(R=R_m, K=K2, q=q_0)
rho_c_A = s.calculate_rho_c(M=M_A_fit, K=K2, n=u.q2n(q_0))
###########################################################################


'''
################ Curve Fit to find K and q ################################
# Initial guesses for K and q
q_0 = 3.0
K_0 = 4e6
K1, q1 = s.white_dwarf_fit(
    M=mass_kg[low_mass_ix], R=R_m[low_mass_ix], K=K_0, q=q_0, type="Kq")
print("Estimated:\nK = {0}\nq = {1}".format(K1, q1))

M_Kq_fit = s.mass_radius_relation(R=R_m, K=K1, q=q1)
rho_c_Kq = s.calculate_rho_c(M = M_Kq_fit, K = K1, n = u.q2n(q1))
###########################################################################

################ Curve Fit to find K, where q is known ####################
# Initial guesses for K and q
q_0 = np.round(q1)
K_0 = K1
K2 = s.white_dwarf_fit(
    M=mass_kg[low_mass_ix], R=R_m[low_mass_ix], K=K_0, q=q_0, type="K")
print("Estimated:\nK = {0}\nq = {1}".format(K2, q_0))

M_K_fit = s.mass_radius_relation(R=R_m, K=K2, q=q_0)
rho_c_K = s.calculate_rho_c(M = M_K_fit, K = K1, n = u.q2n(q_0))
###########################################################################
'''


################ Curve Fit to find D where K and q are known ##############
rho_c_0 = 1e9
# Sampling points for rho_c
#rho_c_list = np.linspace(rho_c_0 / 10, rho_c_0 * 10, 10)
rho_c_list = np.linspace(np.min(rho_c_A), np.max(rho_c_A), 20)
# Initial guess for D
K_0 = K2
q_0 = np.round(q1)
# ! TODO: Play with the initial guess of the D_0
D_0 = 2e9
D1 = s.white_dwarf_fit(M=mass_kg[low_mass_ix], R=R_m[low_mass_ix],
                       K=K_0, q=q_0, rho_c_list=rho_c_list, D=D_0, type="D")

print("Estimated:\nK = {0}\nq = {1}\nD = {2}\nC = {3}".format(K_0, q_0, D1, u.KqD2C(K_0, q_0, D1)))

MR = [s.mass_radius_relation(K=K_0, q=q_0, D=D1, rho_c=rho_c)
      for rho_c in rho_c_A]
MR = np.array([*MR])
M_D_fit = MR[:, 0]
R_D_fit = MR[:, 1]
###########################################################################

########### THEORETICAL ASSESSMENT ########################################
q_theo = q_0
D_theo = c.D
K_theo = u.CqD2K(c.C, q_theo, D_theo)
print("Theoritical:\nK = {0}\nq = {1}\nD = {2}\nC = {3}".format(
    K_theo, q_theo, D_theo, c.C))

MR_theo = [s.mass_radius_relation(K=K_theo, q=q_theo, D=D_theo, rho_c=rho_c)
           for rho_c in rho_c_A]
MR_theo = np.array([*MR_theo])
M_theo = MR_theo[:, 0]
R_theo = MR_theo[:, 1]
###########################################################################
################# PLOTTING ################################################ 
p.figure_()
p.draw(handle="loglog", x_list=[mass], y_list=[rho_c_Aq, rho_c_A], labels=["Curve Fit for $K$ and $q$", "Curve Fit for $K$"],
       xlabel="Mass [$M_\odot$]", ylabel="Central Density [$\dfrac{kg}{m^3}$]", title="Central Density Distributions", ls=" ", savefn="report/figures/1_n_ll_rho_m.png")

p.figure_()
p.draw(handle="loglog", x_list=[mass[low_mass_ix]], y_list=[rho_c_Aq[low_mass_ix], rho_c_A[low_mass_ix]], labels=["Curve Fit for $K$ and $q$", "Curve Fit for $K$"],
       xlabel="Mass [$M_\odot$]", ylabel="Central Density [$\dfrac{kg}{m^3}$]", title="Central Density Distributions\n[Only Low Mass Stars]",
       ls=" ",  savefn="report/figures/2_n_ll_rho_m_.png")

p.figure_()
p.draw(handle="scatter", x_list=[mass], y_list=[rho_c_Aq, rho_c_A], labels=["Curve Fit for $K$ and $q$", "Curve Fit for $K$"],
       xlabel="Mass [$M_\odot$]", ylabel="Central Density [$\dfrac{kg}{m^3}$]", title="Central Density Distributions", savefn = "report/figures/3_n_s_rho_m.png")

p.figure_()
p.draw(handle="scatter", x_list=[mass[low_mass_ix]], y_list=[rho_c_Aq[low_mass_ix], rho_c_A[low_mass_ix]], labels=["Curve Fit for $K$ and $q$", "Curve Fit for $K$"],
       xlabel="Mass [$M_\odot$]", ylabel="Central Density [$\dfrac{kg}{m^3}$]", title="Central Density Distributions\n[Only Low Mass Stars]", savefn = "report/figures/4_n_s_rho_m_.png")

p.figure_()
p.draw(handle="loglog", x_list=[R_aer], y_list=[mass_kg, M_Aq_fit, M_A_fit, M_D_fit, M_theo], labels=["Actual Data", "Curve Fit for $K$ and $q$", "Curve Fit for only $K$", "Curve Fit for only $D$", "Theoretical Result"],
       xlabel="Radius [Average Earth Radius]", ylabel="Mass [kg]", title="Mass Distributions", ls=" ", savefn = "report/figures/5_n_ll_m_r.png")

p.figure_()
p.draw(handle="loglog", x_list=[R_aer[low_mass_ix]], y_list=[mass_kg[low_mass_ix], M_Aq_fit[low_mass_ix], M_A_fit[low_mass_ix], M_D_fit[low_mass_ix], M_theo[low_mass_ix]], labels=["Actual Data", "Curve Fit for $K$ and $q$", "Curve Fit for only $K$", "Curve Fit for only $D$", "Theoretical Result"],
       xlabel="Radius [Average Earth Radius]", ylabel="Mass [kg]", title="Mass Distributions\n[Only Low Mass Stars]", ls=" ", savefn = "report/figures/6_n_ll_m_r_.png")


p.figure_()
p.draw(handle="loglog", x_list=[R_aer], y_list=[mass, M_Aq_fit/c.SOLAR_MASS, M_A_fit/c.SOLAR_MASS, M_D_fit/c.SOLAR_MASS, M_theo/c.SOLAR_MASS], labels=["Actual Data", "Curve Fit for $K$ and $q$", "Curve Fit for only $K$", "Curve Fit for only $D$", "Theoretical Result"],
       xlabel="Radius [Average Earth Radius]", ylabel="Mass [$M_\odot$]", title="Mass Distributions", ls=" ", savefn = "report/figures/7_n_ll_ms_r.png")

p.figure_()
p.draw(handle="loglog", x_list=[R_aer[low_mass_ix]], y_list=[mass[low_mass_ix], M_Aq_fit[low_mass_ix]/c.SOLAR_MASS, M_A_fit[low_mass_ix]/c.SOLAR_MASS, M_D_fit[low_mass_ix]/c.SOLAR_MASS, M_theo[low_mass_ix]/c.SOLAR_MASS], labels=["Actual Data", "Curve Fit for $K$ and $q$", "Curve Fit for only $K$", "Curve Fit for only $D$", "Theoretical Result"],
       xlabel="Radius [Average Earth Radius]", ylabel="Mass [$M_\odot$]", title="Mass Distributions\n[Only Low Mass Stars]", ls=" ", savefn = "report/figures/8_n_ll_ms_r_.png")

p.figure_()
p.draw(handle="scatter", x_list=[R_aer], y_list=[mass_kg, M_Aq_fit, M_A_fit, M_D_fit, M_theo], labels=["Actual Data", "Curve Fit for $K$ and $q$", "Curve Fit for only $K$", "Curve Fit for only $D$", "Theoretical Result"],
       xlabel="Radius [Average Earth Radius]", ylabel="Mass [kg]", title="Mass Distributions", savefn = "report/figures/9_n_s_m_r.png")

p.figure_()
p.draw(handle="scatter", x_list=[R_aer[low_mass_ix]], y_list=[mass_kg[low_mass_ix], M_Aq_fit[low_mass_ix], M_A_fit[low_mass_ix], M_D_fit[low_mass_ix], M_theo[low_mass_ix]], labels=["Actual Data", "Curve Fit for $K$ and $q$", "Curve Fit for only $K$", "Curve Fit for only $D$", "Theoretical Result"],
       xlabel="Radius [Average Earth Radius]", ylabel="Mass [kg]", title="Mass Distributions\n[Only Low Mass Stars]", savefn = "report/figures/10_n_s_m_r_.png")


p.figure_()
p.draw(handle="scatter", x_list=[R_aer], y_list=[mass, M_Aq_fit/c.SOLAR_MASS, M_A_fit/c.SOLAR_MASS, M_D_fit/c.SOLAR_MASS, M_theo/c.SOLAR_MASS], labels=["Actual Data", "Curve Fit for $K$ and $q$", "Curve Fit for only $K$", "Curve Fit for only $D$", "Theoretical Result"],
       xlabel="Radius [Average Earth Radius]", ylabel="Mass [$M_\odot$]", title="Mass Distributions", savefn = "report/figures/11_n_s_ms_r.png")

p.figure_()
p.draw(handle="scatter", x_list=[R_aer[low_mass_ix]], y_list=[mass[low_mass_ix], M_Aq_fit[low_mass_ix]/c.SOLAR_MASS, M_A_fit[low_mass_ix]/c.SOLAR_MASS, M_D_fit[low_mass_ix]/c.SOLAR_MASS, M_theo[low_mass_ix]/c.SOLAR_MASS], labels=["Actual Data", "Curve Fit for $K$ and $q$", "Curve Fit for only $K$", "Curve Fit for only $D$", "Theoretical Result"],
       xlabel="Radius [Average Earth Radius]", ylabel="Mass [$M_\odot$]", title="Mass Distributions\n[Only Low Mass Stars]", savefn = "report/figures/12_n_s_ms_r_.png")


p.show_()
