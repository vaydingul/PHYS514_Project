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

def demo(Kq_optimization_alternative = 1, CD_optimization_alternative = 1, show_fig = True, print_verbose = True, mass_constraint = 0.34):

	"""
	demo:
	demo(Kq_optimization_alternative = 1, CD_optimization_alternative = 1)

	This function perform the following processes:
		- It visualizes the results obtained in Newtonian approach
		- It provides optimization techniques based on different alternatives

	Input:
		Kq_optimization_alternative = Optimization alternative to find K and q
		CD_optimization_alternative = Optimization alternative to find K and q
		show_fig = It determines whether graphs will show up or not
		mass_constraint = It is the low mass region constraint, purely for testing purposes
		
	Output:
			[]

	Example:
			[]
	"""
	# Data path
	DATA_PATH = "../data/white_dwarf_data.csv"

	# Data reading and partitioning
	data = u.read_dw_data(DATA_PATH)
	logg = data[:, 0]
	mass = data[:, 1]

	# Low mass part
	low_mass_ix = mass < mass_constraint
	
	# WD data unit conversions
	mass_kg, R_m, R_aer = u.wd_data_conversion(mass, logg)
	# Here, the variables represent that:
	# mass_kg -> Mass in kilograms
	# R_m     -> Radius in meters
	# R_aer   -> Radius in ´average Earth radius´

	
	######################### CURVE FIT FOR K and q ###############################
	if Kq_optimization_alternative == 1:
		################ Curve Fit to find K and q  concurrently ##################

		# In the first alternative, the proportional form of the 
		# mass-radius relation equation equation is used.
		# Therefore, the proportional constant in called A.
		# Purpose is to find A, and then, convert it to the K.

		# Initial guesses for A and q
		q_0 = 3.0 # Initial guess for q (Educated guess :p)
		A_0 = 9e50 # Initial guess for A (proportionality constant)
		
		# Curve fitting for A and q via using the general ´white_dwarf_fit´ routine
		A1, q1 = s.white_dwarf_fit(
		M=mass_kg[low_mass_ix], R=R_m[low_mass_ix], A=A_0, q=q_0, type="Kq", mode = Kq_optimization_alternative)
		
		# Conversion from A to K
		K1 = u.Aq2K(A1, np.round(q1))
		
		if print_verbose:
			# Verbose printing
			print("============================================================================")
			print("Curve fitting for K and q concurrently, based on Alternative 1")
			print("Estimated:\nK = {0}\nq = {1}\nA = {2}".format(K1, q1, A1))

		# Calculation of mass based on the estimated K and q value
		M_Kq_fit = s.mass_radius_relation(R=R_m, K=K1, q=q1)
		# Calculation of central density based on the estimated K and q value
		rho_c_Kq = s.calculate_rho_c(M=M_Kq_fit, K=K1, n=u.q2n(q1))
		###########################################################################

		################ Curve Fit to find K, where q is fixed ####################

		# Here, the idea is to make further estimation via fixing the one of the
		# variables that is found in previous iteration.

		# Initial guesses for A
		q_0 = np.round(q1) # Fixation of q
		A_0 = A1 # Initial guess for A

		# Curve fitting for A and q via using the general ´white_dwarf_fit´ routine
		A2 = s.white_dwarf_fit(
		M=mass_kg[low_mass_ix], R=R_m[low_mass_ix], A=A_0, q=q_0, type="K")

		# Conversion from A to K
		K2 = u.Aq2K(A2, q_0)

		if print_verbose:
			# Verbose printing
			print("============================================================================")
			print("Curve fitting for K while q is fixed, based on Alternative 1")
			print("Estimated:\nK = {0}\nq = {1}\nA = {2}".format(K2, q_0, A2))

		# Calculation of mass based on the estimated K and q value
		M_K_fit = s.mass_radius_relation(R=R_m, K=K2, q=q_0)
		# Calculation of central density based on the estimated K and q value
		rho_c_K = s.calculate_rho_c(M=M_K_fit, K=K2, n=u.q2n(q_0))
		###########################################################################

	else:
		
		################ Curve Fit to find K and q concurrently ##################

		# In the second alternative, the full form of the 
		# mass-radius relation equation equation is used.
		# Therefore, the fitting process is directly 
		# applied on the K and q

		# Initial guesses for K and q
		q_0 = 3.0 # Initial guess for q
		K_0 = 4e6 # Initial guess for K

		# Curve fitting for K and q via using the general ´white_dwarf_fit´ routine
		K1, q1 = s.white_dwarf_fit(
		M=mass_kg[low_mass_ix], R=R_m[low_mass_ix], K=K_0, q=q_0, type="Kq", mode = Kq_optimization_alternative)

		if print_verbose:
			# Verbose printing
			print("============================================================================")
			print("Curve fitting for K and q concurrently, based on Alternative 2")
			print("Estimated:\nK = {0}\nq = {1}".format(K1, q1))

		# Calculation of mass based on the estimated K and q value
		M_Kq_fit = s.mass_radius_relation(R=R_m, K=K1, q=q1)
		# Calculation of central density based on the estimated K and q value
		rho_c_Kq = s.calculate_rho_c(M = M_Kq_fit, K = K1, n = u.q2n(q1))
		###########################################################################

		################ Curve Fit to find K, where q is known ####################

		# Here, the idea is to make further estimation via fixing the one of the
		# variables that is found in previous iteration.

		# Initial guesses for K and q
		q_0 = np.round(q1) # Fixation of q
		K_0 = K1 # Initial guess for K

		# Curve fitting for K and q via using the general ´white_dwarf_fit´ routine
		K2 = s.white_dwarf_fit(
		M=mass_kg[low_mass_ix], R=R_m[low_mass_ix], K=K_0, q=q_0, type="K", mode = Kq_optimization_alternative)

		if print_verbose:
			# Verbose printing
			print("============================================================================")
			print("Curve fitting for K while q is fixed, based on Alternative 2")
			print("Estimated:\nK = {0}\nq = {1}".format(K2, q_0))

		# Calculation of mass based on the estimated K and q value
		M_K_fit = s.mass_radius_relation(R=R_m, K=K2, q=q_0)
		# Calculation of central density based on the estimated K and q value
		rho_c_K = s.calculate_rho_c(M = M_K_fit, K = K2, n = u.q2n(q_0))
		###########################################################################
		


	################ Curve Fit to find D while K and q are known ##############
	# Sampling points for central density, based on previously calculated emprical
	# central densities
	rho_c_list = np.linspace(np.min(rho_c_K), np.max(rho_c_K), 20)

	K_0 = K2 # Fixation of K
	q_0 = np.round(q1) # Fixation of q
	D_0 = 2e9 # Initial guess for D

	# Curve fitting for D via using the general ´white_dwarf_fit´ routine
	D1 = s.white_dwarf_fit(M=mass_kg[low_mass_ix], R=R_m[low_mass_ix],
						K=K_0, q=q_0, rho_c_list=rho_c_list, D=D_0, type="D", mode = CD_optimization_alternative)

	if print_verbose:
		# Verbose printing
		print("============================================================================")
		print("Curve fitting for D while K and q are known, based on Alternative {0}".format(CD_optimization_alternative))
		print("Estimated:\nK = {0}\nq = {1}\nD = {2}\nC = {3}".format(K_0, q_0, D1, u.KqD2C(K_0, q_0, D1)))

	# Mass and radius calculation based on the known K, q and D vaues
	MR = [s.mass_radius_relation(K=K_0, q=q_0, D=D1, rho_c=rho_c)
	for rho_c in rho_c_K]
	MR = np.array([*MR])
	M_D_fit = MR[:, 0]
	R_D_fit = MR[:, 1]
	###########################################################################

	########### THEORETICAL ASSESSMENT ########################################
	# Theoretical calculation of the constants

	q_theo = q_0 # It is known anyway :p
	D_theo = c.D # Theoretical value of D
	K_theo = u.CqD2K(c.C, q_theo, D_theo) # Conversion from C, q, D to K

	if print_verbose:
		print("============================================================================")
		print("Theoritical:\nK = {0}\nq = {1}\nD = {2}\nC = {3}".format(
		K_theo, q_theo, D_theo, c.C))

	# Mass and radius calculation based on the theoretical K, q and D vaues
	MR_theo = [s.mass_radius_relation(K=K_theo, q=q_theo, D=D_theo, rho_c=rho_c)
			for rho_c in rho_c_K]
	MR_theo = np.array([*MR_theo])
	M_theo = MR_theo[:, 0]
	R_theo = MR_theo[:, 1]
	###########################################################################
	
	################# PLOTTING ################################################ 
	if show_fig: 
		p.figure_()
		p.draw(handle="loglog", x_list=[mass], y_list=[rho_c_Kq, rho_c_K], labels=["Curve Fit for $K$ and $q$", "Curve Fit for $K$"],
				xlabel="Mass [$M_\odot$]", ylabel="Central Density [$\dfrac{kg}{m^3}$]", title="Central Density Distributions", ls=" ")

		p.figure_()
		p.draw(handle="loglog", x_list=[mass[low_mass_ix]], y_list=[rho_c_Kq[low_mass_ix], rho_c_K[low_mass_ix]], labels=["Curve Fit for $K$ and $q$", "Curve Fit for $K$"],
				xlabel="Mass [$M_\odot$]", ylabel="Central Density [$\dfrac{kg}{m^3}$]", title="Central Density Distributions\n[Only Low Mass Stars]",ls=" ")

		
		p.figure_()
		p.draw(handle="scatter", x_list=[mass], y_list=[rho_c_Kq, rho_c_K], labels=["Curve Fit for $K$ and $q$", "Curve Fit for $K$"],
				xlabel="Mass [$M_\odot$]", ylabel="Central Density [$\dfrac{kg}{m^3}$]", title="Central Density Distributions")

		p.figure_()
		p.draw(handle="scatter", x_list=[mass[low_mass_ix]], y_list=[rho_c_Kq[low_mass_ix], rho_c_K[low_mass_ix]], labels=["Curve Fit for $K$ and $q$", "Curve Fit for $K$"],
				xlabel="Mass [$M_\odot$]", ylabel="Central Density [$\dfrac{kg}{m^3}$]", title="Central Density Distributions\n[Only Low Mass Stars]")

		p.figure_()
		p.draw(handle="loglog", x_list=[R_aer], y_list=[mass, M_Kq_fit/c.SOLAR_MASS, M_K_fit/c.SOLAR_MASS, M_D_fit/c.SOLAR_MASS, M_theo/c.SOLAR_MASS], labels=["Actual Data", "Curve Fit for $K$ and $q$", "Curve Fit for only $K$", "Curve Fit for only $D$", "Theoretical Result"],
				xlabel="Radius [Average Earth Radius]", ylabel="Mass [$M_\odot$]", title="Mass Distributions", ls=" ")

		p.figure_()
		p.draw(handle="loglog", x_list=[R_aer[low_mass_ix]], y_list=[mass[low_mass_ix], M_Kq_fit[low_mass_ix]/c.SOLAR_MASS, M_K_fit[low_mass_ix]/c.SOLAR_MASS, M_D_fit[low_mass_ix]/c.SOLAR_MASS, M_theo[low_mass_ix]/c.SOLAR_MASS], labels=["Actual Data", "Curve Fit for $K$ and $q$", "Curve Fit for only $K$", "Curve Fit for only $D$", "Theoretical Result"],
				xlabel="Radius [Average Earth Radius]", ylabel="Mass [$M_\odot$]", title="Mass Distributions\n[Only Low Mass Stars]", ls=" ")

		p.figure_()
		p.draw(handle="scatter", x_list=[R_aer], y_list=[mass, M_Kq_fit/c.SOLAR_MASS, M_K_fit/c.SOLAR_MASS, M_D_fit/c.SOLAR_MASS, M_theo/c.SOLAR_MASS], labels=["Actual Data", "Curve Fit for $K$ and $q$", "Curve Fit for only $K$", "Curve Fit for only $D$", "Theoretical Result"],
				xlabel="Radius [Average Earth Radius]", ylabel="Mass [$M_\odot$]", title="Mass Distributions")

		p.figure_()
		p.draw(handle="scatter", x_list=[R_aer[low_mass_ix]], y_list=[mass[low_mass_ix], M_Kq_fit[low_mass_ix]/c.SOLAR_MASS, M_K_fit[low_mass_ix]/c.SOLAR_MASS, M_D_fit[low_mass_ix]/c.SOLAR_MASS, M_theo[low_mass_ix]/c.SOLAR_MASS], labels=["Actual Data", "Curve Fit for $K$ and $q$", "Curve Fit for only $K$", "Curve Fit for only $D$", "Theoretical Result"],
				xlabel="Radius [Average Earth Radius]", ylabel="Mass [$M_\odot$]", title="Mass Distributions\n[Only Low Mass Stars]")


		p.show_()


if __name__ == "__main__":
	demo(1, 1)










############ DEPRECATED CODE FOR PLOTTING WITH SAVING #################

"""

	################# PLOTTING ################################################ 
	p.figure_()
	p.draw(handle="loglog", x_list=[mass], y_list=[rho_c_Kq, rho_c_K], labels=["Curve Fit for $K$ and $q$", "Curve Fit for $K$"],
			xlabel="Mass [$M_\odot$]", ylabel="Central Density [$\dfrac{kg}{m^3}$]", title="Central Density Distributions", ls=" ", savefn="report/figures/{0}1_n_ll_rho_m.png".format(prefix))

	p.figure_()
	p.draw(handle="loglog", x_list=[mass[low_mass_ix]], y_list=[rho_c_Kq[low_mass_ix], rho_c_K[low_mass_ix]], labels=["Curve Fit for $K$ and $q$", "Curve Fit for $K$"],
			xlabel="Mass [$M_\odot$]", ylabel="Central Density [$\dfrac{kg}{m^3}$]", title="Central Density Distributions\n[Only Low Mass Stars]",
			ls=" ",  savefn="report/figures/{0}2_n_ll_rho_m_.png".format(prefix))

	
	p.figure_()
	p.draw(handle="scatter", x_list=[mass], y_list=[rho_c_Kq, rho_c_K], labels=["Curve Fit for $K$ and $q$", "Curve Fit for $K$"],
			xlabel="Mass [$M_\odot$]", ylabel="Central Density [$\dfrac{kg}{m^3}$]", title="Central Density Distributions", savefn = "report/figures/{0}3_n_s_rho_m.png".format(prefix))

	p.figure_()
	p.draw(handle="scatter", x_list=[mass[low_mass_ix]], y_list=[rho_c_Kq[low_mass_ix], rho_c_K[low_mass_ix]], labels=["Curve Fit for $K$ and $q$", "Curve Fit for $K$"],
			xlabel="Mass [$M_\odot$]", ylabel="Central Density [$\dfrac{kg}{m^3}$]", title="Central Density Distributions\n[Only Low Mass Stars]", savefn = "report/figures/{0}4_n_s_rho_m_.png".format(prefix))
	

	'''
	p.figure_()
	p.draw(handle="loglog", x_list=[R_aer], y_list=[mass_kg, M_Kq_fit, M_K_fit, M_D_fit, M_theo], labels=["Actual Data", "Curve Fit for $K$ and $q$", "Curve Fit for only $K$", "Curve Fit for only $D$", "Theoretical Result"],
			xlabel="Radius [Average Earth Radius]", ylabel="Mass [kg]", title="Mass Distributions", ls=" ", savefn = "report/figures/5_n_ll_m_r.png")

	p.figure_()
	p.draw(handle="loglog", x_list=[R_aer[low_mass_ix]], y_list=[mass_kg[low_mass_ix], M_Kq_fit[low_mass_ix], M_K_fit[low_mass_ix], M_D_fit[low_mass_ix], M_theo[low_mass_ix]], labels=["Actual Data", "Curve Fit for $K$ and $q$", "Curve Fit for only $K$", "Curve Fit for only $D$", "Theoretical Result"],
			xlabel="Radius [Average Earth Radius]", ylabel="Mass [kg]", title="Mass Distributions\n[Only Low Mass Stars]", ls=" ", savefn = "report/figures/6_n_ll_m_r_.png")
	'''

	p.figure_()
	p.draw(handle="loglog", x_list=[R_aer], y_list=[mass, M_Kq_fit/c.SOLAR_MASS, M_K_fit/c.SOLAR_MASS, M_D_fit/c.SOLAR_MASS, M_theo/c.SOLAR_MASS], labels=["Actual Data", "Curve Fit for $K$ and $q$", "Curve Fit for only $K$", "Curve Fit for only $D$", "Theoretical Result"],
			xlabel="Radius [Average Earth Radius]", ylabel="Mass [$M_\odot$]", title="Mass Distributions", ls=" ", savefn = "report/figures/{0}7_n_ll_ms_r.png".format(prefix))

	p.figure_()
	p.draw(handle="loglog", x_list=[R_aer[low_mass_ix]], y_list=[mass[low_mass_ix], M_Kq_fit[low_mass_ix]/c.SOLAR_MASS, M_K_fit[low_mass_ix]/c.SOLAR_MASS, M_D_fit[low_mass_ix]/c.SOLAR_MASS, M_theo[low_mass_ix]/c.SOLAR_MASS], labels=["Actual Data", "Curve Fit for $K$ and $q$", "Curve Fit for only $K$", "Curve Fit for only $D$", "Theoretical Result"],
			xlabel="Radius [Average Earth Radius]", ylabel="Mass [$M_\odot$]", title="Mass Distributions\n[Only Low Mass Stars]", ls=" ", savefn = "report/figures/{0}8_n_ll_ms_r_.png".format(prefix))

	'''
	p.figure_()
	p.draw(handle="scatter", x_list=[R_aer], y_list=[mass_kg, M_Kq_fit, M_K_fit, M_D_fit, M_theo], labels=["Actual Data", "Curve Fit for $K$ and $q$", "Curve Fit for only $K$", "Curve Fit for only $D$", "Theoretical Result"],
			xlabel="Radius [Average Earth Radius]", ylabel="Mass [kg]", title="Mass Distributions", savefn = "report/figures/9_n_s_m_r.png")

	p.figure_()
	p.draw(handle="scatter", x_list=[R_aer[low_mass_ix]], y_list=[mass_kg[low_mass_ix], M_Kq_fit[low_mass_ix], M_K_fit[low_mass_ix], M_D_fit[low_mass_ix], M_theo[low_mass_ix]], labels=["Actual Data", "Curve Fit for $K$ and $q$", "Curve Fit for only $K$", "Curve Fit for only $D$", "Theoretical Result"],
			xlabel="Radius [Average Earth Radius]", ylabel="Mass [kg]", title="Mass Distributions\n[Only Low Mass Stars]", savefn = "report/figures/10_n_s_m_r_.png")
	'''

	p.figure_()
	p.draw(handle="scatter", x_list=[R_aer], y_list=[mass, M_Kq_fit/c.SOLAR_MASS, M_K_fit/c.SOLAR_MASS, M_D_fit/c.SOLAR_MASS, M_theo/c.SOLAR_MASS], labels=["Actual Data", "Curve Fit for $K$ and $q$", "Curve Fit for only $K$", "Curve Fit for only $D$", "Theoretical Result"],
			xlabel="Radius [Average Earth Radius]", ylabel="Mass [$M_\odot$]", title="Mass Distributions", savefn = "report/figures/{0}11_n_s_ms_r.png".format(prefix))

	p.figure_()
	p.draw(handle="scatter", x_list=[R_aer[low_mass_ix]], y_list=[mass[low_mass_ix], M_Kq_fit[low_mass_ix]/c.SOLAR_MASS, M_K_fit[low_mass_ix]/c.SOLAR_MASS, M_D_fit[low_mass_ix]/c.SOLAR_MASS, M_theo[low_mass_ix]/c.SOLAR_MASS], labels=["Actual Data", "Curve Fit for $K$ and $q$", "Curve Fit for only $K$", "Curve Fit for only $D$", "Theoretical Result"],
			xlabel="Radius [Average Earth Radius]", ylabel="Mass [$M_\odot$]", title="Mass Distributions\n[Only Low Mass Stars]", savefn = "report/figures/{0}12_n_s_ms_r_.png".format(prefix))


	p.show_()
"""


############# ANOTHER DEPRECATED PLOT CODE ############################

"""
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
"""