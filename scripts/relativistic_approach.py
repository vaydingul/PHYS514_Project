# Ready-to-use module import

from scipy.interpolate import interp1d
from scipy.constants import speed_of_light
import numpy as np
import sys
sys.path.insert(1, "./")
sys.path.insert(2, "./../")

# Custom module import
from src import Utils as u
from src import Solvers as s
from src import Plotter as p
from src import Constants as c

def demo(N1 = 100, N2 = 200 ):

    """
    demo:
    demo(N1 = 100, N2 = 200)

    This function perform the following processes:
        - It visualizes the results obtained in relativistic approach


    Input:
        N1 = Number of sampling points for polytropic constant (K_SN)
        N2 = Number of sampling points for central pressure (p_c)

    Output:
        []

    Example:
        []
    """


    # Number of sampling points for K_SN
    N1 = 200
    # Number of sampling points for central pressure
    N2 = 200
    # Polytopic constant
    K_NS_list = np.arange(100 - int(N1 / 2),100 + int(N1 / 2))
    K_ix = np.where(K_NS_list == 100) 
    # Array index for K_NS = 100

    # Polytropic index
    n = 1

    # Central pressure sampling in SI units
    p_SI = np.logspace(20, 45, N2)
    # Central pressure sampling in geometric units
    p_gu = u.geometric_SI_conversion(pressure = p_SI, mode = 1)

    # Mass, radius, and rest mass matrix preallocation
    M = np.empty((N1, N2))
    R = np.empty((N1, N2))
    M_p = np.empty((N1, N2))

    # Iteration for different central pressures
    for (ix1, K_NS) in enumerate(K_NS_list):
        print(ix1)
        for (ix2, p_) in enumerate(p_gu):

            # Solution of the TOV equations
            r, soln = s.solve_tov(p_, K_NS, n)

            # Get the latest element of mass, radius and rest mass solutions
            M[ix1, ix2] = soln[0, -1]
            M_p[ix1, ix2] = soln[-1, -1]
            R[ix1, ix2] = r[-1]

    # Calculation of the maximum masses for each polytopic constant setting
    M_max = np.max(M, axis = 1)
    # Fractional binding energy
    Delta = (M_p - M) / M 

    #############################################################
    # In here, it is assumed that the polytropic constant
    # depends on the units that we are working, so this 
    # long path for conversion is chosen.

    # Pressure to density conversion
    rho_gu = u.polytropic_conversion(K_NS_list[K_ix], n, p = p_gu)
    # Geometric unit to SI unit conversion
    rho_SI = u.geometric_SI_conversion(density=rho_gu, mode = -1)
    R_km = u.geometric_SI_conversion(length = R, mode = -1) / 1000
    #############################################################


    # Derivative calculation of the mass w.r.t. central density
    # to calculate the stability characteristics
    dM_drho = np.gradient(M[K_ix, :], rho_SI)
    stable_ixs = dM_drho >= 0.0
    unstable_ixs = dM_drho < 0.0
    # Allocation of the stable and unstable regions
    M_stable, rho_SI_stable, M_unstable, rho_SI_unstable = M[K_ix, stable_ixs], rho_SI[stable_ixs], M[K_ix, unstable_ixs], rho_SI[unstable_ixs]

    M_max_allowable_ix = np.where(np.isclose(M_max, 2.14, atol = 0.01))

    ################################## PLOTTING #####################################
    p.figure_()
    p.draw("scatter", [R_km[K_ix, :]], [M[K_ix, :]], xlabel = "Radius [km]", ylabel = "Mass [$M_\odot$]", title = "Mass-Radius Relation\n$K_{NS} = ${0}".format(K_NS_list[K_ix]), savefn = "report/figures/17_e_s_m_r.png")

    p.figure_()
    p.draw("scatter", [R_km[K_ix, :]], [Delta[K_ix, :]], xlabel = "Radius [km]", ylabel = "$\Delta$", title = "Fractional Binding Energy and Radius Relation\n$K_{NS} = ${0}".format(K_NS_list[K_ix]), savefn = "report/figures/18_e_s_delta_r.png")

    #p.figure_()
    #p.draw("scatter", [rho_SI], [M], xlabel = "Central Density [$\dfrac{kg}{m^3}$]", ylabel = "Mass [$M_\odot$]", title = "Mass and Central Density Relation")

    p.figure_()
    p.draw("loglog", [rho_SI_stable, rho_SI_unstable], [M_stable, M_unstable], labels = ["Stable", "Unstable"], xlabel = "Central Density [$\dfrac{kg}{m^3}$]", ylabel = "Mass [$M_\odot$]", title = "Mass and Central Density Relation\n$K_{NS} = ${0}".format(K_NS_list[K_ix]), savefn = "report/figures/19_e_ll_m_rho.png")

    p.figure_()
    p.draw("scatter", [K_NS_list], [M_max], xlabel = "$K_{NS}$", ylabel = "$M_{max}$ [$M_\odot$]", title = "Maximum Mass and Polytopic Constant Relation")
    p.draw("plot", [np.array([K_NS_list[0], K_NS_list[-1]])], [np.array([M_max[M_max_allowable_ix][0], M_max[M_max_allowable_ix][0]])], color = "red")
    p.draw("plot", [np.array([K_NS_list[M_max_allowable_ix][0], K_NS_list[M_max_allowable_ix][0]])], [np.array([M_max[0], M_max[-1]])], color = "red", savefn = "report/figures/20_e_s_m_k.png")

    p.show_()

    
if __name__ == "__main__":
    demo()