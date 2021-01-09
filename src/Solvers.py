import numpy as np
from numpy.linalg.linalg import solve
from scipy.optimize import curve_fit
from scipy.integrate import solve_ivp
from . import Constants as c
def lane_emden(xi, theta, n):
    """
    docstring
    """
    if xi != 0:
        return np.array([theta[1],
                    (-xi * theta[0] ** n - 2 * theta[1]) / (xi)])
    else:
        return [theta[1], 0]


def solve_lane_emden(n):
    """
    docstring
    """
    print(n)
    # Condition for terminating the iteration when ==> theta(xi_n) = 0
    def is_surface(xi, theta, n): return theta[0]
    is_surface.terminal = True

    theta_0 = [1, 0] # Initial value for theta
    xi_span = [0, 5] # Initial value for xi and possible span

    r = solve_ivp(fun = lane_emden, t_span = xi_span, y0 = theta_0, events = is_surface, args = (n,))

    # If propagation is successful:
    if r.status >= 0:
        return r.t, r.y
    else:
        print("It seems solver is not able to obtain a solution!")

def white_dwarf_fit(M, R):
    """
    docstring
    """
    C_0 = 1.0
    q_0 = 3.0
    D_0 = 1.0
    y_0 = [C_0, q_0, D_0]

    popt, _ = curve_fit(f = mass_radius_relation, xdata = R, ydata = M, p0 = y_0)

    return popt


def mass_radius_relation(R, C, q, D):
    """
    docstring
    """
    # Simplification of the EOS for cold WDs
    n = q/(5-q)
    K = (8 * C) / (5 * D ** (5 / q))

    # Lane Emden solution for given n
    xi, theta = solve_lane_emden(n)
    xi_star = xi[-1]
    theta_dot_xi_star = theta[1, -1]

    # Calculation of mass w.r.t. radius, and etc.
    c1 = (((K * (n + 1)) / (4 * np.pi * c.G)) ** (-0.5 * ((3 - n) / (1 - n))))
    c2 = (((1 / (4 * np.pi)) * (((4 * np.pi * c.G) / (K * (n + 1))) ** 1.5) * (1 / ((xi_star ** 2) * (-theta_dot_xi_star)))) **(-1)) 
    c3 = (xi_star ** ((n - 3) / (1 - n)))
    c4 = (R ** ((3 - n) / (1 - n)))
    M = c1 * c2 * c3 * c4

    return M