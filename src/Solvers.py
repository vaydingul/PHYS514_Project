import numpy as np
from numpy.linalg.linalg import solve
from scipy.optimize import curve_fit
from scipy.integrate import solve_ivp

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
    # Condition for terminating the iteration when ==> theta(xi_n) = 0
    def is_surface(xi, theta, n): return theta[0]
    is_surface.terminal = True

    theta_0 = [1, 0] # Initial value for theta
    xi_span = [0, 5] # Initial value for xi and possible span

    r = solve_ivp(fun = lane_emden, t_span = xi_span, y0 = theta_0, events = is_surface, args = (n,), max_step = 1e-2)

    # If propagation is successful:
    if r.status >= 0:
        return r.t, r.y
    else:
        print("It seems solver is not able to obtain a solution!")

