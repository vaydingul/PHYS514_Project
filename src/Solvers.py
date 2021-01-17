import numpy as np
from numpy.linalg.linalg import solve
from scipy.optimize import curve_fit, minimize, minimize_scalar, root
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp
from . import Constants as c
from . import Utils as u

def mass_equation(r, rho):
    """
    mass_equation:

    mass_equation(r, rho)

    This function perform the following processes:
        - It propagates the mass ODE


    Input:
        r = Radius of the star
        rho = Density of the star 

    Output:

        d m(r) / d r = 4 * pi * r^2 * rho(r)

    Example:
        []



    """
    return 4 * np.pi * (r ** 2) * rho


def pressure_equation(r, m, rho):
    """
    pressure_equation:

    pressure_equation(r, m, rho)

    This function perform the following processes:
        - It propagates the pressure ODE


    Input:
        r = Radius of the star
        m = Mass of the star
        rho = Density of the star 

    Output:

        d p(r) / d r = - (G * m(r) * rho(r)) / (r^2)

    Example:
        []


    """
    if r != 0:
        return -(c.G * m * rho) / (r ** 2)
    else:
        return 0.0



def density_equation(r, rho, m, C, q, D, mode):
    """

    density_equation:

    density_equation(r, rho, m, C, q, D)

    This function perform the following processes:
        - It is basically the rewritten form the pressure equation, in terms of density.
        - It basically envision the following idea:

        dp / dr = (dp / dx) * (dx / drho) * (drho / dr)

    Input:
        r = Radius of the star
        rho = Density of the star 
        m = Mass of the star
        C, q, D = Constants
        mode = It specifies the optimization alternative which is discussed in report

    Output:

        drho / dr = (dp / dr) / ((dp / dx) * (dx / drho))

    Example:
        []



    """
    x = (rho/D) ** (1 / q)

    # Following s1, s2, s3 and s4 contains the derivative of the pressure density relation w.r.t. x
    if mode == 1:
        # Based on alternative 1
        s1 = ((2 * x ** 2 - 3) * ((x ** 2 + 1) ** 0.5))
        s2 = ((x) * (4 * x) * ((x ** 2 + 1) ** 0.5))
        s3 = ((x) * (2 * x ** 2 - 3) * (0.5) * (2 * x) * ((x ** 2 + 1) ** -0.5))
        s4 = ((3) / (np.sqrt(x ** 2 + 1)))
        dpdx = C * (s1 + s2 + s3 + s4)
    else:
        # Based on Alternative 2
        dpdx = 8 * C * (x ** 4)


    # Trivial differentiation
    dxdrho = (D ** (-1/q)) * (rho ** ((1/q) - 1)) / (q)

    # One step calculation of the pressure equation
    dpdr = pressure_equation(r, m, rho)

    # Finally, drho / dr
    if dpdx != 0 and dxdrho != 0:
        drhodr = dpdr / (dpdx * dxdrho)
    else:
        drhodr = 0

    return drhodr




def  mass_density_equation(r, m_rho, C, q, D, mode):
    """

    mass_density_equation:

    mass_density_equation(r, m_rho, C, q ,D)

    This function perform the following processes:
        - System of ODE representation of the mass and density equation, 
            to be able to solved concurrently


    Input:
        r = Radius of the star
        m_rho = Joint representation of the propagated variables 
        C, q, D = Constants
        mode = It specifies the optimization alternative which is discussed in report

    Output:

        d m(r) / d r = 4 * pi * r^2 * rho(r)
        drho / dr = (dp / dr) / ((dp / dx) * (dx / drho))

    Example:
        []

    """
    # Separation of variables
    m = m_rho[0]
    rho = m_rho[1]
    # nan_to_num function is used to surpass irregularities ( but, not quite sure :p)
    return [np.nan_to_num(mass_equation(r, rho)), np.nan_to_num(density_equation(r, rho, m, C, q, D, mode))]
    
    # I could not figure out why, but below option is not working :(
    #return np.nan_to_num(np.array([mass_equation(r, rho), density_equation(r, rho, m, C, q, D)]))



def lane_emden(xi, theta, n):
    """
    lane_emden:

    lane_emden(xi, theta, n)

    This function perform the following processes:
        - It propagates the Lane-Emden equation


    Input:
        xi = ODE independent variable
        theta = ODE dependent variable 

    Output:

        (1 / xi^2) * (d/d xi) * (xi^2 * (d theta / d xi)) + theta^n = 0

    Example:
        []
    """
    if xi != 0:
        return np.array([theta[1],
                         (-xi * np.power(theta[0], n, out=np.zeros_like(theta[0]), where=theta[0] > 0) - 2 * theta[1]) / (xi)])
    else:
        return [theta[1], 0]


def tov_equation(r, s, K, n):
    """
    tov_equation:

    tov_equation(r, mvp, K, n)

    This function perform the following processes:
        - It propagates the Tolman-Oppenheimer-Volkoff equation


    Input:
        r = Radius
        s = Mass - Time Dilation - Pressure - Rest Mass
        K = Polytropic constant
        n = Polytropic index
    Output:

        m' = 4 * pi * r^2 * rho
        v' = 2 * ((m + 4 * pi * r^3 * p) / (r * (r - 2*m)))
        p' = -((m + 4 * pi * r^3 * p) / (r * (r - 2*m))) * (p + rho)
        m_p' = 4 * pi * (1 - ((2 * m)/(r)))^(-0.5)  * r^2 * rho
    Example:
        []
    """
    m, v, p, m_p = s
    #print(r,m,v,p, K, n)

    # Write density in terms of pressure
    rho = u.polytropic_conversion(K, n, p = p)
    
    m_next = 4 * np.pi * (r ** 2) * rho

    n_ = (m + 4 * np.pi * (r ** 3) * p)
    d_ = (r * (r - 2 * m))
    # Check if denominator is zero or not
    if d_ != 0:
        v_next = 2 * (n_ / d_)
    else:
        v_next = 0
    
    p_next = -0.5 * (p + rho) * v_next 

    if r != 0:
        m_p_next = 4 * np.pi * ((1 - ((2 * m) / (r))) ** (-0.5)) * (r ** 2) * rho 
    else:
        m_p_next = 0.0

    return (np.array([m_next, v_next, p_next, m_p_next]))
    

def pressure_density_relation(rho, C, q, D):
    """

    pressure_density_relation:

    pressure_density_relation(rho, C, q , D)

    This function perform the following processes:
        - It calculates the pressure value based on the density cold white dwarfs


    Input:
        rho = Density of the star
        C, q ,D = Constants

    Output:

        P = C [ x * (2 * x^2 - 3) * ((x^2 + 1)^0.5) + 3 * sinh^-1(x)]

        where

        x = (rho / D) ^ (1 / q)

    Example:
        []



    """

    x = (rho/D) ** (1 / q)
    return C * (x * (2 * x ** 2 - 2) * ((x ** 2 + 1) ** 0.5) + 3 * (np.arcsinh(x)))


def calculate_rho_c(M=None, R=None, K=None, n=None):
    """
    calculate_rho_c:

    calculate_rho_c(M, K, n)

    This function perform the following processes:
        - It calculates the central density of a star for a given
        mass, K and n value


    Input:
        M = Mass of the star
        K = Constant
        n = Polytropic constant 

    Output:

        rho_c = Central density value which is calculated by the equation in the 
                project report

    Example:
        []
    """
    # Lane Emden solution for given n
    xi, theta = solve_lane_emden(n)
    xi_star = xi[-1]
    theta_dot_xi_star = theta[1, -1]

    # It determines base of the calculation
    # M or R
    if R is None:
        return ((M / (4 * np.pi)) * (((4 * np.pi * c.G) / (K * (n + 1))) ** 1.5) * (1 / ((xi_star ** 2) * (-theta_dot_xi_star)))) ** ((2 * n) / (3 - n))
    else:
        return ((R) / ((((K * (n + 1)) / (4 * np.pi * c.G)) ** 0.5) * (xi_star))) ** ((2 * n) / (1 - n))


def white_dwarf_fit(M, R, K=None, C=None, q=None, D=None, A = None, rho_c_list=None, type="Kq", mode = 1):
    """
    white_dwarf_fit:

    white_dwarf_fit(M, R, K=None, C=None, q=None, D=None, rho_c_list = None, type = "Kq")

    This function perform the following processes:
        - It calculates the required parameters as a result of the curve fitting problem
        - It handles the problem in a generic way


    Input:
        M = Mass vector to be fitted
        R = Radius vector to be fitted
        K, C, q, D, A = Constant
        rho_c_list = Rho_c sweep for D calculation
        type = It specifies which parameters are optimized
        mode = It specifies the optimization alternative which is discussed in report

    Output:

        popt = Solution ot curve fitting problem

    Example:
        []
    """
    
    if (K is not None) or (A is not None):

        if q is not None:

            if D is None:
                
                if type == "Kq":

                    if mode == 2:
                        
                        # If K and q are to be optimized

                        # Initial guess for K and q
                        y_0 = [K, q]

                        # Monkey patching
                        def mrr_(R, K, q): return mass_radius_relation(
                            R=R, K=K, q=q)

                    else: 
                        # If A and q are to be optimized
                        y_0 = [A, q]
                        def mrr_(R, A, q): return mass_radius_relation(R = R, q = q, A = A)


                elif type == "K":
                    if mode == 2:

                        # If K is to be optimized for a given q value

                        # Initial guess for K
                        y_0 = [K]

                        # Monkey patching
                        def mrr_(R, K): return mass_radius_relation(R=R, K=K, q=q)
                    else:

                        # If A is to be optimized for a given q
                        y_0 = [A]
                        def mrr_(R, A): return mass_radius_relation(R = R, q = q, A = A)
            else:

                if type == "D":

                    # If D is to be optimized for a given K and q value

                    # Initial guess for D
                    y_0 = [D]

                    # Monkey patching - 1
                    def mrr(R, K, q, D, mode):

                        MR = [mass_radius_relation(
                            K=K, q=q, D=D, rho_c=rho_c, mode = mode) for rho_c in rho_c_list]
                        MR = np.array([*MR])
                        M_ = MR[:, 0]
                        R_ = MR[:, 1]
                        sp_f = interp1d(x=R_, y=M_, kind="cubic",
                                        fill_value="extrapolate")
                        M_spline = sp_f(R)
                        return M_spline

                    # Monkey patching to Monkey Patching - 1
                    def mrr_(R, D): return mrr(R=R, K=K, q=q, D=D, mode = mode)
                

    
    # Solution of the curve-fitting problem
    popt, _ = curve_fit(f=mrr_, xdata=R, ydata=M, p0=y_0)

    return popt


def mass_radius_relation(R=None, K=None, C=None, q=None, D=None, rho_c=None, A = None, mode = 1):
    """
    mass_radius_relation:

    mass_radius_relation(R, K=None, C=None, q=None, D=None, rho_c = None, )

    This function perform the following processes:
        - It calculates the massRadius or both of a star for a given constant values
        - It handles the problem in a generic way


    Input:
        R = Radius vector to be fitted
        K, C, q, D, A = Constant
        rho_c = Rho_c
        mode = It specifies the optimization alternative which is discussed in report

    Output:

        M = Mass of a star
        R = Radius of a star

    Example:
        []
    """

    if R is not None:

        if K is not None:

            if q is not None:

                # If K and q are given
                n = q/(5-q)
                # Solution of Lane-Emden equation
                xi, theta = solve_lane_emden(n)
                # The following are required for the mass calculation
                xi_star = xi[-1]
                theta_dot_xi_star = theta[1, -1]

                # Low level function
                return mass_radius_relation_(R, K, n, xi_star, theta_dot_xi_star)

            else:
                
                print("q should be defined!")

        else:

            if A is not None:

                # If K and q are given
                n = q/(5-q)

                return _mass_radius_relation(R, A, n)

                # Solution of Lane-Emden equation
                #xi, theta = solve_lane_emden(n)
                # The following are required for the mass calculation
                #xi_star = xi[-1]
                #theta_dot_xi_star = theta[1, -1]
                # Low level function
                


    else:

        if K is not None:

            if q is not None:

                if D is not None:

                    if rho_c is not None:
                        # If K, q, D and rho_c is given

                        # Calculation of C from K, q, and D
                        C = 5 * K * (D ** (5 / q)) / 8
                        # Mass and Density ODE solution for a given
                        # rho_c and C, q, and D
                        r, m_rho = solve_m_rho(rho_c, C, q, D, mode)

                        # Low level function
                        return mass_radius_relation__(r, m_rho)

                    else:

                        print("rho_c should be defined!")

                else:

                    print("D should be defined!")

            else:

                print("q should be defined!")

        else:

            print("K should be defined!")


def mass_radius_relation_(R, K, n, xi_star, theta_dot_xi_star):
    """
        It calculates the mass based on the Lane-Emden equation
    """
    # Calculation of mass w.r.t. radius, and etc.
    c1 = (((K * (n + 1)) / (4 * np.pi * c.G)) ** (-0.5 * ((3 - n) / (1 - n))))
    c2 = (((1 / (4 * np.pi)) * (((4 * np.pi * c.G) / (K * (n + 1))) ** 1.5)
           * (1 / ((xi_star ** 2) * (-theta_dot_xi_star)))) ** (-1))
    c3 = (xi_star ** ((n - 3) / (1 - n)))
    c4 = (R ** ((3 - n) / (1 - n)))
    M = c1 * c2 * c3 * c4
    return M


def mass_radius_relation__(r, m_rho):
    """
        It calculates the mass and radius
        via using the integrated mass and density ODE,
        until rho(R) = 0
    """

    # It cal
    R = r[-1]
    M = m_rho[0, -1]

    return M, R

def _mass_radius_relation(R, A, n):
    """
        It calculates the mass based on the proportional
        form of mass-radius relation:

        M = A * R ^ ((3 - n) / (1 - n))
    """

    M  = A * (R ** ((3 - n) / (1 - n)))
    return M 

def solve_m_rho(rho_c, C, q, D, mode):
    """
    solve_m_rho:

    solve_m_rho(rho_c, C, q, D)

    This function perform the following processes:
        - It solves the mass density ODE set


    Input:
        rho_c = Initial value of density equation
        C, q, D = Constants arguments
        mode = It specifies the optimization alternative which is discussed in report

    Output:

        r = Propagated radius span
        m_rho = Propagated mass and density span

    Example:
        []
    """
    # Termination criteria for the surface
    def is_surface(r, m_rho, C, q, D, mode): return m_rho[1] - 1
    is_surface.terminal = True
    is_surface.direction = 0

    # Initial value for mass and density:
    m_rho_0 = [0, rho_c]
    # Possible radius span
    r_span = [0, 1e15]

    # Solution of the IVP
    r = solve_ivp(fun=mass_density_equation, t_span=r_span,
                  y0=m_rho_0, events=is_surface, args=(C, q, D, mode))

    if r.status >= 0:
        return r.t, r.y
    else:
        print("It seems solver is not able to obtain a solution!")



def solve_tov(p_c, K, n):
    """
    solve_tov:

    solve_tov(p_c, K, n)

    This function perform the following processes:
        - It solves the TOV equation


    Input:
        p_c = Initial central pressure value
        K, n = Constants arguments

    Output:

        r = Propagated radius span
        s = Solution vector in which includes:
            Mass - Time Dilation - Pressure - Rest Mass

    Example:
        []
    """
    # Termination criteria for the surface
    def is_surface(r, s, K, n): return s[-2] - 1e-50
    is_surface.terminal = True
    is_surface.direction = 0

    # Initial value for Mass - Time Dilation - Pressure - Rest Mass:
    s_0 = [0, 0, p_c, 0]
    # Possible radius span
    r_span = [0, 1e15]

    # Solution of the IVP
    r = solve_ivp(fun=tov_equation, t_span=r_span,
                  y0=s_0, events=is_surface, args=(K, n), max_step = 1)

    if r.status >= 0:
        return r.t, r.y
    else:
        print("It seems solver is not able to obtain a solution!")



def solve_lane_emden(n):
    """
    solve_lane_emden:

    solve_lane_emden(n)

    This function perform the following processes:
        - It solves the Lane-Emden equation


    Input:
        n = Polytropic index

    Output:

        xi = Propagated xi span
        theta = Propagated theta and theta_dot span

    Example:
        []
    """

    # Condition for terminating the iteration when ==> theta(xi_star) = 0
    def is_surface(xi, theta, n): return theta[0]
    is_surface.terminal = True
    is_surface.direction = 0

    theta_0 = [1, 0]  # Initial value for theta
    xi_span = [0, 10]  # Initial value for xi and possible span

    r = solve_ivp(fun=lane_emden, t_span=xi_span,
                  y0=theta_0, events=is_surface, args=(n,))

    # If propagation is successful:
    if r.status >= 0:
        return r.t, r.y
    else:
        print("It seems solver is not able to obtain a solution!")



##### DEPRECATED CODE #######
'''
def white_dwarf_fit(M, R):
    """
    docstring
    """
    C_0 = 6e22
    q_0 = 3.0
    D_0 = 9.82e5
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


    ==> C, q, D = 1.7411410126411223e+18 3.3436716763923866 11154154.638198024
'''
