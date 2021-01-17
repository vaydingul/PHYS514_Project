
from . import Constants as c
from . import Solvers as s
import numpy as np
import csv


def read_dw_data(fname, out_type="np"):
    """
    read_dw_data:
    data = read_dw_data(fname, out_type = "np")

    This function perform the following processes:
        - It reads the ´White Dwarf´data
        - It outputs based on the specified data type 


    Input:
        fname = The directory of the file to be read
        out_type = Output type ==> NumPy array, or Python List

    Output:
        data = Output data, which is the concatenated form of the log(g) and mass

    Example:
        # Data reading and partitioning
        data = u.read_dw_data(DATA_PATH)
        logg = data[:,0]; mass = data[:,1];
    """
    data = []  # Data list initialization
    with open(fname, "r") as f:
        # More accurate control on the file via DictReader
        reader = csv.DictReader(f, delimiter=",")

        for row in reader:
            # Read each line
            data.append([float(row["logg"]), float(row["mass"])])

    # Return based on the desired output
    if out_type == "np":
        # np.reshape(np.array(data), (len(data), 2))
        return np.array(data)

    else:

        return data


def wd_data_conversion(val_mass, val_logg):
    """
    wd_data_conversion:
    mass_kg, R_m, R_aer = wd_data_conversion(mass, logg)

    This function perform the following processes:
        - It converts the log(g) data to radius in ´meters´
            and ´average earth radius´
        - It converts mass in ´solar mass´ to ´kg´

    Input:
        val_mass = Mass input  (in solar mass)
        val_logg = Log(g) input

    Output:
        mass_kg = Mass in kilograms
        R_m = Radius in meters
        R_aer = Radius in average earth radius

    Example:

        # Data reading and partitioning
        data = u.read_dw_data(DATA_PATH)
        logg = data[:,0]; mass = data[:,1]; mean_mass = np.mean(mass)
        mass_kg, R_m, R_aer = u.wd_data_conversion(mass, logg)


    """

    g = 10**(val_logg) * 1e-2  # Gravitational acceleration in SI, where
    # (1e-2) is used as a conversion factor from CGS to SI

    mass_kg = c.SOLAR_MASS * val_mass  # Mass conversion from solar mass to kg

    R_m = np.sqrt(c.G * mass_kg / g)  # Radius of the star in meters

    R_aer = R_m / c.AER  # Radius of the star in terms of Average Earth Radii

    return mass_kg, R_m, R_aer

def polytropic_conversion(K, n, p = None, rho = None):
    """
    polytopic_conversion:
    polytropic_conversion(K, n, p = None, rho = None)

    This function perform the following processes:
        - It calculate pressure or density based on the asked
            quantity
        - If pressure is given, it calculates density based on
            polytropic relation
        - If density is given, it calculates pressure based on 
            polytropic relation

    Input:
        K = Constants
        n = Polytropic index
        p = Pressure
        rho = Density

    Output:
        Pressure or density

    Example:

        []
    """
    
    if p is None:
        return K * ((rho) ** (1 + (1/n)))
    else:
        return (np.abs(p) / K) ** (n / (n+1))


def geometric_SI_conversion(time = None, mass = None, length = None, pressure = None, density = None, mode = 1):
    """
    geometric_SI_conversion:
    geometric_SI_conversion(time = None, mass = None, length = None, pressure = None, density = None, mode = 1)

    This function perform the following processes:
        - It applies unit conversion on the given quantity
        - It can vconvert from SI to geometric unit or
            vice versa
       
    Input:
        time = Time quantity
        mass = Mass quantity
        length = Length quantity
        pressure = Pressure quantity
        density = Density quantity
        mode = Conversion direction
            if mode == 1, it asserts the conversion from SI to geometric unit
            if mode == -1, it asserts the conversion from geometric unit to SI

    Output:
        Converted quantity

    Example:

        []
    """
    if time is not None:
        
        factor = ((c.G * c.SOLAR_MASS) / (c.SPEED_OF_LIGHT ** 3))
        return time / (factor ** mode)

    if mass is not None:

        factor = c.SOLAR_MASS
        return mass / (factor ** mode)  

    if length is not None:

        factor = ((c.G * c.SOLAR_MASS) / (c.SPEED_OF_LIGHT ** 2))
        return length / (factor ** mode)

    if pressure is not None:

        factor = (((c.G ** 3) * (c.SOLAR_MASS ** 2)) / (c.SPEED_OF_LIGHT ** 8))
        return pressure * (factor ** mode)

    if density is not None:

        factor = (((c.G ** 3) * (c.SOLAR_MASS ** 2)) / (c.SPEED_OF_LIGHT ** 6))
        return density * (factor ** mode)




##### SOME HELPER FUNCTIONS ##########

# It converts given ´q value to ´n´ counterpars
def q2n(q): return q / (5-q)

# It converts the given ´C´, ´q´, ´D´ values
# to ´K´ counterpart
def CqD2K(C, q, D): return (8 * C) / (5 * (D ** (5 / q)))

# It converts the given ´K´, ´q´, ´D´ values
# to ´C´ counterpart
def KqD2C(K, q, D): return (5 / 8) * K * (D ** (5 / q))

def Aq2K(A, q):
    """
    It converts given ´A´ and ´q´ values
    to ´K´ counterpart.
    """
    n = q2n(q)
    xi, theta = s.solve_lane_emden(n)
    # The following are required for the mass calculation
    xi_star = xi[-1]
    theta_dot_xi_star = theta[1, -1]

    B = (((n + 1) / (4 * np.pi * c.G)) ** (((n - 3) / (2 - 2 * n)) + 1.5)) * ((xi_star) ** (((n - 3) / (1 - n)) + 2)) * 4 * np.pi * (-theta_dot_xi_star)

    K = (A / B) ** ((n - 1) / (n))

    return K

