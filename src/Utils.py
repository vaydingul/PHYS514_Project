
from . import Constants as c
import numpy as np
import csv

def read_dw_data(fname, out_type = "np"):
    """
    docstring
    """
    data = [] # Data list initialization 
    with open(fname, "r") as f:
        # More accurate control on the file via DictReader
        reader = csv.DictReader(f, delimiter = "," ) 
        
        for row in reader:
            # Read each line
            data.append([float(row["logg"]), float(row["mass"])])
    
    # Return based on the desired output
    if out_type == "np":
        # np.reshape(np.array(data), (len(data), 2))
        return np.array(data)

    else:

        return data


def logg_2_aer(val_mass, val_logg):
    """
    docstring
    """
    
    

    g = 10**(val_logg) * 1e-2 # Gravitational acceleration in SI, where 
    # (1e-2) is used as a conversion factor from CGS to SI 

    mass_kg = c.SOLAR_MASS * val_mass # Mass conversion from solar mass to kg
    
    R_m = np.sqrt(c.G * mass_kg / g) # Radius of the star in meters

    R_aer = R_m / c.AER # Radius of the star in terms of Average Earth Radii

    return mass_kg, R_m, R_aer 


    
    