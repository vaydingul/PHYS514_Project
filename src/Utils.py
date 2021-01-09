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
    
    AER = 6371e3 # Average Earth radii in meters
    SOLAR_MASS = 1.98847e30 # Mass of the Sun
    EARTH_MASS = 5.9722e24 # Mass of the Earth
    G = 6.67408e-11 # Gravitational constant

    g = 10**(val_logg) * 1e-2 # Gravitational acceleration in SI, where 
    # (1e-2) is used as a conversion factor from CGS to SI 

    mass_kg = SOLAR_MASS * val_mass # Mass conversion from solar mass to kg
    
    R_m = np.sqrt(G * mass_kg / g) # Radius of the star in meters

    R_aer = R_m / AER # Radius of the star in terms of Average Earth Radii

    return mass_kg, R_m, R_aer 


    
    