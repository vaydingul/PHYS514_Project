
from . import Constants as c
import numpy as np
import csv

def read_dw_data(fname, out_type = "np"):
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


def general_conversion(val_mass, val_logg):
    """
    general_conversion:
    mass_kg, R_m, R_aer = general_conversion(mass, logg)

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
        mass_kg, R_m, R_aer = u.general_conversion(mass, logg)


    """
    
    

    g = 10**(val_logg) * 1e-2 # Gravitational acceleration in SI, where 
    # (1e-2) is used as a conversion factor from CGS to SI 

    mass_kg = c.SOLAR_MASS * val_mass # Mass conversion from solar mass to kg
    
    R_m = np.sqrt(c.G * mass_kg / g) # Radius of the star in meters

    R_aer = R_m / c.AER # Radius of the star in terms of Average Earth Radii

    return mass_kg, R_m, R_aer 


# It convert given ´q value to ´n´ value :p 
def q2n(q): return q / (5-q)

    
    