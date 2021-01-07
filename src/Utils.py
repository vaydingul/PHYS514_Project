import numpy as np
import csv

def read_dw_data(fname, out_type = "np"):
    """
    docstring
    """
    data = []
    with open(fname, "r") as f:
        reader = csv.DictReader(f, delimiter = "," )
        
        for row in reader:
            
            data.append([float(row["logg"]), float(row["mass"])])
    
    if out_type == "np":

        return np.array(data)

    else:

        return data
    