# Ready-to-use module import

from scipy.interpolate import interp1d
import numpy as np
import sys
sys.path.insert(1, "./")
sys.path.insert(2, "./../")

from scripts import newtonian_approach

if __name__ == "__main__":
    np.warnings.filterwarnings('ignore')
    
    print("For M < 0.34 M")
    newtonian_approach.demo(show_fig = False, mass_constraint= 0.34)
    
    print("\n\n\n")

    print("For M < 0.5 M")
    newtonian_approach.demo(show_fig = False, mass_constraint=0.5)