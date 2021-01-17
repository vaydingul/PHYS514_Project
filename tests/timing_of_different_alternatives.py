



from scipy.interpolate import interp1d
import numpy as np
import sys
from timeit import timeit
sys.path.insert(1, "./")
sys.path.insert(2, "./../")

from scripts import newtonian_approach
from src import Plotter as p

if __name__ == "__main__":

    num_iter = 10
    np.warnings.filterwarnings('ignore')
    
    
    time_11 = timeit(lambda: newtonian_approach.demo(Kq_optimization_alternative = 1, CD_optimization_alternative = 1,
    show_fig = False, print_verbose = False), number = num_iter)
    
    time_12 = timeit(lambda: newtonian_approach.demo(Kq_optimization_alternative = 1, CD_optimization_alternative = 2,
    show_fig = False, print_verbose = False), number = num_iter)

    time_21 = timeit(lambda: newtonian_approach.demo(Kq_optimization_alternative = 2, CD_optimization_alternative = 1,
    show_fig = False, print_verbose = False), number = num_iter)

    time_22 = timeit(lambda: newtonian_approach.demo(Kq_optimization_alternative = 2, CD_optimization_alternative = 2,
    show_fig = False, print_verbose = False), number = num_iter)

    labels = ["Alternative 1 for\n($K-q$)\nAlternative 1 for\n($C-D$)",
    "Alternative 1 for\n($K-q$)\nAlternative 2 for\n($C-D$)",
    "Alternative 2 for\n($K-q$)\nAlternative 1 for\n($C-D$)",
    "Alternative 2 for\n($K-q$)\nAlternative 2 for\n($C-D$)"]

    p.figure_()
    p.draw("bar", [labels], [[time_11/num_iter, time_12/num_iter, time_21/num_iter, time_22/num_iter]], title = "Timings of Different Alternatives", 
        ylabel = "Average Elapsed Time [sec]")
    p.show_()