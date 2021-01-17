import numpy as np
from . import newtonian_approach, relativistic_approach


if __name__ == "__main__":
    np.warnings.filterwarnings('ignore')
    print("NEWTONIAN APPROACH RESULTS:")
    newtonian_approach.demo()

    print("RELATIVISTIC APPROACH RESULTS:")
    relativistic_approach.demo()
