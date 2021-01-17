import numpy as np
import newtonian_approach, relativistic_approach


if __name__ == "__main__":
    
    # This warning suppression is to display
    # clear output, since all the irregularities in the 
    # code are dealt with. However, one can simply ignore
    # that line by commenting it out.
    np.warnings.filterwarnings('ignore')


    # This part takes approximately the 1 minutes
    # to complete. At the end of the run, several
    # graphs will be displayed. If it is not desired, 
    # one can simply alter this option, by providing
    # input argument as ==> ´show_fig = False´
    # For more detailed input list, please
    # refer to the ´newtonian_approach.py´ script.
    print("NEWTONIAN APPROACH RESULTS:")
    newtonian_approach.demo()

    # ! PLEASE CLOSE ALL THE FIGURES THAT ARE OPENED
    # ! SO FAR, FOR CODE ITSELF TO PROCEED.

    # This part takes approximately the 2 minutes
    # to complete. At the end of the run, several
    # graphs will be displayed. In this part,
    # there is no option for suppressing the
    # figures, since there is only 4 figures to 
    # display. To get information about input arguments
    # , please refer to the ´relativistic_approach.py´ script.
    print("RELATIVISTIC APPROACH RESULTS:")
    relativistic_approach.demo()
