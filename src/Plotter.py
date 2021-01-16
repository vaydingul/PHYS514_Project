import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
handle_dict = {"scatter":plt.scatter, "plot":plt.plot, "semilogx":plt.semilogx, "semilogy":plt.semilogy,"loglog":plt.loglog}
# The list of markers to be used for each plot
markers = ["o", "^", "h", "s", ".", "p"]


def draw(handle, x_list, y_list, labels=None, xlabel=None, ylabel=None, title=None, savefn=None, **kwargs):
    """

    draw:
    draw(handle, x_list, y_list, labels = None, xlabel = None, ylabel = None, title = None, savefn = None)

    This function perform the following processes:
        - It executes a general plotting
        - It basically wraps the matplotlib 


    Input:
        x_list = The list of x inputs, it can be 1D or 2D
        y_list = The list of y inputs, it can be 1D or 2D
        labels = The list of labels for each plotting set
        xlabel = X label of the whole 1D set
        ylabel = Y label of the whole 1D set
        title  = Title of the whole 1D set
        savefn = If it is specified, it saves the figure to the specified location 

    Output:
        []

    Example:
        []

    """
    handle = handle_dict[handle]


    # Plot counter
    cnt = 0


    

    # If there are not individual X's for each Y
    if len(x_list) != len(y_list):
        x_list = [x_list[0] for _ in range(len(y_list))]



    for (x_, y_) in zip(x_list, y_list):
        handle(x_, y_, label=labels[cnt], marker=markers[cnt], **kwargs)
        cnt += 1


    if xlabel is not None:
        plt.xlabel(xlabel)
    if ylabel is not None:
        plt.ylabel(ylabel)
    if title is not None:
        plt.title(title)
    if labels is not None:
        plt.legend()
    
    plt.gca().xaxis.set_minor_formatter(NullFormatter())

    if savefn is not None:
        plt.savefig(savefn, bbox_inches="tight")


    plt.margins(0.0)
    plt.tight_layout()




'''


def plot_(x_list, y_list, labels=None, xlabel=None, ylabel=None, title=None, savefn=None):
    """

    plot_:
    plot_(x_list, y_list, labels = None, xlabel = None, ylabel = None, title = None, savefn = None)

    This function perform the following processes:
        - It executes a general plotting
        - It basically wraps the matplotlib 


    Input:
        x_list = The list of x inputs, it can be 1D or 2D
        y_list = The list of y inputs, it can be 1D or 2D
        labels = The list of labels for each plotting set
        xlabel = X label of the whole 1D set
        ylabel = Y label of the whole 1D set
        title  = Title of the whole 1D set
        savefn = If it is specified, it saves the figure to the specified location 

    Output:
        []

    Example:
        []

    """

    # The list of markers to be used for each plot
    markers = ["o", "^", "h", "s", ".", "p"]
    # Plot counter
    cnt = 0

    if labels is not None:
        # If there is individual X's for each Y
        if len(x_list) == len(y_list):
            for (x_, y_) in zip(x_list, y_list):
                plt.plot(x_, y_, label=labels[cnt], marker=markers[cnt])
                cnt += 1
        else:
            for y_ in y_list:
                plt.plot(x_list, y_, label=labels[cnt], marker=markers[cnt])
                cnt += 1

    # If labels are not provided, then do not take into account
    else:

        # If there is individual X's for each Y
        if len(x_list) == len(y_list):
            for (x_, y_) in zip(x_list, y_list):
                plt.plot(x_, y_, marker=markers[cnt])
                cnt += 1
        else:
            for y_ in y_list:
                plt.plot(x_list, y_, marker=markers[cnt])
                cnt += 1

    if xlabel is not None:
        plt.xlabel(xlabel)
    if ylabel is not None:
        plt.ylabel(ylabel)
    if title is not None:
        plt.title(title)

    plt.tight_layout()


def scatter_(x_list, y_list, labels=None, xlabel=None, ylabel=None, title=None, savefn=None):
    """
    scatter_:
    scatter_(x_list, y_list, labels = None, xlabel = None, ylabel = None, title = None, savefn = None)

    This function perform the following processes:
        - It executes a general plotting
        - It basically wraps the matplotlib 


    Input:
        x_list = The list of x inputs, it can be 1D or 2D
        y_list = The list of y inputs, it can be 1D or 2D
        labels = The list of labels for each plotting set
        xlabel = X label of the whole 1D set
        ylabel = Y label of the whole 1D set
        title  = Title of the whole 1D set
        savefn = If it is specified, it saves the figure to the specified location 

    Output:
        []

    Example:
        []
    """
    # The list of markers to be used for each plot
    markers = ["o", "^", "h", "s", ".", "p"]
    # Plot counter
    cnt = 0

    if labels is not None:
        # If there is individual X's for each Y
        if len(x_list) == len(y_list):
            for (x_, y_) in zip(x_list, y_list):
                plt.scatter(x_, y_, label=labels[cnt], marker=markers[cnt])
                cnt += 1
        else:
            for y_ in y_list:
                plt.scatter(x_list[0], y_, label=labels[cnt],
                            marker=markers[cnt])
                cnt += 1

    # If labels are not provided, then do not take into account
    else:
        # If there is individual X's for each Y
        if len(x_list) == len(y_list):
            for (x_, y_) in zip(x_list, y_list):
                plt.scatter(x_, y_, marker=markers[cnt])
                cnt += 1
        else:
            for y_ in y_list:
                plt.scatter(x_list, y_, marker=markers[cnt])
                cnt += 1

    if xlabel is not None:
        plt.xlabel(xlabel)
    if ylabel is not None:
        plt.ylabel(ylabel)
    if title is not None:
        plt.title(title)
    if labels is not None:
        plt.legend(loc="best")
    plt.tight_layout()

'''

# Regular figure initialization and termination


def figure_(*p): plt.figure(*p)
def show_(): plt.show()
