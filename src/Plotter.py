import matplotlib.pyplot as plt



def plot_(x_list, y_list, labels = None, xlabel = None, ylabel = None, title = None, savefn = None):
    """
    docstring
    """
    markers = ["o", "^", "h", "s", ".", "p"]
    cnt = 0

    if labels is not None:
        if len(x_list) == len(y_list):
            for (x_, y_) in zip(x_list, y_list):
                plt.plot(x_, y_, label = labels[cnt], marker = markers[cnt])
                cnt += 1
        else:
            for y_ in y_list:
                plt.plot(x_list, y_, label = labels[cnt], marker = markers[cnt])
                cnt += 1
    else:
        if len(x_list) == len(y_list):
            for (x_, y_) in zip(x_list, y_list):
                plt.plot(x_, y_, marker = markers[cnt])
                cnt += 1
        else:
            for y_ in y_list:
                plt.plot(x_list, y_, marker = markers[cnt])
                cnt += 1
    
    if xlabel is not None: plt.xlabel(xlabel)
    if ylabel is not None: plt.ylabel(ylabel)
    if title is not None: plt.title(title)
    
    plt.tight_layout()


def scatter_(x_list, y_list, labels = None, xlabel = None, ylabel = None, title = None, savefn = None):
    """
    docstring
    """
    markers = ["o", "^", "h", "s", ".", "p"]
    cnt = 0

    if labels is not None:
        if len(x_list) == len(y_list):
            for (x_, y_) in zip(x_list, y_list):
                plt.scatter(x_, y_, label = labels[cnt], marker = markers[cnt])
                cnt += 1
        else:
            for y_ in y_list:
                plt.scatter(x_list[0], y_, label = labels[cnt], marker = markers[cnt])
                cnt += 1
    else:
        if len(x_list) == len(y_list):
            for (x_, y_) in zip(x_list, y_list):
                plt.scatter(x_, y_, marker = markers[cnt])
                cnt += 1
        else:
            for y_ in y_list:
                plt.scatter(x_list, y_, marker = markers[cnt])
                cnt += 1
    
    if xlabel is not None: plt.xlabel(xlabel)
    if ylabel is not None: plt.ylabel(ylabel)
    if title is not None: plt.title(title)
    if labels is not None: plt.legend(loc = "best")
    plt.tight_layout()


def figure_(): plt.figure()
def show_(): plt.show()
