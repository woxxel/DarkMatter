import matplotlib.pyplot as plt

def set_plot_params():
    # plot parameters
    plt.rc('text', usetex=True)
    plt.rc('font', family='sans-serif')
    # plt.rcParams['font.family'] = ['Tahoma','Verdana']
    plt.rcParams['font.size'] = 12        #### how to get proper and same fontsizes with no specified xticklabels and with specified ones?
    plt.rcParams['xtick.labelsize'] = 10
    plt.rcParams['ytick.labelsize'] = 10
