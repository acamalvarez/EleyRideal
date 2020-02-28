import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def set_style(fontSize=14, linewidth=1):
    """ 
    Set the geneal style for the graphs. Do not require entries. 
    The parameters of this function should be change if the graphs do not look good.
    Example:
    For a bigger font size increase plt.rcParams['font.size'] 
    """
#    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Verdana"
    plt.rcParams['font.size'] = fontSize
    plt.rcParams['font.weight'] = 'regular'
    plt.rcParams['mathtext.default'] = 'regular'
    plt.rcParams['savefig.dpi'] = '500'
    plt.rcParams['savefig.transparent'] = True
    plt.rcParams['lines.linewidth'] = linewidth
    plt.rcParams['lines.markersize'] = 6
    plt.rcParams['legend.fontsize'] = fontSize - 2
    plt.rcParams['legend.fancybox'] = False
    plt.rcParams['legend.labelspacing'] = 0.5
    plt.rcParams['legend.columnspacing'] = 1
    plt.rcParams['legend.borderpad'] = 0.5
    plt.rcParams['axes.labelweight'] = 'bold'
    #plt.rcParams['figure.autolayout'] = 'True'

    # Pandas options
    pd.options.display.max_columns = 100
    pd.options.display.max_rows = 200

##############################################################################


def set_labels(x, y=''):
    """ Set the labels of the graph

    x -- set the label of x axis
    y -- set the label of y axis
    """
    plt.xlabel(x)
    plt.ylabel(y)


def save_fig(name, folder='', show=False):
    """
    Save the fig in the folder Graphs.
    """
    plt.tight_layout()
    path = 'Graphs/' + folder + '/' + name
    plt.savefig(path, bbox_inches='tight')
    if show == True:
        plt.show()
    plt.close()


def select_data(data, list):
    return [data[i] for i in list]


def kineticData():
    kineticData = pd.read_csv('kineticData.csv')

    kineticData['mass'] = np.where(
        kineticData['catalyst'] == 'Cu12AEZ', 1.5e-3, 3.0e-3)
    kineticData['flow'] = 450
    kineticData['MC'] = np.where(
        kineticData['catalyst'] == 'Cu12AEZ', 1.2, 1.4)
    kineticData['dispersion'] = 1

    kineticData['ML'] = kineticData['MC'] * \
        kineticData['mass'] / 100

    kineticData['TMS'] = kineticData['mass'] * \
        kineticData['MC'] * kineticData['dispersion'] / 100 / 63.55

    kineticData['N_lost_from_NOx'] = kineticData['NO_in'] + \
        kineticData['NO2_in'] - kineticData['NO_out'] - kineticData['NO2_out']

    kineticData['N_lost_from_NH3'] = kineticData['NH3_in'] - \
        kineticData['NH3_out']

    kineticData['NO_converted'] = kineticData['NO_in'] - kineticData['NO_out']

    kineticData['NH3_converted'] = kineticData['NH3_in'] - \
        kineticData['NH3_out']

    kineticData['NO2_converted'] = kineticData['NO2_in'] - \
        kineticData['NO2_out']

    kineticData['N2_from_NO'] = kineticData['NO_converted'] - \
        kineticData['NO2_converted']

    kineticData['N2_from_NH3'] = kineticData['NH3_converted'] - \
        2 * kineticData['NO2_converted']

    kineticData['X_NO'] = 100 * \
        kineticData['NO_converted'] / kineticData['NO_in']

    kineticData['X_NH3'] = 100 * \
        kineticData['NH3_converted'] / kineticData['NH3_in']

    kineticData['S_N2O'] = 100 * (kineticData['N2O_out'] -
                                  kineticData['N2O_in']) / kineticData['N_lost_from_NOx']

    kineticData['N2_flow'] = kineticData['N2_from_NO'] * \
        kineticData['flow'] / 1e6

    kineticData['rateMeasured'] = kineticData['N2_flow'] * \
        101.325 / (60000 * 8.314 * 298)

    # kineticData['rateMeasured'] = rateMeasured / kineticData['ML']  # mol NO / s

    kineticData['TOR'] = kineticData['rateMeasured'] / kineticData['TMS']

    return kineticData


class HC_effect():
    HCs = ['None', 'C3', 'C12']
    labels = ['Without HC', 'With C$_3$H$_6$', 'With C$_{12}$H$_{26}$']
    symbols = ['o', 's', '^']

    def __init__(self, kineticData, compounds, points):

        self.kineticData = kineticData
        self.compounds = compounds
        self.points = points

    def plot(self, catalyst, x_axis, labels):

        ax = plt.figure().add_subplot(111)

        for i in np.arange(len(HC_effect.HCs)):

            x = self.kineticData[(self.kineticData['catalyst'] == catalyst) & (
                self.kineticData['compound'] == self.compounds[i]) & (self.kineticData['HC'] == HC_effect.HCs[i])][x_axis]

            x = select_data(np.array(x), self.points[i])

            y = self.kineticData[(self.kineticData['catalyst'] == catalyst) & (
                self.kineticData['compound'] == self.compounds[i]) & (self.kineticData['HC'] == HC_effect.HCs[i])]['TOR']

            y = select_data(np.array(y), self.points[i])

            y_pol = np.polyfit(np.log(x), np.log(y), 1)
            y_val = np.polyval(y_pol, np.log(x))

            ax.plot(np.log(x), np.log(y),
                    self.symbols[i], label=self.labels[i], fillstyle='none')
            ax.plot(np.log(x), y_val, '-')

        ax.legend()
        set_labels(labels[0], labels[1])
        plt.show()


def r_SCR(conc, k_SCR, K_a, k_ox, k_o):

    return k_SCR * conc[0] * conc[1] / (1 + K_a * conc[2]) + k_ox * K_a * conc[2] / (1 + K_a * conc[2]) + k_o * conc[1]


def r_SCR2(conc, k1, k2, k3,):

    p0 = conc[0] * 0.082 * (250 + 273) / (22.71108 * 1e6)
    p1 = conc[1] * 0.082 * (250 + 273) / (22.71108 * 1e6)
    p2 = conc[2] * 0.082 * (250 + 273) / (22.71108 * 1e6)

    den = k2 * k3 * p1 * p2 + 4 * k1 * k3 * \
        p0 * p2 + 4 * k1 * k2 * p0 * p1

    theta2 = 4 * k1 * k3 * p0 * p2 / den

    return k2 * theta2 * p1 / den
