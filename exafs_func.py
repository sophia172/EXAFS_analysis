
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

class exafs_func:

    from scipy.optimize import curve_fit

    def __init__(self):
        return

    def interpolation(f,x):
        g = interp1d(f[:,0],f[:,1])
        return g(x.tolist())


    def plot_spec(self,filename,k,column1,column2):
        data = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Statistic Analysis/demeter/{:s}'.format(filename))[:,(column1,column2)]
        data[:,1] = self.interpolation(data,k)
        plt.plot(k,data[:,1] * k**2,label = '{:s}'.format(filename))

    def plot_diff(self,filename,k,y0,column1,column2,):
        data = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Statistic Analysis/demeter/{:s}'.format(filename))[:,(column1,column2)]
        data[:,1] = self.interpolation(data,k)
        plt.plot(k,data[:,1] * k**2 - y0 * k**2,label = '{:s}'.format(filename))
