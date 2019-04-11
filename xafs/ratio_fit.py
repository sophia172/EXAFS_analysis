####  fit Cd-O and Cd-S ratio
from xafs.__init__ import *

def bond_ratio_fit(self, x, N_CdO, N_CdS):
    return self.EXAFS_model_OS(x, N_CdS, self.R_CdS, self.ss_CdS, N_CdO, self.R_CdO, self.ss_CdS)

def ratio_fit(self):
    popt, pcov = curve_fit(self.bond_ratio_fit, self.k.astype(float), (self.experiment * self.k ** 2).astype(float))
    print('N_CdO = ', popt[0])
    print('N_CdS = ', popt[1])
    plt.subplot(2,1,1)
    plt.plot(self.k, self.bond_ratio_fit(self.k, *popt), label='Cd-S/Cd-O bond ratio = {:f}'.format(popt[0]/popt[1]))
    plt.plot(self.k, self.experiment * self.k ** 2, 'k.', markersize=2, label=self.sample)
    plt.subplot(2,1,2)
    plt.plot(self.k, self.bond_ratio_fit(self.k, *popt) - self.experiment * self.k ** 2,
             label='Cd-S/Cd-O bond ratio = {:f}'.format(popt[0]/popt[1]))
    plt.figlegend()
    plt.xlabel('k')
    plt.ylabel('Chi')
    plt.title('Use artemis (bulk+CdO) crystal result to fit Cd-S/Cd-O ratio')
    plt.show()

