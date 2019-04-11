####  fit Cd-O and Cd-S ratio
from xafs import *


def ratio_fit(sample):
    def bond_ratio_fit(x, N_CdO, N_CdS):
        return sample.EXAFS_model_OS(x, N_CdS, sample.R_CdS, sample.ss_CdS, N_CdO, sample.R_CdO, sample.ss_CdS)
    popt, pcov = curve_fit(bond_ratio_fit, sample.k.astype(float), (sample.experiment * sample.k ** 2).astype(float))
    print('N_CdO = ', popt[0])
    print('N_CdS = ', popt[1])

    plt.subplot(2,1,1)
    plt.plot(sample.k, bond_ratio_fit(sample.k, *popt), label='Cd-S/Cd-O bond ratio = {:f}'.format(popt[0]/popt[1]))
    plt.plot(sample.k, sample.experiment * sample.k ** 2, 'k.', markersize=2, label=sample.sample)
    plt.xlabel('k')
    plt.ylabel('Chi')
    plt.title('Use artemis (bulk+CdO) crystal result to fit Cd-S/Cd-O ratio')

    plt.subplot(2,1,2)
    plt.plot(sample.k, bond_ratio_fit(sample.k, *popt) - sample.experiment * sample.k ** 2,
             label='Cd-S/Cd-O bond ratio = {:f}'.format(popt[0]/popt[1]))
    plt.figlegend()
    plt.xlabel('k')
    plt.ylabel('Chi difference')
    plt.show()

