####  fit Cd-O and Cd-S ratio
from EXAFS_analysis.xafs import *
# from ipywidgets import interact

def interpolation(f, x):
    g = interp1d(f[:, 0], f[:, 1])
    return g(x.tolist())

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

def XANES_ratio_fit(sample='M311'):
    CdO_data = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Experiment_Analysis/CdS_CdO_ref.nor')[:,(0,2)]
    CdS_data = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Experiment_Analysis/CdS/CdS_XAFS_CdK_mu_2018.xmu')[:,(0,2)]
    experiment = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Experiment_Analysis/CdS/CdS_XAFS_CdK_mu_2018.xmu')
    energy = experiment[:,0]
    if sample == 'M311':
        mu = experiment[:,(0,1)]
    if sample == 'bulk':
        mu = experiment[:,(0,2)]
    if sample == 'M322':
        mu = experiment[:,(0,3)]
    e = np.linspace(-10,50,400)
    e0 =CdS_data[:,0][np.argmax(CdS_data[:,1])]

    e = e + e0

    CdO = interpolation(CdO_data,e+3)
    CdS = interpolation(CdS_data,e)

    mu = interpolation(mu,e)
    plt.plot(e,mu,'k.',label=sample)
    plt.plot(e,CdO-0.2)
    plt.plot(e,CdS-0.2)
    plt.legend()
    plt.show()

    plt.figure()
    axes = plt.gca()

    def show_ratio(a,b,c1,c2):
        axes.clear()
        axes.plot(e, a*CdO+c1,label='CdO a={:f},c1={:f}'.format(a,c1))
        axes.plot(e, b*CdS+c2,label='Cds a={:f},c1={:f}'.format(a,c1))
        axes.plot(e, mu,'k.', label=sample)
        axes.plot(e, a*CdO+b*CdS+c1+c2, label='fit_result')
    interact(show_ratio,a=(0,1,0.01),b=(0,1,0.01),c1=(0.5,0.5,0.1),c2=(0.5,0.5,0.1))

    def g(e,a,b):
        return a*(CdO-0.2)/e*e+b*(CdS-0.2)/e*e
    popt, pcov = curve_fit(g,e,mu)
    print('ratio for CdO & CdS = ', *popt)
    plt.plot(e,mu,'k.',label=sample)
    plt.plot(e,g(e,*popt),label='fit result')
    plt.legend()
    plt.show()






