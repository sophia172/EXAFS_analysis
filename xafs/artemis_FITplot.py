from EXAFS_analysis.xafs import *

def artemis_fit(sample):
    chi_list = glob.glob('/Users/Sophia/ownCloud/PhD/Statistic Analysis/data/{:s}_Aug18_{:s}*.k'.format(sample.sample,sample.fitted_curve))
    chik2_list = glob.glob('/Users/Sophia/ownCloud/PhD/Statistic Analysis/data/{:s}_Aug18_{:s}*.k2'.format(sample.sample,sample.fitted_curve))

    plt.figure()
    plt.subplot(2,1,1)
    all_fits = []
    fitted_curve1 = np.genfromtxt(chi_list[0])[:,(0,1)]
    fitted_curve2 = np.genfromtxt(chik2_list[0])[:,(0,1)]
    k = sample.k
    chi = sample.interpolation(fitted_curve1,sample.k)
    chik2 = sample.interpolation(fitted_curve2,sample.k)
    for i in chi_list:
        fit = np.genfromtxt(i)[:,(0,2)]
        fit = sample.interpolation(fit,sample.k)
        all_fits.append(fit)
        plt.plot(k,fit,linewidth=0.5,label=i[-8:-2])
    plt.plot(k,chi,'k.',markersize=0.5,label='{:s} {:s}'.format(sample.sample, sample.fitted_curve))
    plt.title('chi fit in Artemis {:s} {:s}'.format(sample.sample, sample.fitted_curve))
    plt.legend(loc='center right')

    plt.subplot(2,1,2)
    for i in all_fits:
        plt.plot(k,np.array(i)-chi,linewidth=0.5)
    plt.title('fit difference in Artemis {:s} {:s}'.format(sample.sample, sample.fitted_curve))
    plt.savefig('/Users/Sophia/ownCloud/PhD/Statistic Analysis/figure/artemis_fit_{:s}_{:s}_chi.pdf'.format(sample.sample,sample.fitted_curve),
                format='pdf')
    plt.show()

    plt.figure()
    plt.subplot(2,1,1)
    all_fits = []
    for i in chik2_list:
        fit = np.genfromtxt(i)[:,(0,2)]
        fit = sample.interpolation(fit,sample.k)
        all_fits.append(fit)
        plt.plot(k,fit,linewidth=0.5,label=i[-9:-3])
    plt.plot(k,chik2,'k.',markersize=0.5,label='{:s} {:s}'.format(sample.sample, sample.fitted_curve))
    plt.title('chik^2 fit in Artemis {:s} {:s}'.format(sample.sample, sample.fitted_curve))
    plt.legend(loc='center right')

    plt.subplot(2,1,2)
    for i in all_fits:
        plt.plot(k,np.array(i)-chik2,linewidth=0.5)
    plt.title('fit difference in Artemis {:s} {:s}'.format(sample.sample, sample.fitted_curve))
    plt.savefig('/Users/Sophia/ownCloud/PhD/Statistic Analysis/figure/artemis_fit_{:s}_{:s}_chik2.pdf'.format(sample.sample,sample.fitted_curve),
                format='pdf')
    plt.show()








