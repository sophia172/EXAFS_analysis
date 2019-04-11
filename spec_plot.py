from xafs.def_sample import *


experiment = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Experiment_Analysis/CdS/CdS_XAFS_CdK_chi_2018.txt')
k = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Experiment_Analysis/CdS/CdS_XAFS_CdK_chi_2018.txt')[:,0]



##############################################
#
#
#    Plot Artemis fit
#
#
##################################################
#
# plt.subplot(2,1,1)
# exafs_func.plot_spec(exafs_func,'M311_CdO.k',k,0,2)
# # plt.plot(Artemis_fit_OS[:,0], Artemis_fit_OS[:,1] * Artemis_fit_OS[:,0]**2,label='Cd-O,Cd-S')
# exafs_func.plot_spec(exafs_func,'M311_CdS.k',k,0,2)
# exafs_func.plot_spec(exafs_func,'M311_CdCd.k',k,0,2)
# exafs_func.plot_spec(exafs_func,'M311_CdOSCd.k',k,0,2)
# plt.plot(experiment[:,0], experiment[:,3] * experiment[:,0]**2,'k.',markersize=2,label='M311')
# plt.legend(loc='upper right')
# plt.xlabel('k')
# plt.ylabel('Chi')
#
# plt.subplot(2,1,2)
# exafs_func.plot_diff(exafs_func,'M311_CdO.k',k,experiment[:,3],0,2)
# exafs_func.plot_diff(exafs_func,'M311_CdS.k',k,experiment[:,3],0,2)
# exafs_func.plot_diff(exafs_func,'M311_CdCd.k',k,experiment[:,3],0,2)
# exafs_func.plot_diff(exafs_func,'M311_CdOSCd.k',k,experiment[:,3],0,2)
# plt.legend(loc='upper right')
# plt.xlabel('k')
# plt.ylabel('Chi')
#
# plt.show()

##############################################
#
#
#    Compare InP model (InP/InPO/InPOC) FEFF result with experiment
#
#
##################################################

# exafs_func.plot_xmus(exafs_func,'CdS40')

