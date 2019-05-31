from xafs import *





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

##############################################
#
#
#    Compare pyspline and athena
#
#
##################################################

M311_pyspline_Andrei = exafs_func(sample='M311', k_min = 2, k_max =17,fitted_curve='pyspline_Andrei')
M311_pyspline_Ying = exafs_func(sample='M311', k_min = 2, k_max =17,fitted_curve='pyspline_Ying')
M311_athena = exafs_func(sample='M311', k_min = 2, k_max =17,fitted_curve='athena')

M322_pyspline_Andrei = exafs_func(sample='M322', k_min = 2, k_max =17,fitted_curve='pyspline_Andrei')
M322_pyspline_Ying = exafs_func(sample='M322', k_min = 2, k_max =17,fitted_curve='pyspline_Ying')
M322_athena = exafs_func(sample='M322', k_min = 2, k_max =17,fitted_curve='athena')

def plot_chi(sample):
    plt.plot(sample.k,sample.experiment,label= sample.sample +' '+ sample.fitted_curve)

def plot_r(sample):
    r,amp,real,img=calcfft(sample.k,sample.experiment)
    plt.plot(r,amp,label= sample.sample +' '+ sample.fitted_curve)

plt.figure()
plt.subplot(2,1,1)
plot_chi(M311_athena)
plot_chi(M311_pyspline_Andrei)
plot_chi(M311_pyspline_Ying)
plt.legend()
plt.xlabel('k ($\AA^{-1}$)')
plt.ylabel('$\chi$')
plt.title('Experiments after background substraction')
plt.tight_layout()


plt.subplot(2,1,2)
plot_r(M311_athena)
plot_r(M311_pyspline_Andrei)
plot_r(M311_pyspline_Ying)
plt.legend()
plt.xlim([0,5])
plt.xlabel('R ($\AA$)')
plt.ylabel('Intensity')
plt.title('FT of experiment')
plt.tight_layout()
plt.savefig('/Users/Sophia/ownCloud/PhD/Statistic Analysis/figure/M311_experiment.pdf',format='pdf')
plt.show()

plt.figure()
plt.subplot(2,1,1)
plot_chi(M322_athena)
plot_chi(M322_pyspline_Andrei)
plot_chi(M322_pyspline_Ying)
plt.legend()
plt.xlabel('k ($\AA^{-1}$)')
plt.ylabel('$\chi$')
plt.title('Experiments after background substraction')
plt.tight_layout()

plt.subplot(2,1,2)
plot_r(M322_athena)
plot_r(M322_pyspline_Andrei)
plot_r(M322_pyspline_Ying)
plt.legend()
plt.xlim([0,5])
plt.xlabel('R ($\AA$)')
plt.ylabel('Intensity')
plt.title('FT of experiment')
plt.tight_layout()
plt.savefig('/Users/Sophia/ownCloud/PhD/Statistic Analysis/figure/M322_experiment.pdf',format='pdf')
plt.show()


