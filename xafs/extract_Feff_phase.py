
from xafs import *

###########################MAIN#######################################

print('\n#########################################################')
print('#\tChoose Chi file to be processed:\t\t#')

path1 = ('/Users/Sophia/ownCloud/PhD/Experiment_Analysis/CdS/CdS_bulk_Aug18_athena_1sh.chiq')
print('#\t   path:  ',path1,'                          \t#')
print('#########################################################\n')
print('#\t      K               XAFS  \t#')
print('#########################################################\n')
data1 = genfromtxt(path1)[:,(0,1)]
# print(data1)


print('\n#########################################################')
print('#\tChoose Chi file to be processed:\t\t#')

path2 = ('/Users/Sophia/ownCloud/PhD/Experiment_Analysis/CdS/CdS_bulk_Aug18_athena_2sh.chiq')
print('#\t   path:  ',path2,'                          \t#')
print('#########################################################\n')
print('#\t      K               XAFS  \t#')
print('#########################################################\n')
data2 = genfromtxt(path2)[:,(0,1)]
# print(data2)
# e = float(input('Define E0 in this analysis E0 = '))
# R = float(input('Define R in this analysis R = '))
# N = float(input('Define N in this analysis N = '))

path3 = ('/Users/Sophia/ownCloud/PhD/Experiment_Analysis/CdS/CdO_ref_1sh.chiq')
print('#\t   path CdO:  ',path3,'                          \t#')
print('#########################################################\n')
print('#\t      K               XAFS  \t#')
print('#########################################################\n')
data3 = genfromtxt(path3)[:,(0,1)]
# print(data3)

e = 0
R1 = 2.5118
N1 = 4
N2 = 12
R2 =4.1016
N3 = 6
R3 = 2.3475


k1 = np.sqrt(data1[:,0]**2+0.2625*e)
chi1 = data1[:,1]
k_min = (np.abs(k1 - 2)).argmin()
k_max = (np.abs(k1 - 16)).argmin()
k_fit = k1[k_min:k_max]


k2 = np.sqrt(data2[:,0]**2+0.2625*e)
chi2 = data2[:,1]
chi1_fit = interpolation(np.vstack((k1,chi1)).transpose(),k_fit)/k_fit**2
chi2_fit = interpolation(np.vstack((k2,chi2)).transpose(),k_fit)/k_fit**2
k3 = np.sqrt(data3[:,0]**2+0.2625*e)
chi3 = data3[:,1]
# print(k3,k_fit)
chi3_fit = interpolation(np.vstack((k3,chi3)).transpose(),k_fit)/k_fit**2

####################################################################
#
#
#  Experiment data to be fitted
#
#
#
##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
def fitted_exp(experiment):
    if experiment == 'pyspline_M311':
        path = '/Users/Sophia/ownCloud/PhD/Experiment_Analysis/CdS/CdS_M311_Aug18_pyspline.dat'
        data = np.genfromtxt(path)[:,(0,1)]
        exp = data[:,(0,1)]
    elif experiment == 'pyspline_M322':
        path = '/Users/Sophia/ownCloud/PhD/Experiment_Analysis/CdS/CdS_M322_Aug18_pyspline.dat'
        data = np.genfromtxt(path)[:,(0,1)]
        exp = data[:,(0,1)]
    elif experiment == 'pyspline_bulk':
        path = '/Users/Sophia/ownCloud/PhD/Experiment_Analysis/CdS/CdS_bulk_Aug18_pyspline.dat'
        data = np.genfromtxt(path)[:,(0,1)]
        exp = data[:,(0,1)]
    elif experiment == 'athena_R':
        path = '/Users/Sophia/ownCloud/PhD/Experiment_Analysis/CdS/CdS_R_Nov17.chik'
        data = np.genfromtxt(path)[:,(0,1)]
        exp = data[:,(0,1)]
    elif experiment == 'athena_bulk_12sh':
        path = '/Users/Sophia/ownCloud/PhD/Experiment_Analysis/CdS/CdS_bulk_Aug18_athena_12sh.chiq'
        data = np.genfromtxt(path)[:,(0,1)]
        data[1:,1]=data[1:,1]/data[1:,0]**2
        exp = data[:,(0,1)]
    elif experiment == 'athena_M311':
        path = '/Users/Sophia/ownCloud/PhD/Experiment_Analysis/CdS/CdS_XAFS_CdK_chi_2018.txt'
        data = np.genfromtxt(path)[:,(0,2,3,4)]
        exp = data[:,(0,1)]
    elif experiment == 'athena_bulk':
        path = '/Users/Sophia/ownCloud/PhD/Experiment_Analysis/CdS/CdS_XAFS_CdK_chi_2018.txt'
        data = np.genfromtxt(path)[:,(0,2,3,4)]
        exp = data[:,(0,2)]
    elif experiment == 'athena_M322':
        path = '/Users/Sophia/ownCloud/PhD/Experiment_Analysis/CdS/CdS_XAFS_CdK_chi_2018.txt'
        data = np.genfromtxt(path)[:,(0,2,3,4)]
        exp = data[:,(0,3)]
    chi_fit = interpolation(exp,k_fit)

    # plt.figure()
    # plt.plot(k_fit,chi1_fit,label='Cd-S shell')
    # plt.plot(k_fit,chi2_fit,label='Cd-Cd shell')
    # plt.plot(k_fit,chi3_fit,label='Cd-O shell')
    # plt.plot(k_fit,chi_fit,label= experiment)
    # plt.legend()
    # plt.savefig('/Users/Sophia/ownCloud/PhD/Experiment_Analysis/CdS/chi_shell.pdf',format='pdf')
    # plt.show()
    return chi_fit

##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#
#
#  Finish  Experiment data to be fitted
#
#
####################################################################

####################################################################
#
#
#  Find lambda for the file
#
#
#
##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

lambda1 = genfromtxt('/Users/Sophia/ownCloud/PhD/Simulation/CdS/FEFF/EXAFS/Cd/F-43m/feff0001.dat')[:,(0,5)]
lambda1 = interpolation(lambda1,k_fit)
Feff_CdS = genfromtxt('/Users/Sophia/ownCloud/PhD/Simulation/CdS/FEFF/EXAFS/Cd/F-43m/feff0001.dat')[:,(0,2)]
Feff_CdS = interpolation(Feff_CdS,k_fit)

lambda2 = genfromtxt('/Users/Sophia/ownCloud/PhD/Simulation/CdS/FEFF/EXAFS/Cd/F-43m/feff0002.dat')[:,(0,5)]
lambda2 = interpolation(lambda2,k_fit)
Feff_CdCd = genfromtxt('/Users/Sophia/ownCloud/PhD/Simulation/CdS/FEFF/EXAFS/Cd/F-43m/feff0002.dat')[:,(0,2)]
Feff_CdCd = interpolation(Feff_CdCd,k_fit)

lambda3 = genfromtxt('/Users/Sophia/ownCloud/PhD/Simulation/CdS/FEFF/EXAFS/CdO/feff0001.dat')[:,(0,5)]
lambda3 = interpolation(lambda3,k_fit)
Feff_CdO = genfromtxt('//Users/Sophia/ownCloud/PhD/Simulation/CdS/FEFF/EXAFS/CdO/feff0001.dat')[:,(0,2)]
Feff_CdO = interpolation(Feff_CdO,k_fit)
##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#
#
#  Finish find lambda for the file
#
#
####################################################################


def extract_CdS():





    amp = np.abs(hilbert(chi1_fit))
    phase = np.unwrap(np.angle(hilbert(chi1_fit)))+np.pi/2 - 2*k_fit*R1
    Feff_ss_CdS = amp/np.exp(-2*R1/lambda1)/N1*k_fit*R1**2




    plt.figure()
    plt.plot(k_fit,Feff_ss_CdS,label='experiment CdS')
    plt.plot(k_fit,Feff_CdS,label='FEFF calculation')
    plt.legend()
    plt.show()


    return Feff_ss_CdS,phase


def extract_CdO():





    amp = np.abs(hilbert(chi3_fit))
    phase = np.unwrap(np.angle(hilbert(chi3_fit)))+np.pi/2 - 2*k_fit*R3
    Feff_ss_CdO = amp/np.exp(-2*R3/lambda3)/N3*k_fit*R3**2



    plt.figure()
    plt.plot(k_fit,Feff_ss_CdO,label='experiment CdO')
    plt.plot(k_fit,Feff_CdO,label='FEFF calculation')
    plt.legend()
    plt.show()


    return Feff_ss_CdO,phase


def extract_CdCd():
    amp = np.abs(hilbert(chi2_fit))
    phase = np.unwrap(np.angle(hilbert(chi2_fit)))+np.pi/2 - 2*k_fit*R2
    Feff_ss_CdCd = amp/np.exp(-2*R2/lambda2)/N2*k_fit*R2**2

    plt.figure()
    plt.plot(k_fit,Feff_ss_CdCd,label='experiment CdCd')
    plt.plot(k_fit,Feff_CdCd,label='FEFF calculation')
    plt.legend()
    plt.show()

    return Feff_ss_CdCd,phase



def fit_exp(chi_fit,seed):

    def model_shell(x,N,R,del_ss):
        if R > 3            :
            shell = 2
        elif R < 2.47:
            shell = 3
        else:
            shell = 1
        return N * x ** 2 \
               * np.abs( eval('Feff_ss'+str(shell))) / (x * R ** 2) \
               * np.sin(2 * x * R +  eval('phase'+str(shell))) \
               * np.exp(-2 * R /  eval('lambda'+str(shell))) \
               * np.exp(-2 * del_ss * x ** 2)

    bounds_dic = {'CdS': np.array([(0, 6), (2.47, 2.55), (-0.00001, 0.01)]),
                  'CdO': np.array([(0, 8), (2.2, 2.4699), (0.0001, 0.03)]),
                  'CdCd': np.array([(0,12), (3.5, 4.2), (-0.001, 0.03)]),
                  'Ge3': np.array([(11.999,12), (4.68, 4.70), (0.002, 0.05)]),
                  'R' : np.array([(2.449, 2.451)]),
                  'e': np.array([(0,20)]),
                  'S0': np.array([(0.6,1)]),
                  'a':np.array([(-2,2)]),
                  'b':np.array([(-20,20)])}
    # bounds = np.vstack((bounds_dic['Ge'+str(shell)], bounds_dic['S0'], bounds_dic['e']))
    bounds = np.vstack((bounds_dic['CdO'],bounds_dic['CdS'],bounds_dic['CdCd'], bounds_dic['S0']))

    def LLE_model(params):
        S0 = params[-1]
        # print('ss',params[-3])
        # print('k',k)
        # print('k_new',k_new)
        y = chi_fit/S0
        model_LLE = model_shell(k_fit, *params[:3]) + model_shell(k_fit, *params[3:6]) + model_shell(k_fit, *params[6:9])
        LL = -np.sum(stats.norm.logpdf(y * k_fit**2, loc=model_LLE))
        return LL


    LLE_fitresult = differential_evolution(LLE_model,bounds=bounds,strategy='currenttobest1bin',maxiter=100000,seed = seed)
    # LLE_fitresult = shgo(LLE_model,bounds=bounds)

    S0,del_ss1,del_ss2=LLE_fitresult.x[-1],LLE_fitresult.x[2],LLE_fitresult.x[-2]
    # print(S0,del_ss1,del_ss2)
    params = LLE_fitresult.x
    # print('N1    R1   SS1    N2   R2   SS2   S0')
    print(params)

    model_LLE = model_shell(k_fit, *params[:3]) + model_shell(k_fit, *params[3:6]) + model_shell(k_fit, *params[6:9])
    LL = LLE_model(params)

    return model_LLE,LL,params[:-1]

# choose from pyspline_M311,pyspline_M322,pyspline_bulk,athena_R,athena_bulk_12sh,athena_M311,athena_bulk,athena_M322
parameter ={}

def gen_fit_result(experiment,seed):
    parameter[experiment+'_chi_fit'] = fitted_exp(experiment)
    parameter[experiment+'_model_LLE'],parameter[experiment+'_LL'],parameter[experiment+'_params'] = fit_exp(parameter[experiment+'_chi_fit'],seed)

def plot_fit(experiment):
    plt.plot(k_fit,parameter[experiment+'_chi_fit']*k_fit**2,'.',markersize = 0.5,label=experiment)
    plt.plot(k_fit,parameter[experiment+'_model_LLE'],linewidth = 0.5,label=experiment+'_fit')
    plt.plot(k_fit,parameter[experiment+'_model_LLE']-parameter[experiment+'_chi_fit']*k_fit**2,linewidth = 0.5,label=experiment+'_difference')

Feff_ss1,phase1 = extract_CdS()
Feff_ss2,phase2 = extract_CdCd()
Feff_ss3,phase3 = extract_CdO()


def loop_min_LL(experiment,loop_num):
    gen_fit_result(experiment,seed = None)
    prev_chi_fit = parameter[experiment+'_chi_fit']
    prev_model_LLE = parameter[experiment+'_model_LLE']
    prev_LL = parameter[experiment+'_LL']
    prev_params = parameter[experiment+'_params']

    for i in range(loop_num):
        gen_fit_result(experiment,seed = i)
        if parameter[experiment+'_LL'] < prev_LL:
            print('keep current run, LL = ',parameter[experiment+'_LL'], 'previous_LL = ',prev_LL)
            prev_chi_fit = parameter[experiment+'_chi_fit']
            prev_model_LLE = parameter[experiment+'_model_LLE']
            prev_LL = parameter[experiment+'_LL']
            prev_params = parameter[experiment+'_params']

        else:
            print('keep previous run, LL = ',prev_LL)
            continue
    parameter[experiment+'_chi_fit'] = prev_chi_fit
    parameter[experiment+'_model_LLE'] = prev_model_LLE
    parameter[experiment+'_LL'] = prev_LL
    prev_params = parameter[experiment+'_params']
    return


# gen_fit_result('athena_bulk_12sh')
# gen_fit_result('athena_R')
# gen_fit_result('athena_bulk')
# gen_fit_result('athena_M311')
# gen_fit_result('athena_M322')

experiment = 'athena_M322'
loop_min_LL(experiment,200)

plt.figure()
# plot_fit('athena_bulk_12sh')
# plot_fit('athena_R')
# plot_fit('athena_bulk')
plot_fit(experiment)
# plot_fit('athena_M322')
plt.legend(loc='lower right')
plt.xlabel('k ($\AA^{-1}$)')
plt.ylabel('$k^2\chi$')
print(*parameter[experiment+'_params'][:3])
plt.text(10,0.3,'LL = '+ str(parameter[experiment+'_LL'])+'\n'
         + '            N        R      del_ss\n'
         + 'Cd-O   {:.2f}  {:.2f}  {:.5f}\n'.format(*parameter[experiment+'_params'][:3])
         + 'Cd-S   {:.2f}  {:.2f}  {:.5f}\n'.format(*parameter[experiment+'_params'][3:6])
         + 'Cd-Cd {:.2f}  {:.2f}  {:.5f}\n'.format(*parameter[experiment+'_params'][6:9]))
plt.savefig('/Users/Sophia/ownCloud/PhD/Statistic Analysis/figure/fit_'+experiment[-4:]+'_(extracted_feature).pdf',format='pdf')
plt.show()
