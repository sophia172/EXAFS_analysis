import numpy as np
from scipy.interpolate import interp1d
from xafs import *
from scipy.optimize import curve_fit

def interpolation(f, x):
    g = interp1d(f[:, 0], f[:, 1])
    return g(x.tolist())


################################################
#
##
#  Test if Ge each peak adds up together
#
#
####################################################

chiq1 = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Experiment_Analysis/Ge/c_ge12_chi1.chiq')
chiq2 = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Experiment_Analysis/Ge/c_ge12_chi2.chiq')
chiq3 = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Experiment_Analysis/Ge/c_ge12_chi3.chiq')
#
# plt.figure()
# plt.plot(k,exp*k**2)
# plt.plot(chiq1[:,0],chiq1[:,1]+chiq2[:,1]+chiq3[:,1])
# plt.xlim([0,20])
# plt.show()

experiments = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Experiment_Analysis/Ge/c_ge12_chi.dat',skip_header=5)
k = experiments[:,0]
weight=2

k_min = (np.abs(k - 3)).argmin()
k_max = (np.abs(k - 20)).argmin()
k_fit = k[k_min:k_max]

# exp = exp[k_min:k_max]

FEFF_1 = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Simulation/Ge/feff0001.dat',
                         skip_header=15)
FEFF_2 = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Simulation/Ge/feff0002.dat',
                       skip_header=15)
FEFF_3 = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Simulation/Ge/feff0005.dat',
                       skip_header=15)

# data_excurve = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Experiment_Analysis/Ge/c_ge12_pot_pha.dat',skip_header=5)
# Feff1=data_excurve[:,(0,1)]
# Feff2=Feff1
# Feff3=Feff1
#
# phase1=data_excurve[:,(0,2)]
# phase1[:,1] += data_excurve[:,3]
# phase2=phase1
# phase3 = phase1



Feff1 = FEFF_1[:, (0, 2)]
Feff2 = FEFF_2[:, (0, 2)]
Feff3 = FEFF_3[:, (0, 2)]



phase1 = FEFF_1[:, (0, 1)]
         #+0.9
phase1[:, 1] += FEFF_1[:, 3]

phase2 = FEFF_2[:, (0, 1)]
         #+1.468
phase2[:, 1] += FEFF_2[:, 3]

phase3 = FEFF_3[:, (0, 1)]
         #+1.727
phase3[:, 1] += FEFF_3[:, 3]


lambda1 = FEFF_1[:, (0, 5)]
lambda2 = FEFF_2[:, (0, 5)]
lambda3 = FEFF_3[:, (0, 5)]

Feff1 =  interpolation(Feff1,  k_fit)
Feff2 =  interpolation(Feff2,  k_fit)
Feff3 =  interpolation(Feff3,  k_fit)
lambda1 =  interpolation(lambda1,  k_fit)
lambda2 =  interpolation(lambda2,  k_fit)
lambda3 =  interpolation(lambda3,  k_fit)
# phase1 =  interpolation(phase1,  k_fit)
phase2 =  interpolation(phase2,  k_fit)
phase3 =  interpolation(phase3,  k_fit)

###############################################################
##
#
#  Calculate phase from experimental data
#
#  1. Use e0 alligned Chi result
#  2. Withdraw knot position and calculate phase
#
###################################################################
def calc_phase(shell):
    data1 = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Experiment_Analysis/Ge/ge12_1sh_ev.dat',skip_header=2)
    data2 = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Experiment_Analysis/Ge/ge12_2sh_ev.dat',skip_header=2)
    data3 = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Experiment_Analysis/Ge/ge12_3sh_ev.dat',skip_header=2)
    k_fine = np.linspace(3,20,3401)
    # interpolate each file with k_fine
    k1 = np.sqrt((data1[:,0]+13.8)/3.81)
    chi1= data1[:,1]
    exp = np.vstack((k1,chi1)).transpose()
    chi1 = interpolation(exp,k_fine)

    k2 = np.sqrt((data2[:,0]+13.8)/3.81)
    chi2= data2[:,1]
    exp = np.vstack((k2,chi2)).transpose()
    chi2 = interpolation(exp,k_fine)

    k3 = np.sqrt((data3[:,0]+13.8)/3.81)
    chi3= data3[:,1]
    exp = np.vstack((k3,chi3)).transpose()
    chi3 = interpolation(exp,k_fine)

    if shell == 1:
        chi = chi1
        R = 2.45
    elif shell == 2:
        chi = chi2-chi1
        R = 4.0
    elif shell == 3:
        chi = chi3-chi2
        R = 4.69
        correction = 2* pi *10



    peak_index = find_peaks(chi)[0]
    valley_index = find_peaks(-1*chi)[0]
    zero_index = np.where(np.diff(np.sign(chi)))[0]

    # find k, chi and phase for each -1,1,0 knot point
    chi_peak = np.array([])
    chi_valley = np.array([])
    chi_zero = np.array([])
    k_fine_peak = np.array([])
    k_fine_valley = np.array([])
    k_fine_zero = np.array([])
    # for i in peak_index:
    #     chi_peak = np.append(chi_peak,chi[i])
    #     k_fine_peak = np.append(k_fine_peak,k_fine[i])
    #
    # for i in valley_index:
    #     chi_valley = np.append(chi_valley,chi[i])
    #     k_fine_valley = np.append(k_fine_valley,k_fine[i])
    # for i in zero_index:
    #     chi_zero = np.append(chi_zero,chi[i])
    #     k_fine_zero= np.append(k_fine_zero,k_fine[i])

    # Figure out where the curve start
    index = np.hstack((peak_index,valley_index,zero_index))
    index = np.sort(index)
    # print(index)
    phase_knot = np.array([])
    phase_peak = np.array([])
    if index[0] in peak_index:
        print('curve start with peak',index)
        for n,i in enumerate(index):
            phase_knot = np.append(phase_knot, (pi/2 *(n+1)+ correction - 2 * k_fine[i] *R))
            print('2kR: ',2 * k_fine[i] *R)
    elif index[0] in valley_index:
        print('curve start with valley')
        for n,i in enumerate(index):
            phase_knot = np.append(phase_knot, pi/2 *(n+3) - 2 * k_fine[i] *R + correction)
    elif index[1] in peak_index:
        print('curve start with 0 then peak')
        for n,i in enumerate(index):
            phase_knot = np.append(phase_knot, pi/2 *n - 2 * k_fine[i] *R + correction)
    elif index[1] in valley_index:
        print('curve start with 0 then valley')
        for n,i in enumerate(index):
            phase_knot = np.append(phase_knot, pi/2 *(n+2) - 2 * k_fine[i] *R + correction)
    k_knot = k_fine[index]
    for i in peak_index:
        phase_peak = np.append(phase_peak,phase_knot(np.where(index == i)))
    k_valley = k_fine[valley_index]

    # phase_peak = pi/2-2*k_fine_peak*R
    # phase_valley = pi*3/2-2*k_fine_valley*R
    # phase_zero = np.array([])
    # for n,i in enumerate(k_fine_zero):
    #     if n%2 == 0:
    #         phase_zero = np.append(phase_zero,0-2*i*R)
    #     elif n%2 == 1:
    #         phase_zero = np.append(phase_zero,0-2*i*R)
    #
    # # arrange the data in right sequence
    # index = np.hstack((peak_index,valley_index,zero_index))
    # k_knot = np.hstack((k_fine_peak,k_fine_valley,k_fine_zero))
    # chi_knot = np.hstack((chi_peak,chi_valley,chi_zero))
    # phase_knot = np.hstack((phase_peak,phase_valley,phase_zero))
    # all_knot=np.vstack((index,k_knot,chi_knot,phase_knot)).transpose()
    # all_knot = all_knot[all_knot[:,0].argsort()]
    #
    # for n,i in enumerate(phase_peak):
    #
    #     i = i + 2*pi*n
    #     phase_peak[n]=i
    #
    # for n,i in enumerate(phase_valley):
    #     i = i + 2*pi*n
    #     phase_valley[n]=i
    # for n,i in enumerate(phase_zero):
    #     i = i + pi*n
    #     phase_zero[n]=i



    def model_phase(x,a,b,c,d):
        return a*x**3 + b*x**2 + c*x + d

    popt, pcov = curve_fit(model_phase,k_knot,phase_knot)
    global phase
    phase = model_phase(k_fit,*popt[:4])
    print(*popt)
    # plt.plot(all_knot[:,1],all_knot[:,3])
    # plt.plot(k_fine,chi*100-20)
    plt.plot(k_knot,phase_knot,'.',label='knot')
    # plt.plot(k_fine_peak,phase_peak,'.',label='peak')
    # plt.plot(k_fine_valley,phase_valley,'.',label='valley')
    # plt.plot(k_fine_zero,phase_zero,'.',label='zero')
    plt.plot(k_fit,phase,label='phase')
    # plt.xlim([17,19])
    plt.legend()
    plt.show()
    return k_fine,chi,phase


# phase2 = calc_phase()
# phase3 = calc_phase()
###############################################################
##
#
#  Withdraw FEFF, lambda and phase from FEFF output file
#
#
###################################################################
def save_FEFF_calc():
    Ge_Feff_phase = np.vstack((k_fit,Feff1,lambda1,phase1,Feff2,lambda2,phase2,Feff3,lambda3,phase3)).transpose()
    np.savetxt('/Users/Sophia/ownCloud/PhD/Simulation/Ge/Ge_Feff_phase.dat',Ge_Feff_phase,header='# k_fit,Feff1,lambda1,phase1,Feff2,lambda2,phase2,Feff3,lambda3,phase3')
    plt.figure(figsize=(8,6))

    plt.subplot(3,1,1)
    plt.plot(k_fit,Feff1,label='Feff1')
    plt.plot(k_fit,Feff2,label='Feff2')
    plt.plot(k_fit,Feff3,label='Feff3')
    plt.xlabel('K $(\AA^{-1})$')
    plt.legend()
    plt.tight_layout()

    plt.subplot(3,1,2)
    plt.plot(k_fit,lambda1,label='lambda1')
    plt.plot(k_fit,lambda2,label='lambda2')
    plt.plot(k_fit,lambda3,label='lambda3')
    plt.xlabel('K $(\AA^{-1})$')
    plt.legend()
    plt.tight_layout()

    plt.subplot(3,1,3)
    plt.plot(k_fit,phase1,label='phase1')
    plt.plot(k_fit,phase2,label='phase2')
    plt.plot(k_fit,phase3,label='phase3')
    plt.xlabel('K $(\AA^{-1})$')
    plt.legend()
    plt.tight_layout()

    plt.savefig('/Users/Sophia/ownCloud/PhD/Simulation/Ge/Ge_Feff_phase.pdf',format='pdf')

# save_FEFF_calc()
#
# plt.plot(k,phase3)
# plt.show()
def EXAFS_model(x,N1,R1,ss1,N2,R2,ss2,N3,R3,ss3):
    print(N1,R1,ss1,N2,R2,ss2,N3,R3,ss3)
    return N1 * x **weight\
           * abs( Feff1.astype(float)) / (x.astype(float) * R1 ** 2) \
           * np.sin(2 * x.astype(float) * R1 +  phase1.astype(float)) \
           * np.exp(-2 * R1 /  lambda1.astype(float)) \
           * np.exp(-2 * ss1 * x.astype(float) ** 2) \
           + N2 * x **weight\
           * abs( Feff2.astype(float)) / (x.astype(float) * R2 ** 2) \
           * np.sin(2 * x.astype(float) * R2 +  phase2.astype(float)) \
           * np.exp(-2 * R2 /  lambda2.astype(float)) \
           * np.exp(-2 * ss2 * x.astype(float) ** 2)\
           + N3 * x **weight\
           * abs( Feff3.astype(float)) / (x.astype(float) * R3 ** 2) \
           * np.sin(2 * x.astype(float) * R3 +  phase3.astype(float)) \
           * np.exp(-2 * R3 /  lambda3.astype(float)) \
           * np.exp(-2 * ss3 * x.astype(float) ** 2)
#
# def EXAFS_shell(x,N,R,ss):
#     return N * x **weight \
#            * abs( Feff1.astype(float)) / (x.astype(float) * R ** 2) \
#            * np.sin(2 * x.astype(float) * R +  phase1.astype(float)) \
#            * np.exp(-2 * R /  lambda1.astype(float)) \
#            * np.exp(-2 * ss * x.astype(float) ** 2)
#
# plt.figure()
# plt.plot(k,EXAFS_shell(k,3.324,2.447,0.00236),label='path')
# plt.plot(chiq1[:,0],chiq1[:,1],label='RFT')
# plt.legend()
# plt.show()


params_dic = {'Ge1': [6, 2.35, 0.00234], 'Ge2': [12, 4, 0.005], 'Ge3': [12, 4.8, 0.006], 'e':[3], 'sd': [1]}
bounds_dic = {'Ge1': np.array([(3.999, 4), (2.44, 2.46), (0.0001, 0.005)]),
              'Ge2': np.array([(11.999,12), (3.99, 4.01), (0.002, 0.01)]),
              'Ge3': np.array([(11.999,12), (4.68, 4.70), (0.002, 0.01)]),
              'e': np.array([(0,30)]),
              'S0': np.array([(0.8,1)]),
              'a':np.array([(-2,2)]),
              'b':np.array([(-20,20)])}


params = params_dic['Ge1'] + params_dic['Ge2'] + params_dic['Ge3'] + params_dic['e']
bounds = np.vstack((bounds_dic['Ge1'], bounds_dic['Ge2'], bounds_dic['Ge3'], bounds_dic['S0'], bounds_dic['e']))
# bounds_NL = np.vstack((bounds_dic['Ge1'], bounds_dic['Ge2'], bounds_dic['Ge3'], bounds_dic['e']))
# popt, pcov = curve_fit(EXAFS_model, k.astype(float),(exp * k **2).astype(float), p0=params,bounds=bounds_NL.transpose())
# fit = EXAFS_model(k, *popt)
# print(*popt)

# def LLE_model(params):
#
#     e = params[-1]
#
#     exp = experiments[:,(0,1)]
#     model = EXAFS_model(k_fit, *params[:-1])
#     k_new=np.sqrt(k_fit**2-0.2625*e)
#     exp = interpolation(exp,k_new)
#     LL = -np.sum(stats.norm.logpdf(exp * k_fit**weight, loc=model))
#     return LL
#
#
#
# # Local minimum
# # LLE_fitresult = minimize(LLE_model, params, method='TNC', bounds=bounds,options={'maxiter': 12000})
# # Global minimum
# # LLE_fitresult = basinhopping(LLE_model, params,niter=200, take_step=mytakestep)
# # Global minimu rate 4/5
# LLE_fitresult = differential_evolution(LLE_model,bounds=bounds,strategy='randtobest1exp', popsize=100,maxiter=10000)
#
# # LLE_fitresult = shgo(LLE_model,bounds=bounds,iters=5)
#
# # LLE_fitresult = dual_annealing(LLE_model,bounds=bounds)
# fit_LLE = EXAFS_model(k_fit, *LLE_fitresult.x[:-1])
# print(*LLE_fitresult.x)
####################################
#
#       Import fit result from Artemis
#           and plot
#
######################################################
#
# chi_data = np.genfromtxt("/Users/Sophia/ownCloud/PhD/Simulation/Ge/c_ge12_chi.k2")[:,(0,2)]
#
#
#
# e=LLE_fitresult.x[-1]
# exp = experiments[:,(0,1)]
# k_new=np.sqrt(k_fit**2-0.2625*e)
# exp = interpolation(exp,k_new)
# e_artemis = 12.16917
# k_artemis = np.sqrt(k_fit**2-0.2625*e_artemis)
# chi_data = interpolation(chi_data, k_artemis)
#
#
# plt.figure()
# plt.subplot(2,1,1)
# plt.plot(k_fit, exp *k_fit**weight,'k', linewidth = 1, label = 'experiment')
# # plt.plot(k, fit, linewidth=1, label='nonlinear curve fit (scipy) of Germanium')
# plt.plot(k_fit,fit_LLE,linewidth=1, label='Python fit of Ge')
# # plt.plot(k, fit, linewidth=1, label='nonlinear curve fit (scipy) of Germanium')
# plt.plot(k_fit,chi_data,linewidth=1, label='Artemis fit of Ge')
# plt.xlim([2,18])
# plt.ylim([-3,3])
# plt.xlabel('k ($\AA^{-1}$)')
# plt.ylabel('$k^2\chi$')
# plt.legend()
# plt.tight_layout()
# data = np.vstack((k_fit,exp*k_fit**weight,fit_LLE,chi_data)).transpose()
# np.savetxt('/Users/Sophia/ownCloud/PhD/Simulation/Ge/Ge_python_artemis_fit.dat',data)
#
# plt.subplot(2,1,2)
# r,amp,real,imag = calcfft(k_fit,exp*k_fit**3,kmin=3,kmax=18)
# plt.plot(r,amp,'k', linewidth = 1, label = 'experiment')
# # r,amp,real,imag = calcfft(k,fit*k,kmin=3,kmax=18)
# # plt.plot(r,amp, linewidth=1, label='nonlinear curve fit (scipy) of Germanium')
# r,amp,real,imag = calcfft(k_fit,fit_LLE*k_fit,kmin=3,kmax=18)
# plt.plot(r,amp,linewidth=1, label='Python fit of Ge')
#
# r,amp,real,imag = calcfft(k_fit,chi_data*k_fit,kmin=3,kmax=18)
# plt.plot(r,amp,linewidth=1, label='Artemis fit of Ge')
# plt.xlabel('R ($\AA$)')
# plt.ylabel('Intensity')
# plt.legend()
# plt.xlim([0,5])
# plt.tight_layout()
# plt.savefig('/Users/Sophia/ownCloud/PhD/Statistic Analysis/figure/Ge_analysis.pdf',format='pdf')
# plt.show()

#####################################################
#
#      Import chi from each shell
#           and fit
#
######################################################



def Ge_shell():
    filelist = glob.glob('/Users/Sophia/ownCloud/PhD/Experiment_Analysis/Ge/ge12_*')
    print(filelist)
    global exp
    plt.figure(figsize=(10,8))
    for i,file in enumerate(filelist):
        data = np.genfromtxt(file,skip_header=2)
        k = np.sqrt(data[:,0]/3.81)
        chi = data[:,1]
        exp = np.vstack((k,chi)).transpose()
        # plt.plot(k,chi)
        # plt.show()


        def model_shell1(x,N,R,ss):
            return N * x **weight \
                   * abs( Feff1.astype(float)) / (x.astype(float) * R ** 2) \
                   * np.sin(2 * x.astype(float) * R +  phase1.astype(float)) \
                   * np.exp(-2 * R /  lambda1.astype(float)) \
                   * np.exp(-2 * ss * x.astype(float) ** 2)
        def model_shell2(x,N,R,ss):
            return N * x **weight \
                   * abs( Feff2.astype(float)) / (x.astype(float) * R ** 2) \
                   * np.sin(2 * x.astype(float) * R +  phase2.astype(float)) \
                   * np.exp(-2 * R /  lambda2.astype(float)) \
                   * np.exp(-2 * ss * x.astype(float) ** 2)
        def model_shell3(x,N,R,ss):
            return N * x **weight \
                   * abs( Feff3.astype(float)) / (x.astype(float) * R ** 2) \
                   * np.sin(2 * x.astype(float) * R +  phase3.astype(float)) \
                   * np.exp(-2 * R /  lambda3.astype(float)) \
                   * np.exp(-2 * ss * x.astype(float) ** 2)
        global model
        if i == 0:
            def model(x,N1,R1,ss1):
                return model_shell1(x,N1,R1,ss1)
            bounds = np.vstack((bounds_dic['Ge1'], bounds_dic['S0'], bounds_dic['e']))
        elif i == 1:
            def model(x,N1,R1,ss1,N2,R2,ss2):
                return model_shell1(x,N1,R1,ss1) + model_shell2(x,N2,R2,ss2)
            bounds = np.vstack((bounds_dic['Ge1'], bounds_dic['Ge2'], bounds_dic['S0'], bounds_dic['e']))
        elif i == 2:
            def model(x,N1,R1,ss1,N2,R2,ss2,N3,R3,ss3):
                return model_shell1(x,N1,R1,ss1) + model_shell2(x,N2,R2,ss2) + model_shell3(x,N3,R3,ss3)
            bounds = np.vstack((bounds_dic['Ge1'], bounds_dic['Ge2'], bounds_dic['Ge3'], bounds_dic['S0'], bounds_dic['e']))


        def LLE_model(params):
            global exp,model
            # print(exp)
            e = params[-1]
            S0 = params[-2]
            k_new = np.sqrt(k_fit**2-0.2625*e)
            exp_new = interpolation(exp,k_new)/S0
            model_LLE = model(k_fit, *params[:-2])
            LL = -np.sum(stats.norm.logpdf(exp_new * k_fit**weight, loc=model_LLE))
            return LL
        # print(bounds)
        # print(LLE_model((2,3,1,2)))

        LLE_fitresult = differential_evolution(LLE_model,bounds=bounds,strategy='randtobest1exp', popsize=200,maxiter=10000)
        fit_LLE = model(k_fit, *LLE_fitresult.x[:-2])
        LL = LLE_model(LLE_fitresult.x)
        e=LLE_fitresult.x[-1]
        S0=LLE_fitresult.x[-2]
        k_new=np.sqrt(k_fit**2-0.2625*e)
        exp_new = interpolation(exp,k_new)
        params = LLE_fitresult.x[:-2]
        params = [ '%.5f' % elem for elem in params]
        params = np.reshape(params,(-1,3))
        label_result = ''
        for a in range(i+1):
            label_result = label_result + str(params[a]) + '\n'
        plt.subplot(3,1,i+1)
        plt.plot(k_fit, exp_new*k_fit**weight,'k',linewidth=0.5)
        plt.plot(k_fit,fit_LLE*S0,linewidth=0.5,label='shell{:d}\n'.format(i+1)+'e= {:.2f}   S0 = {:.2f}  LL = {:.2F}\n\n'.format(e,S0,LL) + '     N           R          ss  \n'+ label_result )
        plt.plot(k_fit,fit_LLE*S0-exp_new*k_fit**weight -2,linewidth=0.5,label='difference')
        plt.xlim([2,26])
        plt.xlabel('K $(\AA^{-1})$')
        plt.ylabel('$\chi(k)k^2$')
        plt.legend()
        plt.tight_layout()
        data = np.vstack((k_fit,exp_new*k_fit**weight,fit_LLE*S0,exp_new*k_fit**weight-fit_LLE*S0)).transpose()
        np.savetxt('/Users/Sophia/ownCloud/PhD/Simulation/Ge/Ge_fit_shell{:d}.dat'.format(i+1),data,header='k  experiment_k^2 fit_k^2 difference')
    plt.savefig('/Users/Sophia/ownCloud/PhD/Statistic Analysis/figure/Ge_shellfit.pdf',format='pdf')
    plt.show()

# Ge_shell()

#####################################################
#
#      Fix N,R,ss,e, S0
#           and fit Feff and phase
#
######################################################


def calc_Feff_phase():


    def phase_model(b,c,d):
        return phase1 + (b*k_fit**2 + c*k_fit +d)
    def Feff_model(a,b,c,d):
        return Feff1 + a*k_fit**3 +b*k_fit**2 + c*k_fit +d

    def model_FEFF(x,b2,c2,d2):
        N = 4
        R = 2.45
        ss = 0.00305
        return N* x **weight \
               * np.abs(Feff1) / (x * R ** 2) \
               * np.sin(2 * x * R +  phase_model(b2,c2,d2)) \
               * np.exp(-2 * R /  lambda1) \
               * np.exp(-2 * ss * x ** 2)
    global exp
    data = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Experiment_Analysis/Ge/ge12_1sh_ev.dat',skip_header=2)
    k = np.sqrt(data[:,0]/3.81)
    chi = data[:,1]
    exp = np.vstack((k,chi)).transpose()



    def LLE_model_FEFF(params):
        global exp
        # print(exp)
        e =params[-1]
        S0 = 0.96
        k_new = np.sqrt(k_fit**2-0.2625*e)
        # params = params.tolist()
        print(params)
        exp_new = interpolation(exp,k_new)/S0

        model_LLE = model_FEFF(k_fit,*params[:-1])
        LL = -np.sum(stats.norm.logpdf(exp_new * k_fit**weight, loc=model_LLE))
        return LL

    bounds_FEFF =np.vstack((bounds_dic['a'],bounds_dic['a'],bounds_dic['b'],bounds_dic['e']))
    LLE_fitresult = differential_evolution(LLE_model_FEFF,bounds=bounds_FEFF,strategy='randtobest1exp', popsize=10,maxiter=10000)
    params = LLE_fitresult.x

    e =params[-1]
    print(e)
    fit_LLE = model_FEFF(k_fit, *params[:-1])
    LL = LLE_model_FEFF(params)
    params = [ '%.5f' % elem for elem in params]
    label_result = ''
    for a in params:
        label_result = label_result + a + '\n'

    S0 = 0.96
    k_new = np.sqrt(k_fit**2-0.2625*e)
    exp_new = interpolation(exp,k_new)/S0


    plt.figure()
    plt.plot(k_fit,exp_new*k_fit**2,'k',linewidth=0.5)
    plt.plot(k_fit,fit_LLE,linewidth=0.5,label = 'LL = {:.2f} \n'.format(LL) +'parameter = {:s}'.format(label_result))
    plt.legend()
    plt.xlabel('K $(\AA^{-1})$')
    plt.ylabel('$\chi(k)k^2$')
    plt.show()
# calc_Feff_phase()


###############################################################
##
#
#  used calculated phase from calc_phase() to plot the curve
#
#
###################################################################
def use_calc_4_fit(shell):
    k_fine,chi,phase = calc_phase(shell)
    if shell == 1:
        R = 2.45
        N = 4
    elif shell == 2:
        R = 4.0
        N = 12
    elif shell == 3:
        R = 4.69
        N = 12

    fit = N*np.abs(eval('Feff'+str(shell))) / (k_fit * R ** 2) \
          * np.sin(2 * k_fit * R + phase) \
          * np.exp(-2 * R /  eval('lambda'+str(shell))) \
          * np.exp(-2 * 0.0054 * k_fit ** 2)
    plt.plot(k_fit,fit*k_fit**weight)
    plt.plot(k_fine,chi*k_fine**weight,label='experiment',linewidth='0.5')
    plt.show()


use_calc_4_fit(3)




