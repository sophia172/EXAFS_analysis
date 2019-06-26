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

k_min = (np.abs(k - 2.5)).argmin()
k_max = (np.abs(k - 18)).argmin()
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
phase1 =  interpolation(phase1,  k_fit)
phase2 =  interpolation(phase2,  k_fit)
phase3 =  interpolation(phase3,  k_fit)
e = 13.8
ss = 0.0054
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
    data = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Experiment_Analysis/Ge/c_ge12_chi.xmu')
    data1 = data[:,1]
    data2 = data[:,2]
    data3 = data[:,3]
    k_fine = np.linspace(3,18,25601)
    # interpolate each file with k_fine
    k = np.sqrt((data[:,0]+e)/3.81)
    chi1= data1
    exp = np.vstack((k,chi1)).transpose()
    chi1 = interpolation(exp,k_fine)


    chi2= data2
    exp = np.vstack((k,chi2)).transpose()
    chi2 = interpolation(exp,k_fine)

    chi3= data3
    exp = np.vstack((k,chi3)).transpose()
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
    correction = 2* pi *5



    peak_index = find_peaks(chi,prominence=0.0001,distance=10)[0]
    valley_index = find_peaks(-1*chi,prominence=0.0001,distance=10)[0]
    zero_index = np.where(np.diff(np.sign(chi)))[0]


    # Figure out where the curve start
    index = np.hstack((peak_index,valley_index,zero_index))
    index = np.sort(index)
    # print(index)
    phase_knot = np.array([])
    phase_peak = np.array([])
    if index[0] in peak_index:
        # print('curve start with peak',index)
        for n,i in enumerate(index):
            phase_knot = np.append(phase_knot, (pi/2 *(n+1)+ correction - 2 * k_fine[i] *R))
            # print('2kR: ',2 * k_fine[i] *R)
    elif index[0] in valley_index:
        # print('curve start with valley')
        for n,i in enumerate(index):
            phase_knot = np.append(phase_knot, pi/2 *(n+3) - 2 * k_fine[i] *R + correction)
    elif index[1] in peak_index:
        # print('curve start with 0 then peak')
        for n,i in enumerate(index):
            phase_knot = np.append(phase_knot, pi/2 *n - 2 * k_fine[i] *R + correction)
    elif index[1] in valley_index:
        # print('curve start with 0 then valley')
        for n,i in enumerate(index):
            phase_knot = np.append(phase_knot, pi/2 *(n+2) - 2 * k_fine[i] *R + correction)
    k_knot = k_fine[index]
    for i in peak_index:
        # print(np.where(index == i))

        phase_peak = np.append(phase_peak,phase_knot[np.where(index == i)])
    # print(index)
    k_peak = k_fine[peak_index]
    chi_peak = chi[peak_index]
    k_valley = k_fine[valley_index]
    chi_valley = chi[valley_index]
    k_zero = k_fine[zero_index]
    chi_zero = chi[zero_index]



    phase_FEFF = np.vstack((k_fit,eval('phase'+ str(shell)))).transpose()
    phase_FEFF = interpolation(phase_FEFF,k_knot)

    # def model_phaseWOFEFF(x,a,b,c,d):
    #     return a*x**3 + b*x**2 + c*x + d
    # def model_phase(x,a,b,c,d):
    #     return model_phaseWOFEFF(x,a,b,c,d) + phase_FEFF
    #
    # popt, pcov = curve_fit(model_phase,k_knot,phase_knot)
    # global phase
    # phase = model_phaseWOFEFF(k_fit,*popt[:4]) + eval('phase'+ str(shell))
    # print('phase: ',*popt)

    spline_knot = int(len(k_knot)/6)
    t= [k_knot[spline_knot],k_knot[spline_knot*2],k_knot[spline_knot*3],k_knot[spline_knot*4],k_knot[spline_knot*5]]
    spl = LSQUnivariateSpline(k_knot,phase_knot/phase_FEFF,t)
    phase = spl(k_fit)*eval('phase'+str(shell))
    ############################################# show the experiemnt point we used to fit
    plt.figure()
    plt.plot(k_fine,chi,label='exp')
    plt.plot(k_peak,chi_peak,'.',label='peak')
    plt.plot(k_zero,chi_zero,'.',label='zero')
    plt.plot(k_valley,chi_valley,'.',label='valley')
    plt.legend()
    plt.xlabel('k ($\AA^{-1}$)')
    plt.title('find knot point in shell %d'%shell)
    plt.show()
    ###
    #       Plot calculated phase compared with FEFF calculated one
    ###
    ############################################# show the difference between calcuated phase and FEFF phase
    plt.figure()
    plt.plot(k_knot,phase_knot,'.')
    # plt.plot(k_fine_peak,phase_peak,'.',label='peak')
    # plt.plot(k_fine_valley,phase_valley,'.',label='valley')
    # plt.plot(k_fine_zero,phase_zero,'.',label='zero')
    plt.plot(k_fit,phase,label='phase')
    plt.plot(k_fit,eval('phase'+str(shell)),label='previous calculated phase')
    # plt.xlim([3,6])
    # plt.ylim([30,40])
    plt.xlabel('k ($\AA^{-1}$)')
    plt.title('phase fit in shell %d'%shell)
    plt.legend()
    plt.show()

    ###
    #       Plot chi with calculated phase and FEFF result
    ###
    return phase

###############################################################
##
#
#  Calculate Feff from experimental data
#
#  1. Use e0 alligned Chi result
#  2. Withdraw knot position and calculate phase
#
###################################################################
def calc_Feff(shell):
    data = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Experiment_Analysis/Ge/c_ge12_chi.xmu')
    data1 = data[:,1]
    data2 = data[:,2]
    data3 = data[:,3]
    k_fine = np.linspace(3,18,25601)
    # interpolate each file with k_fine
    k = np.sqrt((data[:,0]+e)/3.81)
    chi1= data1
    exp = np.vstack((k,chi1)).transpose()
    chi1_fine = interpolation(exp,k_fine)


    chi2= data2
    exp = np.vstack((k,chi2)).transpose()
    chi2_fine = interpolation(exp,k_fine)

    chi3= data3
    exp = np.vstack((k,chi3)).transpose()
    chi3_fine = interpolation(exp,k_fine)

    if shell == 1:
        chi_fine = chi1_fine
        chi = chi1
        R = 2.45
        N = 4

    elif shell == 2:
        chi_fine = chi2_fine-chi1_fine
        chi = chi2-chi1
        R = 4.0
        N = 12

    elif shell == 3:
        chi_fine = chi3_fine-chi2_fine
        chi = chi3-chi2
        R = 4.69
        N = 12

    correction = 2* pi *5
    # print(k1,k_fit)
    chi_fit = interpolation(np.vstack((k,chi)).transpose(),k_fit)
    peak_index = find_peaks(abs(chi_fine),prominence=0.0001,distance=10)[0]
    k_peak = k_fine[peak_index]
    chi_peak = chi_fine[peak_index]

    phase_FEFF = np.vstack((k_fit,eval('phase'+ str(shell)))).transpose()
    phase_FEFF = interpolation(phase_FEFF,k_peak)
    lambda_FEFF = np.vstack((k_fit,eval('lambda'+ str(shell)))).transpose()
    lambda_FEFF = interpolation(lambda_FEFF,k_peak)
    Feff_FEFF = np.vstack((k_fit,eval('Feff'+ str(shell)))).transpose()
    Feff_FEFF = interpolation(Feff_FEFF,k_peak)

    def EXAFSmodelWOfeff(k_fit,phase,lambdaa):
        return N * 1 / (k_fit * R ** 2) \
               * np.sin(2 * k_fit * R +  phase) \
               * np.exp(-2 * R /  lambdaa) \
               * np.exp(-2 * ss * k_fit ** 2)
    # print(shape(k_fit),shape(phase1),shape(lambda1))
    # print(chi_fit)
    spline_knot = int(len(k_peak)/6)
    t= [k_peak[spline_knot],k_peak[spline_knot*2],k_peak[spline_knot*3],k_peak[spline_knot*4],k_peak[spline_knot*5]]
    spl = LSQUnivariateSpline(k_peak,chi_peak/EXAFSmodelWOfeff(k_peak,phase_FEFF,lambda_FEFF),t)
    Feff = spl(k_fit)

    print('Feff:  ',Feff)


    ############################################# show the experiemnt point we used to fit
    plt.figure()
    plt.plot(k_fine,chi_fine,label='exp')
    plt.plot(k_peak,chi_peak,'.',label='peak')

    plt.legend()
    plt.xlabel('k ($\AA^{-1}$)')
    plt.title('find Feff point in shell %d'%shell)
    plt.show()

    ############################################# show the difference between calcuated Feff and FEFF feff
    plt.figure()
    plt.plot(k_peak,chi_peak/EXAFSmodelWOfeff(k_peak,phase_FEFF,lambda_FEFF),'.')
    plt.plot(k_fit,Feff,label='calculated Feff')
    plt.plot(k_fit,eval('Feff'+str(shell)),label='previous calcuated Feff')
    # plt.xlim([3,6])
    # plt.ylim([30,40])
    plt.xlabel('k ($\AA^{-1}$)')
    plt.title('Feff fit in shell %d'%shell)
    plt.legend()
    plt.show()

    ############################################# show the result from calculated Feff
    plt.figure()
    plt.plot(k,chi,label='experiment')
    plt.plot(k_fit,EXAFSmodelWOfeff(k_fit,eval('phase'+str(shell)),eval('lambda'+str(shell)))*Feff)
    plt.xlabel('k ($\AA^{-1}$)')
    plt.title('Feff fit in shell %d'%shell)
    plt.legend()
    plt.show()

    return Feff

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
bounds_dic = {'Ge1': np.array([(3.999, 4), (2.44, 2.46), (0.0001, 0.05)]),
              'Ge2': np.array([(11.999,12), (3.99, 4.01), (0.002, 0.05)]),
              'Ge3': np.array([(11.999,12), (4.68, 4.70), (0.002, 0.05)]),
              'e': np.array([(0,20)]),
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



def Ge_shell(shell):
    data = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Experiment_Analysis/Ge/c_ge12_chi.xmu')
    data1 = data[:,1]
    data2 = data[:,2]
    data3 = data[:,3]

    k_fine = np.linspace(3,18,12801)
    # interpolate each file with k_fine

    k = np.sqrt(data[:,0]/3.81)
    if shell == 1:
        chi = data1
        exp = np.vstack((k,chi)).transpose()
    elif shell == 2:
        chi = data2-data1
        exp = np.vstack((k,chi)).transpose()
    elif shell == 3:
        chi = data3-data2
        exp = np.vstack((k,chi)).transpose()


    def model_shell(x,N,R,ss):
        return N * x ** weight \
               * np.abs( eval('Feff'+str(shell))) / (x * R ** 2) \
               * np.sin(2 * x * R +  eval('phase'+str(shell))) \
               * np.exp(-2 * R /  eval('lambda'+str(shell))) \
               * np.exp(-2 * ss * x ** 2)



    bounds = np.vstack((bounds_dic['Ge'+str(shell)], bounds_dic['S0'], bounds_dic['e']))


    def LLE_model(params):
        global exp,model
        e = params[-1]
        S0 = params[-2]
        # print('ss',params[-3])
        # print('k',k)
        k_new = np.sqrt(k**2+0.2625*e)
        # print('k_new',k_new)
        g = interp1d(k_new,chi)
        exp_fit = g(k_fit.tolist())/S0
        model_LLE = model_shell(k_fit, *params[:-2])
        LL = -np.sum(stats.norm.logpdf(exp_fit * k_fit**weight, loc=model_LLE))
        return LL


    LLE_fitresult = differential_evolution(LLE_model,bounds=bounds,strategy='best1exp', popsize=200,maxiter=10000)
    fit_LLE = model_shell(k_fit, *LLE_fitresult.x[:-2])
    LL = LLE_model(LLE_fitresult.x)
    e=LLE_fitresult.x[-1]
    S0=LLE_fitresult.x[-2]
    ss=LLE_fitresult.x[-3]
    print(e,S0,ss)
    k_new=np.sqrt(k**2+0.2625*e)
    g = interp1d(k_new,chi)
    exp_fit = g(k_fit.tolist())


    plt.plot(k_fit, exp_fit * k_fit **weight,'k',linewidth=0.5)
    plt.plot(k_fit,fit_LLE*S0,linewidth=0.5,label='shell' + str(shell) + '\n' + 'e= {:.2f}   S0 = {:.2f}   ss = {:.5f}   LL = {:.2F}\n\n'.format(e,S0,ss,LL) )
    plt.plot(k_fit,fit_LLE*S0-exp_fit*k_fit**weight -2,linewidth=0.5,label='difference')
    plt.xlim([2,26])
    plt.xlabel('K $(\AA^{-1})$')
    plt.ylabel('$\chi(k)k^2$')
    plt.legend()
    plt.tight_layout()
    data = np.vstack((k_fit,exp_fit*k_fit**weight,fit_LLE*S0,exp_fit*k_fit**weight-fit_LLE*S0)).transpose()
    np.savetxt('/Users/Sophia/ownCloud/PhD/Simulation/Ge/Ge_fit_shell'+str(shell)+'.dat',data,header='k  experiment_k^2 fit_k^2 difference')

    plt.savefig('/Users/Sophia/ownCloud/PhD/Statistic Analysis/figure/Ge_shell'+str(shell)+'fit.pdf',format='pdf')
    plt.show()
    return e,ss


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
    # LLE_fitresult = differential_evolution(LLE_model_FEFF,bounds=bounds_FEFF,strategy='randtobest1exp', popsize=10,maxiter=10000)
    LLE_fitresult = differential_evolution(LLE_model_FEFF,bounds=bounds_FEFF)
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
    plt.figure()
    plt.plot(k_fit,fit*k_fit**weight)
    plt.plot(k_fine,chi*k_fine**weight,label='experiment',linewidth='0.5')
    plt.xlabel('k ($\AA^{-1}$)')
    plt.title('$\chi$ result with calculated phase in shell %d'%shell)
    plt.show()


# use_calc_4_fit(3)

#

for i in range(2):
    if i ==0:
        phase1 = phase1  # To start the convergence, I need to set FEFF output to start the loop
        Feff1 = calc_Feff(1)
        e,ss = Ge_shell(1)
    else:
        phase1 = calc_phase(1)
        Feff1 = calc_Feff(1)
        e,ss = Ge_shell(1)
data = np.vstack((k_fit,Feff1,phase1,lambda1)).transpose()
np.savetxt('/Users/Sophia/ownCloud/PhD/Simulation/Ge/FEFF1fromEXP.dat',data,header='k,Feff,phase,lambda')



