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
exp = experiments[:,1]
weight=2

k_min = (np.abs(k - 3)).argmin()
k_max = (np.abs(k - 18)).argmin()
k = k[k_min:k_max]
exp = exp[k_min:k_max]

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

Feff1 =  interpolation(Feff1,  k)
Feff2 =  interpolation(Feff2,  k)
Feff3 =  interpolation(Feff3,  k)
lambda1 =  interpolation(lambda1,  k)
lambda2 =  interpolation(lambda2,  k)
lambda3 =  interpolation(lambda3,  k)
phase1 =  interpolation(phase1,  k)
phase2 =  interpolation(phase2,  k)
phase3 =  interpolation(phase3,  k)
#
# plt.plot(k,phase3)
# plt.show()
print('N1              R1                  ss1                    N2                    R2             ss2                  N3                       R3                  ss3')
def EXAFS_model(x,N1,R1,ss1,N2,R2,ss2,N3,R3,ss3,e):
    print(N1,R1,ss1,N2,R2,ss2,N3,R3,ss3)
    return e + N1 * x **weight\
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

def EXAFS_shell(x,N,R,ss):
    return e + N * x **weight \
           * abs( Feff1.astype(float)) / (x.astype(float) * R ** 2) \
           * np.sin(2 * x.astype(float) * R +  phase1.astype(float)) \
           * np.exp(-2 * R /  lambda1.astype(float)) \
           * np.exp(-2 * ss * x.astype(float) ** 2)


params_dic = {'Ge1': [6, 2.35, 0.00234], 'Ge2': [12, 4, 0.005], 'Ge3': [12, 4.8, 0.006], 'e':[3], 'sd': [1]}
bounds_dic = {'Ge1': np.array([(3.3, 3.4), (2.1, 3), (0.0023, 0.0024)]),
              'Ge2': np.array([(9.9,10), (3., 4.3), (0.0031, 0.0032)]),
              'Ge3': np.array([(9.9,10), (4.3, 5.2), (0.0031, 0.0032)]),
              'e': np.array([(-6,6)]),
              'sd': np.array([(0,1)])}
bounds_shell= np.vstack((bounds_dic['Ge1'], bounds_dic['e']))

params = params_dic['Ge1'] + params_dic['Ge2'] + params_dic['Ge3'] + params_dic['e']
bounds = np.vstack((bounds_dic['Ge1'], bounds_dic['Ge2'], bounds_dic['Ge3'], bounds_dic['e']))
# bounds_NL = np.vstack((bounds_dic['Ge1'], bounds_dic['Ge2'], bounds_dic['Ge3'], bounds_dic['e']))
# popt, pcov = curve_fit(EXAFS_model, k.astype(float),(exp * k **2).astype(float), p0=params,bounds=bounds_NL.transpose())
# fit = EXAFS_model(k, *popt)
# print(*popt)

def LLE_model(params):
    model = EXAFS_model(k, *params)
    LL = -np.sum(stats.norm.logpdf(exp * k**weight, loc=model))
    return LL

def LLE_shell(params):
    model = EXAFS_shell(k, *params)
    LL = -np.sum(stats.norm.logpdf(exp * k**weight, loc=model))
    return LL

class MyTakeStep(object):
    def __init__(self, stepsize=0.001):
        self.stepsize = stepsize
    def __call__(self, x):
        s = self.stepsize
        x[0] += np.random.uniform(-2.*s, 2.*s)
        x[1:] += np.random.uniform(-s, s, x[1:].shape)
        return x






mytakestep = MyTakeStep()
minimizer_kwargs = {"method":"L-BFGS-B", "jac":True}
# Local minimum
# LLE_fitresult = minimize(LLE_model, params, method='TNC', bounds=bounds,options={'maxiter': 12000})
# Global minimum
# LLE_fitresult = basinhopping(LLE_model, params,niter=200, take_step=mytakestep)
# Global minimu rate 4/5
LLE_fitresult = differential_evolution(LLE_model,bounds=bounds,strategy='randtobest1exp', popsize=200,maxiter=10000)

# LLE_fitresult = shgo(LLE_model,bounds=bounds,iters=5)

# LLE_fitresult = dual_annealing(LLE_model,bounds=bounds)

fit_LLE = EXAFS_model(k, *LLE_fitresult.x)
print(*LLE_fitresult.x)
####################################
#
#       Import fit result from Artemis
#           and plot
#
######################################################

chi_data = np.genfromtxt("/Users/Sophia/ownCloud/PhD/Simulation/Ge/c_ge12_chi.k2")


plt.figure()
plt.subplot(2,1,1)
plt.plot(k,exp *k**weight,'k', linewidth = 1, label = 'experiment')
# plt.plot(k, fit, linewidth=1, label='nonlinear curve fit (scipy) of Germanium')
plt.plot(k,fit_LLE,linewidth=1, label='Python fit of Ge')
# plt.plot(k, fit, linewidth=1, label='nonlinear curve fit (scipy) of Germanium')
plt.plot(chi_data[:,0],chi_data[:,2],linewidth=1, label='Artemis fit of Ge')
plt.xlim([2,18])
plt.ylim([-3,3])
plt.xlabel('k ($\AA^{-1}$)')
plt.ylabel('$\chi$')
plt.legend()

plt.subplot(2,1,2)
r,amp,real,imag = calcfft(k,exp*k**3,kmin=3,kmax=18)
plt.plot(r,amp,'k', linewidth = 1, label = 'experiment')
# r,amp,real,imag = calcfft(k,fit*k,kmin=3,kmax=18)
# plt.plot(r,amp, linewidth=1, label='nonlinear curve fit (scipy) of Germanium')
r,amp,real,imag = calcfft(k,fit_LLE*k,kmin=3,kmax=18)
plt.plot(r,amp,linewidth=1, label='Python fit of Ge')

r,amp,real,imag = calcfft(chi_data[:,0],chi_data[:,2]*chi_data[:,0],kmin=3,kmax=18)
plt.plot(r,amp,linewidth=1, label='Artemis fit of Ge')
plt.xlabel('R ($\AA$)')
plt.ylabel('Intensity')
plt.legend()
plt.xlim([0,5])
plt.savefig('/Users/Sophia/ownCloud/PhD/Statistic Analysis/figure/Ge_analysis.pdf',format='pdf')
plt.show()