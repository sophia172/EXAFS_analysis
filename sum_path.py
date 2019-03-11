import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
# define model function

def model(k,S0,feff, lambdaN, phase,N,R,ss):
    return S0**2 * N * abs(feff.astype(float))/(k.astype(float)*R**2)  *  np.sin(2*k.astype(float)*R + phase.astype(float))  *  np.exp(-2*R/lambdaN.astype(float))  *  np.exp(-2* ss * k.astype(float)**2)

def FEFF_model(k,S0,feff, lambdaN, phase,N,R,ss):
    return S0**2 * N * abs(feff.astype(float))/(k.astype(float)*R**2)  *  np.sin(2*k.astype(float)*R + phase.astype(float))  *  np.exp(-2*R/lambdaN.astype(float))  *  np.exp(-2* ss * k.astype(float)**2)

def interpolation(f,x):
    g = interp1d(f[:,0],f[:,1])
    return g(x.tolist())


##################################
# Import data
######################################
# import  path data
data = np.array(pd.read_excel('/Users/Sophia/ownCloud/PhD/Statistic Analysis/Parameter fro EXAFS anaylysis.xlsx'))
data1 = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Simulation/CdS/FEFF/EXAFS/Cd/F-43m/feff0001.dat',skip_header=15)
data2 = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Simulation/CdS/FEFF/EXAFS/Cd/F-43m/feff0002.dat',skip_header=15)
# CdO path data
data3 = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Simulation/CdS/FEFF/XANES/CdO/feff0001.dat',skip_header=15)
#  FEFF calculation data using bulk crystal
FEFF_xmu = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Simulation/CdS/FEFF/EXAFS/Cd/F-43m/xmu.dat')
#  experiment data
experiment = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Experiment_Analysis/CdS/CdS_XAFS_CdK_chi_2018.txt')
# artemis data fiT
artemis_bulk = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Statistic Analysis/data/CdS_bulk_FEFF.k')
artemis_M311 = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Statistic Analysis/data/CdS_M311_FEFF13.k')





# feed data files to variable
Feff1 = data1[:,(0,2)]
Feff2 = data2[:,(0,2)]
Feff3 = data3[:,(0,2)]

phase1 = data1[:,(0,1)]
phase1[:,1] += data1[:,3]

phase2 = data2[:,(0,1)]
phase2[:,1] += data2[:,3]

phase3 = data3[:,(0,1)]
phase3[:,1] += data3[:,3]

lambda1 = data1[:,(0,5)]
lambda2 = data2[:,(0,5)]
lambda3 = data3[:,(0,5)]

S01 = data1[:,(0,4)]
S02 = data2[:,(0,4)]
S03 = data3[:,(0,4)]

k = data[1:361,11]

# interpolate given parameter files
Feff1 = interpolation(Feff1,k)
Feff2 = interpolation(Feff2,k)
Feff3 = interpolation(Feff3,k)
lambda1 = interpolation(lambda1,k)
lambda2 = interpolation(lambda2,k)
lambda3 = interpolation(lambda3,k)
phase1 = interpolation(phase1,k)
phase2 = interpolation(phase2,k)
phase3 = interpolation(phase3,k)
S01 = interpolation(S01,k)
S02 = interpolation(S02,k)
S03 = interpolation(S03,k)



# import experiment data and FEFF result

M311 = data[1:361,12]
M322= data[1:361,13]
bulk = experiment[:,(0,2)]
FEFF = FEFF_xmu[:,(2,5)]  # include both x and y axis
ar_bulk = artemis_bulk[:,(0,2)]
ar_M311 = artemis_M311[:,(0,2)]



                         # bulk
N1 = 4                 #   4
N2 = 12                  #   12
N3 = 3.5
R1 = 2.6164             #    2.54      Cd-S bond
R2 = 4.2726               #   4.18     Cd-Cd bond
R3 = 2.45           # This is Cd-O bond
ss1 = 0.00              #    0.005
ss2 = 0.0              #     0.014
ss3 = 0.015
S01 = S02 = S03 = 0.864                 #     1



# Non-linear curve fitting
# def model_fit1(k,N,R,ss):
#     return N * abs(Feff1.astype(float))/(k.astype(float)*R**2)  *  np.sin(2*k.astype(float)*R + phase1.astype(float))  *  np.exp(-2*R/lambda1.astype(float))  *  np.exp(-2* ss * k.astype(float)**2)
#
# def model_fit3(k,N,R,ss):
#     return N * abs(Feff3.astype(float))/(k.astype(float)*R**2)  *  np.sin(2*k.astype(float)*R + phase3.astype(float))  *  np.exp(-2*R/lambda3.astype(float))  *  np.exp(-2* ss * k.astype(float)**2)
#
# popt, pcov = curve_fit(model_fit1(k,3.2, 2.51, 0.005)+model_fit3(k,3.5, 2.45, 0.015),k.astype(float),M311.astype(float))


space = 0
plt.figure()
plt.plot(k,M311,'k.',markersize=2,label = 'experiment M311')
# plt.plot(bulk[:,0],bulk[:,1],'k.',markersize = 2,label = 'experiment bulk')
# plt.plot(FEFF[:,0],FEFF[:,1]-space,'r--',label= 'FEFF')
# plt.plot(ar_bulk[:,0],ar_bulk[:,1],'g-',label = 'artemis_fit')
plt.plot(ar_M311[:,0],ar_M311[:,1],'g-',label = 'artemis_fit')
# plt.plot(k,model(k,S01,Feff1,lambda1,phase1,N1,R1,ss1)-2*space,'r-',label='model_1')
# plt.plot(k,model(k,S02,Feff2,lambda2,phase2,N2,R2,ss2)-2*space,'g-',label='model_2')
# plt.plot(k,model(k,S01,Feff1,lambda1,phase1,N1,R1,ss1) + model(k,S02,Feff2,lambda2,phase2,N2,R2,ss2)-3*space,'g-',label='model_1 + model_2')
# plt.plot(k,model(k,S01,Feff1,lambda1,phase1,N1,R1,ss1) + model(k,S03,Feff3,lambda3,phase3,N3,R3,ss3)-3*space,'r-',label='model_1 + model_3')
# plt.plot(k,model(k,S01,Feff1,lambda1,phase1,N1,R1,ss1) + model(k,S02,Feff2,lambda2,phase2,N2,R2,ss2) + model(k,S03,Feff3,lambda3,phase3,N3,R3,ss3)-3*space,'g-',label='model_1 + model_2 + model_3')

plt.xlabel('k')
plt.ylabel('chi')
plt.title('model not includ S0(k) & real[p]')
plt.ylim([-0.1,0.1])
plt.legend()
plt.show()