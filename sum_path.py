import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
# define model function

def model(k,S0,feff, lambdaN, phase,N,R,ss):
    return k**2 * S0**2 * N \
           * abs(feff.astype(float))/(k.astype(float)*R**2)  \
           * np.sin(2*k.astype(float)*R + phase.astype(float))  \
           * np.exp(-2*R/lambdaN.astype(float))  \
           * np.exp(-2* ss * k.astype(float)**2)

def FEFF_model(k,S0,feff, lambdaN, phase,N,R,ss):
    return k**2 * S0**2 * N \
           * abs(feff.astype(float))/(k.astype(float)*R**2)  \
           *  np.sin(2*k.astype(float)*R + phase.astype(float))  \
           *  np.exp(-2*R/lambdaN.astype(float))  \
           *  np.exp(-2* ss * k.astype(float)**2)

def interpolation(f,x):
    g = interp1d(f[:,0],f[:,1])
    return g(x.tolist())

def k_range(a,b):
    idx_a = (np.abs(experiment[:,0] - a)).argmin()
    idx_b = (np.abs(experiment[:,0] - b)).argmin()
    return idx_a,idx_b

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


k_min,k_max = k_range(2,4)
k = experiment[k_min:k_max,0]
print('k = ',k)
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

M311 = k**2 * experiment[k_min:k_max,3]
M322= k**2 * experiment[k_min:k_max,4]
bulk = k**2 * experiment[k_min:k_max,2]
FEFF = FEFF_xmu[:,(2,5)]  # include both x and y axis
ar_bulk = artemis_bulk[:,(0,2)]
ar_M311 = artemis_M311[:,(0,2)]

# print(k,M311)

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
def model_fit123(k,N1,R1,ss1,N2,R2,ss2,N3,R3,ss3):
     return N1 * k**2 \
            * abs(Feff1.astype(float))/(k.astype(float)*R1**2)  \
            * np.sin(2 * k.astype(float) * R1 + phase1.astype(float))  \
            * np.exp(-2 * R1 / lambda1.astype(float))  \
            * np.exp(-2 * ss1 * k.astype(float)**2) \
            + N2 *  k**2 \
            * abs(Feff2.astype(float)) / (k.astype(float)*R2**2) \
            * np.sin(2 * k.astype(float)*R2 + phase2.astype(float)) \
            * np.exp(-2 * R2 / lambda2.astype(float)) \
            * np.exp(-2 * ss2 * k.astype(float)**2) \
            + N3 *  k**2 \
            * abs(Feff3.astype(float)) / (k.astype(float)*R3**2) \
            * np.sin(2 * k.astype(float)*R3 + phase3.astype(float)) \
            * np.exp(-2 * R3 / lambda3.astype(float)) \
            * np.exp(-2 * ss3 * k.astype(float)**2)

def model_fit12(k,N1,R1,ss1,N2,R2,ss2):
    return N1 * k**2 \
           * abs(Feff1.astype(float))/(k.astype(float)*R1**2) \
           * np.sin(2 * k.astype(float) * R1 + phase1.astype(float)) \
           * np.exp(-2 * R1 / lambda1.astype(float)) \
           * np.exp(-2 * ss1 * k.astype(float)**2) \
           + N2 *  k**2 \
           * abs(Feff2.astype(float)) / (k.astype(float)*R2**2) \
           * np.sin(2 * k.astype(float)*R2 + phase2.astype(float)) \
           * np.exp(-2 * R2 / lambda2.astype(float)) \
           * np.exp(-2 * ss2 * k.astype(float)**2)

def model_fit13(k,N1,R1,ss1,N3,R3,ss3):
    return N1 * k**2 \
           * abs(Feff1.astype(float))/(k.astype(float)*R1**2) \
           * np.sin(2 * k.astype(float) * R1 + phase1.astype(float)) \
           * np.exp(-2 * R1 / lambda1.astype(float)) \
           * np.exp(-2 * ss1 * k.astype(float)**2) \
           + N3 *  k**2 \
           * abs(Feff3.astype(float)) / (k.astype(float)*R3**2) \
           * np.sin(2 * k.astype(float)*R3 + phase3.astype(float)) \
           * np.exp(-2 * R3 / lambda3.astype(float)) \
           * np.exp(-2 * ss3 * k.astype(float)**2)


#
# def model_fit3(k,N3,R3,ss3):
#     return N3 * abs(Feff3.astype(float))/(k.astype(float)*R3**2)  *  np.sin(2*k.astype(float)*R3 + phase3.astype(float))  *  np.exp(-2*R3/lambda3.astype(float))  *  np.exp(-2* ss3 * k.astype(float)**2)
# #
# popt, pcov = curve_fit(model_fit12,k.astype(float),
#                        M311.astype(float),bounds=([1,2.3,0.002,0,3.8,0.01], [4,2.7,0.02,12,4.3,0.2]))
#
# popt, pcov = curve_fit(model_fit13,k.astype(float),
#                        M311.astype(float),bounds=([1,2.3,0.002,0,2.3,0.002], [4,2.7,0.02,6,2.6,0.02]))
#
popt, pcov = curve_fit(model_fit123,k.astype(float),
                       M311.astype(float),bounds=([1,2.3,0.002,0,3.8,0.01,0,2.3,0.002], [4,2.7,0.02,12,4.3,0.2,6,2.6,0.02]))
# #
#



perr = np.sqrt(np.diag(pcov))
print('parameters are', popt)
print('standard deviation error ',perr)
space = 0
plt.figure()
plt.plot(k,M311,'k.',markersize=2,label = 'experiment M311')
# plt.plot(bulk[:,0],bulk[:,1],'k.',markersize = 2,label = 'experiment bulk')
# plt.plot(FEFF[:,0],FEFF[:,1]-space,'r--',label= 'FEFF')
# plt.plot(ar_bulk[:,0],ar_bulk[:,1],'g-',label = 'artemis_fit')
# plt.plot(ar_M311[:,0],ar_M311[:,0]**2*ar_M311[:,1],'r-',label = 'artemis_fit')
plt.plot(k,model_fit123(k,*popt),'g-',label = 'curve_fit(model 1+2+3)')
# plt.plot(k,model(k,S01,Feff1,lambda1,phase1,N1,R1,ss1)-2*space,'r-',label='model_1')
# plt.plot(k,model(k,S02,Feff2,lambda2,phase2,N2,R2,ss2)-2*space,'g-',label='model_2')
# plt.plot(k,model(k,S01,Feff1,lambda1,phase1,N1,R1,ss1) + model(k,S02,Feff2,lambda2,phase2,N2,R2,ss2)-3*space,'g-',label='model_1 + model_2')
# plt.plot(k,model(k,S01,Feff1,lambda1,phase1,N1,R1,ss1) + model(k,S03,Feff3,lambda3,phase3,N3,R3,ss3)-3*space,'r-',label='model_1 + model_3')
# plt.plot(k,model(k,S01,Feff1,lambda1,phase1,N1,R1,ss1) + model(k,S02,Feff2,lambda2,phase2,N2,R2,ss2) + model(k,S03,Feff3,lambda3,phase3,N3,R3,ss3)-3*space,'g-',label='model_1 + model_2 + model_3')

plt.xlabel('k')
plt.ylabel('chi')
plt.title('model not includ S0(k) & real[p]')
# plt.ylim([-0.1,0.1])
plt.legend()
plt.show()

# 1+3  parameters are [4.         2.50032439 0.00596682 4.91423657 2.35566526 0.02      ]
# standard deviation error  [0.30151865 0.00147638 0.00038155 0.31684124 0.01283126 0.00241257]

# 1+2+3   parameters are [4.         2.50039311 0.00596238 2.64859949 3.96858409 0.01859781
# 4.88750536 2.35582555 0.02      ]
# standard deviation error  [2.93560624e-01 1.43814854e-03 3.71448389e-04 1.02322068e+00
# 1.91052239e-02 4.67219681e-03 3.09343969e-01 1.25697120e-02