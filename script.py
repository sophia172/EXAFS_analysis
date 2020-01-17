import numpy as np
from matplotlib import pyplot as plt
from xafs import *

# mu_data = np.genfromtxt('/Users/Sophia/Downloads/drive-download-20190501T161928Z-001/mu.nor')
# mu_energy = mu_data[:,0]
# mu_pt_centre = mu_data[:,1]
# mu_pt_edge = mu_data[:,2]
# mu_pto2 = mu_data[:,3]
# mu_pt = mu_data[:,4]
# mu_ptcl2 = mu_data[:,5]
# mu_pth2_centre = mu_data[:,6]
# mu_pth2_edge = mu_data[:,7]
# plt.figure()
# plt.plot(mu_energy[340:600],mu_data[340:600,1],linewidth=0.5,label='Pt centre')
# plt.plot(mu_energy[340:600],mu_data[340:600,2],linewidth=0.5,label='Pt edge')
# plt.plot(mu_energy[340:600],mu_data[340:600,3],linewidth=0.5,label='PtO_2')
# plt.plot(mu_energy[340:600],mu_data[340:600,4],linewidth=0.5,label='Pt foil')
# plt.plot(mu_energy[340:600],mu_data[340:600,5],linewidth=0.5,label='PtCl_2')
# plt.plot(mu_energy[340:600],mu_data[340:600,6],linewidth=0.5,label='Pt H_2 centre')
# plt.plot(mu_energy[340:600],mu_data[340:600,7],linewidth=0.5,label='Pt H_2 edge')
# plt.xlabel('Energy (eV)')
# plt.ylabel('mu')
# plt.legend()
# plt.title('normalized mu data')
# plt.savefig('/Users/Sophia/Downloads/drive-download-20190501T161928Z-001/Pt_xanes.pdf',format='pdf')
# plt.show()
#
#
#
#
#
# chi_data = np.genfromtxt('/Users/Sophia/Downloads/drive-download-20190501T161928Z-001/chi.chik2')
# chi_k = chi_data[:,0]
# chi_pt_centre = chi_data[:,2]
# chi_pt_edge = chi_data[:,3]
# chi_pto2 = chi_data[:,4]
# chi_pt = chi_data[:,5]
# chi_ptcl2 = chi_data[:,6]
# chi_pth2_centre = chi_data[:,7]
# chi_pth2_edge = chi_data[:,8]
# plt.figure()
# # plt.plot(chi_k,chi_data[:,2],linewidth=0.5,label='Pt centre')
# # plt.plot(chi_k,chi_data[:,3],linewidth=0.5,label='Pt edge')
# plt.plot(chi_k,chi_data[:,4],linewidth=0.5,label='PtO_2')
# plt.plot(chi_k,chi_data[:,5],linewidth=0.5,label='Pt foil')
# plt.plot(chi_k,chi_data[:,6],linewidth=0.5,label='PtCl_2')
# # plt.plot(chi_k,chi_data[:,7],linewidth=0.5,label='Pt H_2 centre')
# # plt.plot(chi_k,chi_data[:,8],linewidth=0.5,label='Pt H_2 edge')
# plt.legend()
# plt.xlabel('K')
# plt.ylabel('chi')
# plt.savefig('/Users/Sophia/Downloads/drive-download-20190501T161928Z-001/Pt_chi_reference',format='pdf')
# plt.show()
#
#
#
# FT_data = np.genfromtxt('/Users/Sophia/Downloads/drive-download-20190501T161928Z-001/FT.chir_mag')
# FT_r = FT_data[:,0]
# FT_pt_centre = FT_data[:,1]
# FT_pt_edge = FT_data[:,2]
# FT_pto2 = FT_data[:,3]
# FT_pt = FT_data[:,4]
# FT_ptcl2 = FT_data[:,5]
# FT_pth2_centre = FT_data[:,6]
# FT_pth2_edge = FT_data[:,7]
# plt.figure()
# plt.plot(FT_r,FT_data[:,1],linewidth=0.5,label='Pt centre')
# plt.plot(FT_r,FT_data[:,2],linewidth=0.5,label='Pt edge')
# plt.plot(FT_r,FT_data[:,3],linewidth=0.5,label='PtO_2')
# plt.plot(FT_r,FT_data[:,4],linewidth=0.5,label='Pt foil')
# plt.plot(FT_r,FT_data[:,5],linewidth=0.5,label='PtCl_2')
# plt.plot(FT_r,FT_data[:,6],linewidth=0.5,label='Pt H_2 centre')
# plt.plot(FT_r,FT_data[:,7],linewidth=0.5,label='Pt H_2 edge')
# plt.legend()
# plt.xlabel('K')
# plt.ylabel('chi')
# plt.savefig('/Users/Sophia/Downloads/drive-download-20190501T161928Z-001/Pt_ft',format='pdf')
# plt.show()


# #########################################
#
#     Compare Feff and excurve calculation
#
#############################################

#
# data_feff = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Simulation/Ge/feff0001.dat',skip_header=15)
# data_excurve = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Experiment_Analysis/Ge/c_ge12_pot_pha.dat',skip_header=5)
#
# plt.figure()
# plt.plot(data_feff[:,0],data_feff[:,2],label='FEFF')
# plt.plot(data_excurve[:,0],data_excurve[:,1],label='EXCURVE')
# plt.legend()
# plt.xlabel('k $\AA^{-1}$')
# plt.title('scattering amplitude')
# plt.show()
#
#
# plt.figure()
# plt.plot(data_feff[:,0],data_feff[:,1]+data_feff[:,3],label='FEFF ')
# plt.plot(data_excurve[:,0],data_excurve[:,3]+data_excurve[:,2],label='EXCURVE')
# plt.legend()
# plt.xlabel('k $\AA^{-1}$')
# plt.title('Phase')
# plt.show()




# #########################################
#
#     Fit gaussian to R space EXAFS peaks
#
#############################################

data = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Experiment_Analysis/CdS/CdS_M322_Aug18_pyspline.dat',skip_header=1)[:,(2,7)]
k = data[:,0]
exp = data[:,1]
r,amp,real,imag = calcfft(k,exp*k**3,kmin=2.4,kmax=16)
r_min = (np.abs(r - 1.6)).argmin()  #2.8
r_max = (np.abs(r - 2.9)).argmin() # 4.5
r1 = r[r_min:r_max]
amp1 = amp[r_min:r_max]

r_min = (np.abs(r - 2.8)).argmin()  #2.8
r_max = (np.abs(r - 4.3)).argmin() # 4.5
r2 = r[r_min:r_max]
amp2 = amp[r_min:r_max]

def model_gaussian(x,a,b,c):
    return a * np.exp(-(x-b)**2/2*c**2)


popt1,pcov = curve_fit(model_gaussian,r1,amp1)
popt2,pcov = curve_fit(model_gaussian,r2,amp2)

fit1 = model_gaussian(r,*popt1)
fit2 = model_gaussian(r,*popt2)

plt.figure()
plt.plot(r,amp,'k.')
plt.plot(r,fit1+fit2)
# plt.plot(r,fit2)
plt.xlim([0,5])
plt.show()

chiq = calcifft(fit1+fit2)

krange = np.linspace(0,k[-1],len(chiq))
print(krange)
print(len(chiq))
plt.figure()
plt.plot(k,exp*k**3,'k.')
plt.plot(krange,chiq)
plt.show()