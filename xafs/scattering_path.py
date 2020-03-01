import os
import glob
import numpy as np
from scipy.signal import hilbert
import matplotlib.pyplot as plt
from EXAFS_analysis.xafs import interpolation



class scattering_path():
    def __init__(self,linked_atoms):
        self.linked_atoms = linked_atoms
        try:
            self.FEFF_path = os.path.join(os.getcwd(), 'source', 'feff_' + linked_atoms + '.dat')
            self.FEFF = np.genfromtxt(self.FEFF_path)
            print(os.path.basename(self.FEFF_path), '      Found \n')
        except:
            print('FEFF file not found for %s, this step is required for mean free path ... \n' %linked_atoms)
            print(self.FEFF_path, '\n')
        try:
            self.exp_path = os.path.join(os.getcwd(), 'source', linked_atoms + '.chiq')
            self.exp_chi = np.genfromtxt(self.exp_path)[:,(0,1)]
            print(os.path.basename(self.exp_path), '      Found \n')
        except:
            print('reference file not found for single scattering path %s, this step is required for amplitude and phase extraction ....\n' %linked_atoms)
            print(self.exp_path,'\n')
        self.k = np.linspace(2,17,60000)

        self.N = 4
        self.R = 2.5

    def set_fit_range(self,start=2,end=17):
        print('fit range reset ...')
        self.k = np.linspace(start-1, end +0.5, int((end-start) * 4000))


    def set_params(self, N=4, R=2.5):
        print('parameter reset ...')
        self.N = N
        self.R = R

    def FEFF_emfp(self):
        self.emfp = self.FEFF[:,(0,5)]
        self.emfp = interpolation(self.emfp, self.k)
        print('electron mean free path extracted from FEFF package for %s', self.linked_atoms)

    def cal_exp_params(self):
        exp_chi = interpolation(self.exp_chi,self.k) / self.k**2
        amp = np.abs(hilbert(exp_chi))
        self.exp_phase = np.unwrap(np.angle(hilbert(exp_chi)))+np.pi/2 - 2 * self.k * self.R
        print('Phase calculated for %s' %self.linked_atoms)
        self.FEFF_emfp()
        self.exp_feff = amp/np.exp(-2 * self.R / self.emfp) / self.N* self.k * self.R**2
        print('Effective backscattering amplitude calculated for %s' %self.linked_atoms)

    def cal_FEFF_params(self):
        FEFF = np.genfromtxt(self.FEFF_path)
        self.FEFF_feff = np.vstack((FEFF[:,0],FEFF[:,2])).T
        self.FEFF_feff = interpolation(self.FEFF_feff, self.k)
        self.FEFF_phase = np.vstack((FEFF[:, 0], FEFF[:, 1]+ FEFF[:, 3])).T
        self.FEFF_phase = interpolation(self.FEFF_phase, self.k)
        print('FEFF package Feff and phase calculated for %s' %self.linked_atoms)

    def plot(self):

        self.cal_exp_params()
        self.cal_FEFF_params()

        plt.figure()
        plt.subplot(2,1,1)
        plt.plot(self.k, self.FEFF_feff, label='FEFF calculated Feff')
        plt.plot(self.k, self.exp_feff, label='experiment Feff  including MSRD')
        plt.title('foundamental experiement Feff extraction')
        plt.legend()
        plt.tight_layout()

        plt.subplot(2, 1, 2)
        plt.plot(self.k, self.FEFF_phase, label='FEFF calculated phase')
        plt.plot(self.k, self.exp_phase, label='experiment phase')
        plt.title('foundamental experiement phase extraction')
        plt.legend()
        plt.tight_layout()

        fig_path = os.path.join(os.getcwd() + '/result/')
        plt.savefig(fig_path +'extracted_{:s}.pdf'.format(self.linked_atoms),
                    format='pdf')
        plt.close()
        data = np.vstack((self.k,self.exp_feff,self.exp_phase)).transpose()
        np.savetxt(fig_path +'extracted_{:s}.dat'.format(self.linked_atoms),data,
                   header='k, Feff including MSRD, phase')



