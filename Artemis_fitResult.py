from matplotlib import pyplot as plt
import numpy as np
import glob
from matplotlib import cm
from xafs import *


def plot(sample):
    plt.figure()
    plt.subplot(1, 2, 1)
    filename_chi = glob.glob(
        '/Users/Sophia/ownCloud/PhD/Statistic Analysis/demeter/{:s}'.format(sample) + '*.k2')
    experiment = np.genfromtxt(filename_chi[0])
    experiment = experiment[:, (0, 1)]
    print(filename_chi)
    cmap = cm.get_cmap('winter')
    for file in filename_chi:
        data = np.genfromtxt(file)[:, 2]
        if 'OSCdCd' in file[-12:-3]:
            label = 'O S Cd Cd'
            i = 4
        elif 'OSCd' in file[-12:-3]:
            label = 'O S Cd'
            i = 3
        elif 'CdCd' in file[-12:-3]:
            label = 'S Cd'
            i = 2
        elif 'OS' in file[-12:-3]:
            label = 'O S'
            i = 1
        elif 'CdS' in file[-12:-3]:
            label = 'S'
            i = 0
        print(np.shape(experiment[:, 0]), np.shape(data))
        plt.plot(experiment[:, 0], experiment[:, 1] -
                 1.5 * i, 'ko', markersize=0.5)
        plt.plot(experiment[:, 0], data - 1.5 * i, label=label,
                 color=cmap((float(i) * 2 + 2) / 10), linewidth=0.5)
        plt.plot(experiment[:, 0], data - experiment[:, 1] - 7,
                 color=cmap((float(i) * 2 + 2) / 10), linewidth=0.5)
        r, amp = FT_chi(experiment[:, 0], data - experiment[:, 1], dx=2)
        # plt.plot(r*3,amp-5,color = cmap((float(i)*2+2)/10),linewidth=0.5)

    plt.xlim([2, 18])
    plt.ylim([-8, 1])
    plt.xlabel('k ($\AA$)')
    plt.ylabel('$k^2\chi(k)$')
    plt.legend(loc='upper right')

    plt.subplot(1, 2, 2)
    filename_r = glob.glob(
        '/Users/Sophia/ownCloud/PhD/Statistic Analysis/demeter/{:s}'.format(sample) + '*.rmag')
    experiment = np.genfromtxt(filename_r[0])
    experiment = experiment[:, (0, 1)]
    cmap = cm.get_cmap('winter')
    for file in filename_r:
        data = np.genfromtxt(file)[:, 2]
        if 'OSCdCd' in file[-12:-3]:
            label = 'O S Cd Cd'
            i = 4
        elif 'OSCd' in file[-12:-3]:
            label = 'O S Cd'
            i = 3
        elif 'CdCd' in file[-12:-3]:
            label = 'S Cd'
            i = 2
        elif 'OS' in file[-12:-3]:
            label = 'O S'
            i = 1
        elif 'CdS' in file[-12:-3]:
            label = 'S'
            i = 0
        plt.plot(experiment[:, 0], experiment[:, 1] -
                 1.5 * i, 'ko', markersize=0.5)
        plt.plot(experiment[:, 0], data - 1.5 * i, label=label,
                 color=cmap((float(i) * 2 + 2) / 10), linewidth=0.5)

    plt.xlim([0, 5])
    plt.xlabel('Radial distance ($\AA$)')
    plt.yticks([])
    plt.legend(loc='upper right')
    plt.savefig(
        '/Users/Sophia/ownCloud/PhD/Statistic Analysis/figure/artemis_fit_{:s}.pdf'.format(sample), format='pdf')


plot('M311')
plot('M322')

plt.figure()
plt.subplot(1, 2, 1)
M311_chi_OS = np.genfromtxt(
    '/Users/Sophia/ownCloud/PhD/Statistic Analysis/demeter/M311_CdOS.k2')
M311_chi_OSCd = np.genfromtxt(
    '/Users/Sophia/ownCloud/PhD/Statistic Analysis/demeter/M311_CdOSCd.k2')
plt.plot(M311_chi_OS[:, 0], M311_chi_OS[:, 1],
         'ko', markersize=0.5, label='M311')
plt.plot(M311_chi_OS[:, 0], M311_chi_OS[:, 2], linewidth=0.5, label='O S')
plt.plot(M311_chi_OSCd[:, 0], M311_chi_OSCd[:, 2],
         linewidth=0.5, label='O S Cd')

plt.xlim([0, 18])
plt.xlabel('k ($\AA$)')
plt.ylabel('$k^2\chi(k)$')
plt.legend(loc='upper right')


plt.subplot(1, 2, 2)
M322_chi_OS = np.genfromtxt(
    '/Users/Sophia/ownCloud/PhD/Statistic Analysis/demeter/M322_CdOS.k2')
M322_chi_OSCdCd = np.genfromtxt(
    '/Users/Sophia/ownCloud/PhD/Statistic Analysis/demeter/M322_CdOSCdCd.k2')
plt.plot(M322_chi_OS[:, 0], M322_chi_OS[:, 1],
         'ko', markersize=0.5, label='M322')
plt.plot(M322_chi_OS[:, 0], M322_chi_OS[:, 2], linewidth=0.5, label='O S')
plt.plot(M322_chi_OSCdCd[:, 0], M322_chi_OSCdCd[:, 2],
         linewidth=0.5, label='O S Cd Cd')

plt.xlim([0, 18])
plt.xlabel('k ($\AA$)')
plt.yticks([])
plt.legend(loc='upper right')

plt.savefig(
    '/Users/Sophia/ownCloud/PhD/Statistic Analysis/figure/artemis_fit_compare.pdf', format='pdf')
