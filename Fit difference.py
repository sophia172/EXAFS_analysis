import numpy as np
from matplotlib import pyplot as plt
import glob
from xafs.FT import *
from matplotlib import cm
# data_path = '/Users/sophia/ownCloud/PhD/Statistic Analysis/CdS/'
# figure_path = '/Users/sophia/ownCloud/PhD/Statistic Analysis/CdS/figure/'
# filelist = glob.glob(data_path+'athena*_fit.dat')
cmap = cm.get_cmap('tab10')
data_dic={}
filelist = ['bulk','R343','M311','M322']
def sample_file():
    data_path = '/Users/sophia/owncloud/PhD/Experiment_Analysis/CdS/'
    figure_path = '/Users/sophia/owncloud/PhD/Experiment_Analysis/CdS/'
    data = np.genfromtxt(data_path+'CdS_XAFS_CdK_chi_2018.chi')
    data_bulk = np.genfromtxt(data_path+'CdS_bulk_Aug18_athena_12sh.chiq')
    data_dic['bulk'] = data_bulk[:,(0,1)]
    data_dic['R343'] = data[:,(0,4)]
    data_dic['M311'] = data[:,(0,2)]
    data_dic['M322']  = data[:,(0,3)]
    return figure_path


def plot_chi(file,i):
    data = data_dic[file]
    if file == 'bulk':
        k_weight =0
    else:
        k_weight =2
    plt.plot(data[:,0],data[:,1]*data[:,0]**k_weight,color = cmap(i*2/10),linewidth=1,label=file)
    return data[:,0],data[:,1]*data[:,0]**k_weight

def plot_FT(file,i):
    data = data_dic[file]
    if file == 'bulk':
        k_weight =0
    else:
        k_weight =2
    r,amp,real,imag = calcfft(data[:,0],data[:,1]*data[:,0]**k_weight,kmin=2.0,kmax=17)
    plt.plot(r,amp,color = cmap(i*2/10),linewidth = 1,label=file)
    return r,amp


def plot_fit_diff(file_list):
    plt.figure()

    plt.subplot(2, 1, 1)
    for file in filelist:
        plot_chi(file[57:-8])
    plt.xlim([1, 17])
    plt.ylim([-0.2, 0.2])
    plt.xlabel('k ($\AA$)')
    plt.ylabel('$k^2\chi(k)$')
    plt.legend(loc='upper right')

    plt.subplot(2, 1, 2)
    for file in filelist:
        plot_FT(file[57:-8])
    plt.xlim([0, 5])
    plt.xlabel('Radial distance ($\AA$)')
    plt.yticks([])
    plt.legend(loc='upper right')
    plt.tight_layout()

def plot_chi_r(filelist):
    fig = plt.figure(figsize=(5,3.5))
    fig.patch.set_alpha(0)



    ax = fig.add_subplot(2,1,1)
    for i,file in enumerate(filelist):
        plot_chi(file,i)
    ax.spines['bottom'].set_color('white')
    ax.spines['top'].set_color('white')
    ax.spines['right'].set_color('white')
    ax.spines['left'].set_color('white')

    ax.tick_params(axis='x', colors='white')
    ax.tick_params(axis='y', colors='white')
    ax.xaxis.label.set_color('white')
    ax.yaxis.label.set_color('white')

    plt.xlabel('k ($\AA^{-1}$)')
    plt.ylabel('$k^2\chi(k)$')
    legend = plt.legend(loc='upper right',framealpha=0)
    plt.setp(legend.get_texts(), color='w')
    ax.patch.set_alpha(0)
    plt.tight_layout()

    ax = fig.add_subplot(2,1,2)
    for i, file in enumerate(filelist):
        plot_FT(file, i)
    plt.xlim([0, 7])

    ax.spines['bottom'].set_color('white')
    ax.spines['top'].set_color('white')
    ax.spines['right'].set_color('white')
    ax.spines['left'].set_color('white')

    ax.tick_params(axis='x', colors='white')
    ax.tick_params(axis='y', colors='white')
    ax.xaxis.label.set_color('white')
    ax.yaxis.label.set_color('white')

    plt.xlabel('Radial distance ($\AA$)')
    plt.yticks([])
    legend = plt.legend(loc='upper right',framealpha=0)
    plt.setp(legend.get_texts(), color='w')
    # legend.get_frame().set_framealpha(0)
    ax.patch.set_alpha(0)

    plt.tight_layout()
    fig.savefig(figure_path+'FT_CdS_2018.pdf',format='pdf')
    plt.show()
    plt.close()

dic = {}
def plot_rs(filelist):
    plt.figure(figsize=(5,3.5))
    # plt.subplot(2,1,1)
    # for i,file in enumerate(filelist):
    #     plot_chi(file,i)
    # plt.xlabel('k ($\AA^{-1}$)')
    # plt.ylabel('$k^2\chi(k)$')
    # plt.legend(loc='upper right')
    # plt.tight_layout()

    # plt.subplot(2,1,2)
    for i, file in enumerate(filelist):
        dic['r_'+file],dic['amp_'+file] = plot_FT(file, i)

    # plt.plot(dic['r_M311'],dic['amp_M322'] - dic['amp_M311'],'k-',label='residule of M311 and M322',linewidth=0.5)
    plt.xlim([0, 7])
    plt.xlabel('Radial distance ($\AA$)')
    plt.yticks([])
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig(figure_path+'FT_CdS_2018.pdf',format='pdf')
    plt.show()
    plt.close()

def plot_chis(filelist):
    plt.figure(figsize=(5,3.5))
    # plt.subplot(2,1,1)
    for i,file in enumerate(['M311','M322']):
        print(file)
        dic['k_'+file],dic['chik2_'+file] = plot_chi(file,i+2)
    plt.plot(dic['k_M311'],dic['chik2_M322'] - dic['chik2_M311'],'k-',label='residual  of \nM311&M322',linewidth=0.5)


    plt.xlabel('k ($\AA^{-1}$)')
    plt.ylabel('$k^2\chi(k)$')
    plt.legend(loc='upper right',fontsize=8)
    plt.tight_layout()

    # plt.subplot(2,1,2)
    # for i, file in enumerate(filelist):
    #     plot_FT(file, i)
    # plt.xlim([0, 7])
    # plt.xlabel('Radial distance ($\AA$)')
    # plt.yticks([])
    # plt.legend(loc='upper right')
    # plt.tight_layout()
    plt.savefig(figure_path+'chi_CdS_2018.pdf',format='pdf')
    plt.show()
    plt.close()

def plot_UV():

    data_path = '/Users/sophia/owncloud/PhD/Experiment_Analysis/PL_UV/'
    figure_path = '/Users/sophia/owncloud/PhD/Experiment_Analysis/PL_UV/'

    data_dic['bulk'] = np.genfromtxt(data_path + 'CdSBulk.txt')[:, ]
    data_dic['R343'] = np.genfromtxt(data_path + 'CdSBulk.txt')[:, ]
    data_dic['M311'] = np.genfromtxt(data_path + 'CdS311M.txt')[:, ]
    data_dic['M322'] = np.genfromtxt(data_path + 'CdS322M.txt')[:,]

    plt.figure(1)
    for file in filelist:
        plt.plot( data_dic[file][:,0], data_dic[file][:,1], label = file)

    plt.show()

# plot_UV()


def plot_fits():
    data_path = '/Users/sophia/ownCloud/PhD/Statistic Analysis/CdS/result/CdO_CdS_CdCd/'
    fig_path = '/Users/sophia/ownCloud/PhD/Statistic Analysis/CdS/figure/'
    data_dic['bulk'] = np.genfromtxt(data_path + 'athena_bulk_12sh_fit.dat')[:,]
    data_dic['R343'] = np.genfromtxt(data_path + 'athena_R_fit.dat')[:, ]
    data_dic['M311'] = np.genfromtxt(data_path + 'athena_M311_fit.dat')[:, ]
    data_dic['M322'] = np.genfromtxt(data_path + 'athena_M322_fit.dat')[:, ]

    fig = plt.figure(figsize=(5,4))
    fig.patch.set_alpha(0)
    ax  = fig.add_subplot(1, 1, 1)
    ax.patch.set_alpha(0)
    for i,file in enumerate(filelist):
        plt.plot(data_dic[file][:,0], data_dic[file][:,1] * data_dic[file][:,0]**2 - i*1.5,'w.',markersize = 1)
        plt.plot(data_dic[file][:,0], data_dic[file][:,2] * data_dic[file][:,0]**2 - i*1.5, color=cmap(i*2/10), linewidth = 1, label=file)
        plt.plot(data_dic[file][:, 0], data_dic[file][:, 3] * data_dic[file][:,0]**2 - i *1.5,color=cmap(i*2/10), linewidth=1)
        plt.text(14, -i*1.5+0.3, file ,
                verticalalignment='bottom', horizontalalignment='left',
                color=cmap(i*2/10), fontsize=12)
    plt.xlabel('k ($\AA^{-1}$)')
    plt.ylabel('$k^2\chi(k)$')

    ax.spines['bottom'].set_color('white')
    ax.spines['top'].set_color('white')
    ax.spines['right'].set_color('white')
    ax.spines['left'].set_color('white')

    ax.tick_params(axis='x', colors='white')
    ax.tick_params(axis='y', colors='white')
    ax.xaxis.label.set_color('white')
    ax.yaxis.label.set_color('white')
    plt.tight_layout()
    plt.savefig(fig_path+'fit.pdf',format='pdf')
    plt.show()
    plt.close()

figure_path = sample_file()
# plot_chi_r(filelist)
# plot_rs(filelist[1:])
plot_chis(filelist)
# plot_fits()

