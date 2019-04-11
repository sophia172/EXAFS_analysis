from xafs.__init__ import *

def plot_FEFF_average(self,model):
    filename = glob.glob('/Users/Sophia/ownCloud/PhD/Matlab figure/CdS/EXAFS/EXAFSsimulation_{:s}*'.format(model))
    for i in filename:
        data = np.genfromtxt(i)
        plt.plot(data[:,2],data[:,5] * data[:,2]**2,label = i[67:-8])
    plt.plot(self.k, self.experiment * self.k**2,'k.',markersize=2,label=self.sample)
    plt.legend()
    plt.xlabel('k')
    plt.ylabel('Chi')
    plt.show()

def plot_xmus(self, folder):
    filelist = glob.glob('/Users/Sophia/ownCloud/PhD/Simulation/CdS/FEFF/EXAFS/Cd/{:s}/*xmu.dat'.format(folder))
    for file in filelist:
        data = np.genfromtxt(file)
        plt.plot(data[:, 2], data[:, 5] * data[:, 2] ** 2)
    plt.plot(self.k, self.experiment * self.k ** 2, 'k.', markersize=2, label=self.sample)
    plt.legend()
    plt.xlabel('k')
    plt.ylabel('Chi')
    plt.title('FEFF InP model (SCd shell)')
    plt.show()
    return

def compare_output(bond,column,shell,i):
    data = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Statistic Analysis/data/feff000{:s}_{:s}.dat'.format(shell,bond),skip_header=15)
    k = data[:,0]
    if column == 'phase':
        y = data[:,1]+data[:,3]
    elif column == 'Feff':
        y = data[:,2]
    elif column == 'lambda':
        y = data[:,5]
    plt.plot(k,y,label = 'Cd-S bond length {:.0f}% of 2.6164'.format(float(bond)/2.6164*100),color = (i*20/255,(i*10+50)/255,200/255))
    plt.xlabel('k')
    plt.ylabel('{:s}'.format(column))
    plt.title('{:s} in Shell {:s}'.format(column,shell))

    print('plot {:s}{:s} is ready to show'.format(column,bond))
    return

def plot_feff(num_file,shell,variance):

    if variance == 2:
        gap = 0.052328
    elif variance == 10:
        gap = 0.26164
    plt.figure(figsize=(6,8))
    plt.subplot(3,1,1)
    for i in range(1,num_file):
        bond = 2.6164 + gap * (i-6)    # for 2% change, it is 0.052328, for 10% change, it is 0.26164
        compare_output('{:.4f}'.format(bond),'phase','{:d}'.format(shell),i)

    plt.subplot(3,1,2)
    for i in range(1,num_file):
        bond = 2.6164 + gap * (i-6)
        compare_output('{:.4f}'.format(bond),'Feff','{:d}'.format(shell),i)

    plt.subplot(3,1,3)
    for i in range(1,num_file):
        bond = 2.6164+ gap * (i-6)
        compare_output('{:.4f}'.format(bond),'lambda','{:d}'.format(shell),i)

    plt.figlegend()
    plt.show()
    return

