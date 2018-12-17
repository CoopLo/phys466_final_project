from matplotlib import pyplot as plt
import numpy as np

def analyze(num_atoms, model):
    #betas = np.array([100, 10, 5, 4, 3, 2, 1, 0.5, 0.25, 0.1, 0.05, 0.01, 0.005, 0.001])
    betas = np.array([50, 10, 5., 4., 3., 2., 1., 0.75, 0.7, 0.6, 0.5, 0.4, 0.3, .25, 0.1])
    #betas = np.array([0.001])

    #fig, ax = plt.subplots()
    fig, ax = plt.subplots()
    #fig1, ax1 = plt.subplots(2)
    print(num_atoms)
    for n in num_atoms:
        print(n)
        av_occ = []
        for b in betas:
            #print(b)
            en = np.load("./good_run_dat/en_trace_{}_{}_{}.npy".format(int(n), b,
                                                                       model))
            occ = np.load("./good_run_dat/occ_trace_{}_{}_{}.npy".format(int(n), b,
                                                                       model))
            acc = np.load("./good_run_dat/acc_rate_{}_{}_{}.npy".format(int(n), b,
                                                                       model))
            #print(b, acc)
            av_occ.append(np.average(occ[500:]))
            #if(n==100):
            #    ax1[0].plot(en, label=1/b, linewidth=0.3)
            #    ax1[1].plot(occ, label=1/b, linewidth=0.3)
            #    print(av_occ[0]*n)
            #    #plt.show()
            #    #exit(1)
            #print(np.average(occ[50:]))

        #print(av_occ)
        ax.plot(1/(betas * (n**(1/3))), av_occ, label="{} Particles".format(n))
    ax.set(xlim=(0, 1.5), xlabel=r"$T/N^{1/3}$", ylabel=r"$\left<N_0\right>/N$",
           title="Ground State Occupancy vs. Temperature")
    ax.legend(loc='best')
    #ax1[0].legend(loc='best')
    #ax1[1].legend(loc='best')
    plt.show()

if __name__ == '__main__':
    num_atoms = [5, 10, 50, 100, 200, 500]
    model = 'harmonic'
    analyze(num_atoms, model)

