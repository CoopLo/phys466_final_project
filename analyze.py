from matplotlib import pyplot as plt
import numpy as np
import tools

def canonical_analyze(num_atoms, model):
    #betas = np.array([100, 10, 5, 4, 3, 2, 1, 0.5, 0.25, 0.1, 0.05, 0.01, 0.005, 0.001])
    if(model == "harmonic"):
        betas = np.array([50, 10, 5., 4., 3., 2., 1., 0.75, 0.7, 0.6, 0.5, 0.4, 0.3, .25, 0.1])

    elif(model == "free"):
        betas = np.array([50, 10, 5., 4., 3., 2., 1., 0.75, 0.7, 0.6, 0.5, 0.4, 0.3, .25, 0.1,
                          0.08, 0.05, 0.03, 0.01, 0.008, 0.005])
    else:
        print("INVALID MODEL SELECTED")
        exit(1)

    fig, ax = plt.subplots()
    #fig1, ax1 = plt.subplots(2)
    print(num_atoms)
    colors = ['#a6cee3', '#b2df8a', '#e31a1c', '#cab2d6', '#b15928', '#e6ab02']
    marker = ['s', 'p', 'P', 'o', 'd', "D"]
    for idx, n in enumerate(num_atoms):
        print(n)
        av_occ = []
        av_std = []
        for b in betas:
            en = np.load("./bad_run_dat/en_trace_{}_{}_{}.npy".format(int(n), b,
                                                                       model))
            occ = np.load("./bad_run_dat/occ_trace_{}_{}_{}.npy".format(int(n), b,
                                                                       model))
            acc = np.load("./bad_run_dat/acc_rate_{}_{}_{}.npy".format(int(n), b,
                                                                       model))
            av_occ.append(np.average(occ[500:]))
            #av_std.append(np.std(occ[500:]))
            av_std.append(tools.std_error(occ[500:]))
            #if(n==100):
            #    ax1[0].plot(en, label=1/b, linewidth=0.3)
            #    ax1[1].plot(occ, label=1/b, linewidth=0.3)
                #print(av_occ[0]*n)
                #plt.show()
                #exit(1)
            #print(np.average(occ[50:]))

        #print(av_occ)

        #ax.plot(1/(betas * (n**(1/3))), av_occ, label="{} Particles".format(n), linewidth=0.7)
        ax.scatter(1/(betas * (n**(1/3))), av_occ, marker=marker[idx], color=colors[idx],
                   edgecolor='k', label='{} Particles'.format(n))
        ax.errorbar(1/(betas * (n**(1/3))), av_occ, av_std, capsize=7, color=colors[idx])
    ax.set(xlabel=r"$T/N^{1/3}$", ylabel=r"$\left<N_0\right>/N$",
           title="Boltzmanon Ground State Occupancy vs. Temperature", xlim=(0, 1.5))
    ax.legend(loc='best')
    #ax1[0].legend(loc='best')
    #ax1[1].legend(loc='best')
    plt.show()


def canonical_dynamic(num_atoms, model):
    #betas = np.array([100, 10, 5, 4, 3, 2, 1, 0.5, 0.25, 0.1, 0.05, 0.01, 0.005, 0.001])
    if(model == "harmonic"):
        betas = np.array([50, 10, 5., 4., 3., 2., 1., 0.75, 0.7, 0.6, 0.5, 0.4, 0.3, .25, 0.1])

    elif(model == "free"):
        betas = np.array([50, 10, 5., 4., 3., 2., 1., 0.75, 0.7, 0.6, 0.5, 0.4, 0.3, .25, 0.1,
                          0.08, 0.05, 0.03, 0.01, 0.008, 0.005])
    else:
        print("INVALID MODEL SELECTED")
        exit(1)

    #fig, ax = plt.subplots()
    fig, ax = plt.subplots()
    fig1, ax1 = plt.subplots(2)
    print(num_atoms)
    for n in num_atoms:
        print(n)
        av_occ = []
        #print(b)
        en = np.load("./dynamic_run_dat/en_trace_{}_{}.npy".format(int(n),
                                                                   model))
        occ = np.load("./dynamic_run_dat/occ_trace_{}_{}.npy".format(int(n),
                                                                   model))
        acc = np.load("./dynamic_run_dat/acc_rate_{}_{}.npy".format(int(n),
                                                                   model))
        #av_occ.append(np.average(occ[500:]))
        #if(n==100):
        #    ax1[0].plot(en, label=1/b, linewidth=0.3)
        #    ax1[1].plot(occ, label=1/b, linewidth=0.3)
            #print(av_occ[0]*n)
            #plt.show()
            #exit(1)
        #print(np.average(occ[50:]))

        #print(av_occ)
        ax1[0].plot(en, label="{} Particles".format(n), linewidth=0.7)
        ax1[1].plot(occ, label="{} Particles".format(n), linewidth=0.7)
    ax.set(xlabel=r"$T/N^{1/3}$", ylabel=r"$\left<N_0\right>/N$",
           title="Ground State Occupancy vs. Temperature")#, xlim=(0, 1.5))
    ax.legend(loc='best')
    #ax1[0].legend(loc='best')
    #ax1[1].legend(loc='best')
    plt.show()


def grand_canonical_analyze(num_atoms, model):
    #betas = np.array([100, 10, 5, 4, 3, 2, 1, 0.5, 0.25, 0.1, 0.05, 0.01, 0.005, 0.001])
    if(model == "harmonic"):
        betas = np.array([50, 10, 5., 4., 3., 2., 1., 0.75, 0.7, 0.6, 0.5, 0.4, 0.3, .25, 0.1])

    elif(model == "free"):
        betas = np.array([50, 10, 5., 4., 3., 2., 1., 0.75, 0.7, 0.6, 0.5, 0.4, 0.3, .25, 0.1,
                          0.08, 0.05, 0.03, 0.01, 0.008, 0.005])
    else:
        print("INVALID MODEL SELECTED")
        exit(1)

    #fig, ax = plt.subplots()
    fig, ax = plt.subplots()
    fig1, ax1 = plt.subplots(2)
    print(num_atoms)
    for n in num_atoms:
        print(n)
        av_occ = []
        for b in betas:
            #print(b)
            en = np.load("./gc_run_dat/en_trace_{}_{}_{}.npy".format(int(n), b,
                                                                       model))
            occ = np.load("./gc_run_dat/occ_trace_{}_{}_{}.npy".format(int(n), b,
                                                                       model))
            acc = np.load("./gc_run_dat/acc_rate_{}_{}_{}.npy".format(int(n), b,
                                                                       model))
            #print(b, acc)
            av_occ.append(np.average(occ[500:]))
            if(n==100):
                ax1[0].plot(en, label=1/b, linewidth=0.3)
                ax1[1].plot(occ, label=1/b, linewidth=0.3)
                print(av_occ[0]*n)
                #plt.show()
                #exit(1)
            #print(np.average(occ[50:]))

        #print(av_occ)
        ax.plot(1/(betas * (n**(1/3))), av_occ, label="{} Particles".format(n))
    ax.set(xlabel=r"$T/N^{1/3}$", ylabel=r"$\left<N_0\right>/N$",
           title="Ground State Occupancy vs. Temperature", xlim=(0, 1.5))
    ax.legend(loc='best')
    #ax1[0].legend(loc='best')
    #ax1[1].legend(loc='best')
    plt.show()


if __name__ == '__main__':
    num_atoms = [5, 10, 50, 100, 200, 500]
    model = 'harmonic'
    #canonical_analyze(num_atoms, model)
    canonical_dynamic(num_atoms, model)
    #grand_canonical_analyze(num_atoms, model)

