from matplotlib import pyplot as plt
import numpy as np

def analyze(num_atoms, model):
    betas = np.array([100, 10, 5, 4, 3, 2, 1, 0.5, 0.25, 0.1, 0.05, 0.01, 0.005, 0.001])

    av_occ = []
    fig, ax = plt.subplots()
    for b in betas:
        #print(b)
        en = np.load("./run_dat/en_trace_{}_{}_{}.npy".format(int(num_atoms), b, model))
        occ = np.load("./run_dat/occ_trace_{}_{}_{}.npy".format(int(num_atoms), b, model))

        acc = np.load("./run_dat/acc_rate_{}_{}_{}.npy".format(int(num_atoms), b, model))
        #print(occ)
        av_occ.append(np.average(occ[50:]))
        
        #if(b == 1 or b == 0.005):
        #    ax.plot(occ, label=b)

    fig, ax = plt.subplots()
    ax.plot(1/(betas*num_atoms), av_occ)
    ax.set(xlim=(0, 5))
    ax.legend(loc='best')
    plt.show()

if __name__ == '__main__':
    num_atoms = 50
    model = 'free'
    analyze(num_atoms, model)
    pass
