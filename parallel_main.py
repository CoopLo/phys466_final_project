import numpy as np
from time import clock
import scipy
from scipy import integrate as integrate
from matplotlib import pyplot as plt
from multiprocessing import Pool


# Initialize lattice by samlping bose-einstein distribution
# In chemical potential = 0 limit?
def be_lattice(size, beta):
    lattice = np.zeros((size, size, size, 3))
    # inverse bose-einstein distribution

    def be_dist(x):
        return x*(x*(x+3)+1)/(2*(np.exp(-beta*x))+1)

    def boltz_dist(x):
        return x*np.exp(-beta*x)

    #print(be_dist(1))
    #e = integrate.quad(be_dist, 0, np.inf)
    #print(e)

    e = integrate.quad(boltz_dist, 0, np.inf)
    if(e[0] < 1):
        for i in range(lattice.shape[0]):
            for j in range(lattice.shape[1]):
               for k in range(lattice.shape[2]):
                   for l in range(lattice.shape[3]):
                       if(np.random.random() < e[0]):
                           lattice[i][j][k][l] = 1
    else:
        lattice.fill(int(e[0]/3))
                
    #return np.full((size, size, size, 3), 10)
    return np.array(lattice)


def get_energy(lattice, mass):
   ''' 
       parameters:
         lattice: the lattice
         mass: mass of each particle
   '''
   return np.sum(lattice+0.5)
    

###
# parameters:
#    lattice: the lattice we're working with
#    site: the site on the lattice being changed
#    e_coord: the energy coordinate being changed e.g. (0,0,1) for z, or (1,0,0) for x
#
# returns:
#    Difference in energy from changing specifiec site
###
def energy_diff(lattice, site, e_coord):
    old_en = np.sum((lattice[site[0]][site[1]][site[2]]))
    new_en = np.sum(lattice[site[0]][site[1]][site[2]]+e_coord)
    return old_en - new_en


#def mc(nsweeps, betas, i, idx):
def mc(size):

    # Temperatures over which to simulate
    betas = np.array([50., 40., 30., 20., 10.,
                      5., 4., 3., 2., 1., 0.5, 0.25, 0.1])

    # Number of sweeps varies for each system size
    nsweeps = [10000, 20000, 20000, 50000, 100000, 100000, 5000000]

    # Fraction of data to throw out for equilibrium
    eq_frac = [0.4, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5]

    # Sizes of system, used here for index
    sizes = [2, 3, 4, 5, 7, 10, 22]
    i = sizes.index(size)

    # Because I guess I use this somewhere else
    mass = 48

    runs = 20
    nsweep = nsweeps[i]
    num_atoms = size**3

    # Do for each temperature
    for idx, beta in enumerate(betas):
        av_en = np.zeros((len(betas), runs))
        av_gs_occ = np.zeros((len(betas), runs))

        # More than one run so we can average later
        for t in range(runs):
            lattice = be_lattice(size, betas[idx])
            energy = []
            gs_occ = []

            # Number of sweeps
            for isweep in range(nsweeps[i]):

                # Each sweep could change every spin
                for k in range(num_atoms):

                    # Select site to change, coordinate, and up or down
                    site = np.random.randint(0, size, (1, 3))[0]
                    coord = np.random.randint(0, 3)
                    e_coord = [1,0,0] if coord==0 else [0,1,0] if coord==1 else [0,0,1] 
                    change = np.random.randint(0, 2)
                    change = -1 if change==0 else 1

                    # change in energy from site change
                    delta_e = energy_diff(lattice, site, np.array(e_coord)*change)
                    
                    # transition probability
                    prob = 1/2 * 1./np.cosh(delta_e*beta/2) * np.exp(change*delta_e*beta/2)

                    # Transition if you should
                    if(prob > np.random.random()):
                        lattice[site[0]][site[1]][site[2]][coord] += change

                    # calculation of ground state occupancy
                    gs_occ.append(gsoccupancy/num_atoms)
                    
                    # calculation of total energy
                    energy.append(get_energy(lattice, mass))

            # Average energy and gs occupancy for temperature
            av_en[idx][t] = np.average(energy[int(eq_frac[i]*nsweeps[i]):])
            av_gs_occ[idx][t] = np.average(gs_occ[int(eq_frac[i]*nsweeps[i]):])
            
        # Averageing and making plot
        fig, ax = plt.subplots()
        av_en = np.average(av_en, axis=1)
        av_gs_occ = np.average(av_gs_occ, axis=1)
        ax.plot(1/betas, av_gs_occ)
        plt.savefig("./graphs/{}_parallel.png".format(size))


if __name__ == '__main__':

    # Size of simulation
    sizes = [2, 3, 4, 5, 7, 10, 22]
    with Pool(len(sizes)) as p:
        print(p.map(mc, sizes))
                
