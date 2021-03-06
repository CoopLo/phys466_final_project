import numpy as np
from time import clock
import scipy
from matplotlib import pyplot as plt

def initialize_lattice(num_atoms, beta, lim):
    ''' oparameters:
            num_atoms: number of atoms we're simulating with. Should be cubic number.
            beta: temperature parameter, hbar * omega / T
            lim: highest energy state occupied
        returns:
               lattice: (num_atoms**1/3, num_atoms**1/3, num_atoms**1/3, 3) array with energy levels distributed according to BE distribution
      NOTE: Assumes energy levels of each particle same as that of quantum harmonic oscillator hbar * omega (n + 1/2)         
      '''
    
    nx = int(round(num_atoms ** 1/3))
    lattice = np.zeros([nx, nx, nx, 3], dtype=float)
    lat = np.random.rand(nx, nx, nx, 3)   # generates (nx, nx, nx, 3) array of random integers
    prob = []
    for i in range(lim):
        compare = 1 / (np.exp(beta * (i+0.5))-1) # calculates probability of being in energy level i
        prob.append(compare)
    prob = np.cumsum(prob/np.sum(prob))    # calculates cumulative probability for energy levels
    for i in range(lim - 1):   # places lattice sites in higher energy levels given probabilities
        sel = (lat < prob[i+1]) & (lat > prob[i])
        lattice[sel] = i
    return lattice          # returns (nx, nx, nx, 3) array of energy levels 


# Initialize lattice by samlping bose-einstein distribution
# In chemical potential = 0 limit?
def be_lattice(size, beta):
    
    lattice = np.zeros((size, size, size, 3))
    # inverse bose-einstein distribution

    def be_dist(x):
        return (x*(x+3)+1)/(2*(np.exp(-beta*x))+1)

    #return np.full((size, size, size, 3), 10)
    return np.array(lattice)



def get_energy(lattice, mass):
   ''' 
       parameters:
         lattice: the lattice
         mass: mass of each particle
   '''
   return mass/2 * np.sum(np.square(lattice))
    

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


if __name__ == '__main__':

    # Initial conditions in reduced units
    #betas = np.array([0.01, 0.05, 0.1, 0.5, 1.0, 5.0, 10.0, 1000.0])
    betas = np.array([10000., 1000., 100., 75., 50., 40., 30., 20., 10., 9., 8., 7., 6.,
                      5., 4., 3., 2., 1.])
    #betas = np.array([1])
    mass = 48.
    lim = 10
    nsweep = 20000
    lat_dat = 'lat.dat'
    scalar_dat = 'scalar.dat'
    nconf = 10
    ntherm = 1
    therm_fmt = '%8d %14.8f %14.8f %14.8f\n'
    therm_header = '# isweep  m2  mavg  pe\n'

    runs = 10
    av_en = np.zeros((len(betas), runs))
    av_gs_occ = np.zeros((len(betas), runs))

    sizes = [2, 3, 4, 5, 7, 10]
    for size in sizes:
        num_atoms = size**3
        for t in range(runs):
            for idx, beta in enumerate(betas):
                energy = []
                gs_occ = []
                lattice = be_lattice(size, beta)
                msg = 'starting %s sweeps using temperature %s and number of particles %s' % \
                  (nsweep, beta, num_atoms)
                print(msg)
                for isweep in range(nsweep):
                    if (isweep % ntherm == 0):

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
                        if(prob > np.random.random()):
                            lattice[site[0]][site[1]][site[2]][coord] += change

                        # calculation of ground state occupancy
                        #gsoccupancy = sum(np.shape(~lattice.any(axis=3))) / num_atoms    
                        gsoccupancy = np.count_nonzero(np.sum(lattice, axis=3)==0)
                        gs_occ.append(gsoccupancy/num_atoms)
                        # generates array of Booleans (True if it exists in ground state,
                        # False otherwise)
                        # and sums for ground state occupancy
                        
                        # calculation of total energy
                        energy.append(get_energy(lattice, mass))
                        #print(lattice)
                #print(energy)
                #fig, ax = plt.subplots(2)
                #ax[0].plot(energy)
                #ax[1].plot(gs_occ)
                #plt.show()
                #exit(1)

                av_en[idx][t] = np.average(energy[4000:])
                av_gs_occ[idx][t] = np.average(gs_occ[4000:])

        fig, ax = plt.subplots()
        av_en = np.average(av_en, axis=1)
        av_gs_occ = np.average(av_gs_occ, axis=1)
        print(1/betas)
        #ax[0].plot(1/betas, av_en)
        #ax[0].set(title="Average Energy", xlabel="Temperature", ylabel="Energy")
        ax.plot(1/betas, av_gs_occ)
        ax.set(title="Ground State Occupancy for {} Particles".format(num_atoms),
               xlabel="Temperature", ylabel="Occupancy")
        plt.savefig("./graphs/{}.png".format(num_atoms))
