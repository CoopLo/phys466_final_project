import numpy as np
from time import clock
import scipy
from scipy import integrate as integrate
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


if __name__ == '__main__':

    # Initial conditions in reduced units
    #betas = np.array([0.01, 0.05, 0.1, 0.5, 1.0, 5.0, 10.0, 1000.0])
    betas = np.array([1000., 100., 50., 40., 30., 20., 10.,
                      5., 4., 3., 2., 1., 0.5, 0.25, 0.1])
    betas = np.array([10.0])
    mass = 48.
    lim = 10
    nsweeps = [10000, 20000, 20000, 50000, 100000]#, 100000, 3000000]
    eq_frac = [0.4, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5]
    nsweeps = [500000]
    lat_dat = 'lat.dat'
    scalar_dat = 'scalar.dat'
    nconf = 10
    ntherm = 1
    therm_fmt = '%8d %14.8f %14.8f %14.8f\n'
    therm_header = '# isweep  m2  mavg  pe\n'

    runs = 10

    sizes = [4, 3, 4, 5, 7]#, 10, 22]
    #sizes = [22]
    lattice_snapshots = []
    for i, size in enumerate(sizes):
        num_atoms = size**3
        av_en = np.zeros((len(betas), runs))
        av_gs_occ = np.zeros((len(betas), runs))
        for t in range(runs):
            for idx, beta in enumerate(betas):
                nsweep = nsweeps[i]
                energy = []
                gs_occ = []
                lattice = be_lattice(size, beta)
                msg = 'starting %s sweeps using temperature %s and number of particles %s' % \
                  (nsweep, beta, num_atoms)
                print(msg)
                print(nsweeps[i])
                acc = 0
                for isweep in range(nsweeps[i]):
                    #if (isweep % ntherm == 0):

                    # Select site to change, coordinate, and up or down
                    site = np.random.randint(0, size, (1, 3))[0]
                    coord = np.random.randint(0, 3)
                    e_coord = [1,0,0] if coord==0 else [0,1,0] if coord==1 else [0,0,1] 
                    change = np.random.randint(0, 2)
                    change = -1 if change==0 else 1

                    # change in energy from site change
                    delta_e = energy_diff(lattice, site, np.array(e_coord)*change)
                    #print(delta_e)
                    
                    # transition probability
                    prob = 1/2 * 1./np.cosh(delta_e*beta/2) * np.exp(beta/2)
                    if(prob > np.random.random()):
                        #print("TRANSITIONING")
                        #print(prob)
                        lattice[site[0]][site[1]][site[2]][coord] += change
                        acc += 1

                    # calculation of ground state occupancy
                    #gsoccupancy = sum(np.shape(~lattice.any(axis=3))) / num_atoms    
                    gsoccupancy = np.count_nonzero(np.sum(lattice, axis=3)==0)
                    gs_occ.append(gsoccupancy/num_atoms)
                    # generates array of Booleans (True if it exists in ground state,
                    # False otherwise)
                    # and sums for ground state occupancy
                    
                    # calculation of total energy
                    energy.append(get_energy(lattice, mass))
                    lattice_snapshots.append(np.average(lattice, axis=3))
                    #print(lattice)
                #print(energy)
                fig, ax = plt.subplots(2)
                ax[0].plot(energy)
                ax[1].plot(gs_occ)
                with open('./data/lattice_{}_{}.txt'.format(beta, size), 'w') as fout:
                    fout.write(str(lattice_snapshots))
                #plt.show()
                print("ACCEPTANCE RATIO: {}".format(acc/nsweeps[i]))
                exit(1)

                av_en[idx][t] = np.average(energy[int(eq_frac[i]*nsweeps[i]):])
                av_gs_occ[idx][t] = np.average(gs_occ[int(eq_frac[i]*nsweeps[i]):])

        fig, ax = plt.subplots()
        av_en = np.average(av_en, axis=1)
        av_gs_occ = np.average(av_gs_occ, axis=1)
        print(1/betas)
        #ax[0].plot(1/betas, av_en)
        #ax[0].set(title="Average Energy", xlabel="Temperature", ylabel="Energy")
        #ax.plot(1/betas, av_gs_occ)
        #ax.set(title="Ground State Occupancy for {} Particles".format(num_atoms),
        #       xlabel="Temperature", ylabel="Occupancy")
        #plt.savefig("./graphs/{}.png".format(num_atoms))
