import numpy as np
import scipy
from matplotlib import pyplot as plt

def initialize_lattice(num_atoms, beta, lim):
    ''' parameters:
            num_atoms: number of atoms we're simulating with. Should be cubic number.
            beta: temperature parameter, hbar * omega / T
            lim: highest energy state occupied
        returns:
               lattice: (num_atoms**1/3, num_atoms**1/3, num_atoms**1/3, 3) array with energy levels distributed according to BE distribution
      NOTE: Assumes energy levels of each particle same as that of quantum harmonic oscillator hbar * omega (n + 1/2)         
      '''
    
    nx = int(round(num_atoms ** 1/3))
    lattice = np.zeros([nx, nx, nx, 3], dtype=float)
    lat = np.random.rand(nx, nx, nx, 3)                            # generates (nx, nx, nx, 3) array of random integers
    prob = []
    for i in range(lim):
        compare = 1 / (np.exp(beta * (i+0.5))-1)                   # calculates probability of being in energy level i
        prob.append(compare)
    prob = np.cumsum(prob/np.sum(prob))                            # calculates cumulative probability for energy levels
    for i in range(lim - 1):                                       # places lattice sites in higher energy levels given probabilities
        sel = (lat < prob[i+1]) & (lat > prob[i])
        lattice[sel] = i
    return lattice                                                 # returns (nx, nx, nx, 3) array of energy levels 


def energy(lattice):
   ''' parameters:
   '''
    


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
    old_en = np.sum(np.square(lattice[site]))
    new_en = np.sum(np.square(lattice[site]+e_coord))
    return old_en - new_en


if __name__ == '__main__':

    # Initial conditions in reduced units
    beta = 0.25
    mass = 48.
    lim = 10
    num_atoms = 64
    nsweep = 100
    lat_dat = 'lat.dat'
    scalar_dat = 'scalar.dat'
    nconf = 10
    ntherm = 1
    therm_fmt = '%8d %14.8f %14.8f %14.8f\n'
    therm_header = '# isweep  m2  mavg  pe\n'

    # Initialize the lattice at specified temperature
    lattice = initialize_lattice(num_atoms, beta, lim)
    nx = int(num_atoms**1/3)
    restart = False
    if restart:
        lat = np.loadtxt('last_lat.dat')
        nx = lat.shape[0]
    else:
        with open(lat_dat, 'w') as flat:
            flat.write('# nx=%d\n' % nx)
        with open(scalar_dat, 'w') as fsca:
            fsca.write(therm_header)
    flat = open(lat_dat, 'a')
    fsca = open(scalar_dat, 'a')

    msg = 'starting %s sweeps using temperature %s and number of particles %s' % (nsweep, beta, num_atoms)
    if restart:
        msg = 're'+msg
    #etot = np.sum(energy(lattice))   #needs more work
    from time import clock
    start = clock()
    print(msg)
    print(np.shape(lattice))
    naccept = 0
    nattempt = 0
    for isweep in range(nsweep):
        if (isweep % ntherm == 0):
            #energy(lattice) ** 2        # calculation of specific heat
            # calculation of ground state occupancy
            # calculation of total energy
            pass
